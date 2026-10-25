# SPDX-License-Identifier: Apache-2.0
"""Offline archive fixtures; no Docker client or service is invoked."""
import hashlib
import importlib.util
import io
import json
import lzma
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import types
import unittest
from unittest import mock


SCRIPTS = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("image_archive", SCRIPTS / "archive_docker_image.py")
ARCHIVE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ARCHIVE)
REVISION = "a" * 40
IDENTITY = {"id": "sha256:" + "b" * 64, "os": "linux", "architecture": "amd64", "revision": REVISION}
PAYLOAD = b"retained image layer fixture\n" * 100


class ExportProcess:
    def __init__(self):
        self.stdout = io.BytesIO()
        self.returncode = None
        self.terminated = self.killed = False

    def poll(self):
        return self.returncode

    def terminate(self):
        self.terminated = True
        self.returncode = -15

    def kill(self):
        self.killed = True
        self.returncode = -9

    def wait(self, timeout):
        return self.returncode


class ImageArchiveTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.args = types.SimpleNamespace(image="smsd:release-check", revision=REVISION,
                                         architecture="amd64", output=self.root / "image.tar.xz", timeout=60)

    def run_export(self, chunks=None, identities=None, limit=ARCHIVE.ASSET_LIMIT):
        process = ExportProcess()
        def export(child, deadline, interrupted):
            self.assertIs(child, process)
            self.assertFalse(interrupted())
            yield from (chunks if chunks is not None else [PAYLOAD])
            child.returncode = 0
        with mock.patch.object(ARCHIVE, "inspect_image", side_effect=identities or [IDENTITY, IDENTITY]), \
                mock.patch.object(ARCHIVE.subprocess, "Popen", return_value=process) as spawn, \
                mock.patch.object(ARCHIVE, "export_chunks", side_effect=export), \
                mock.patch.object(ARCHIVE, "ASSET_LIMIT", limit):
            receipt = ARCHIVE.archive_image(self.args)
        self.assertEqual(spawn.call_args.args[0], ["docker", "image", "save", self.args.image])
        self.assertFalse(process.killed)
        return receipt, process

    def test_round_trip_hashes_platform_and_unmodified_image_identity(self):
        previous = signal.getsignal(signal.SIGTERM)
        receipt, process = self.run_export(chunks=[PAYLOAD[:31], PAYLOAD[31:]])
        self.assertEqual(lzma.decompress(self.args.output.read_bytes()), PAYLOAD)
        self.assertTrue(receipt["accepted"])
        self.assertTrue(receipt["archive_complete"])
        self.assertFalse(receipt["runtime_validation_performed"])
        self.assertEqual(receipt["archive_sha256"], hashlib.sha256(self.args.output.read_bytes()).hexdigest())
        self.assertEqual(receipt["input_sha256"], hashlib.sha256(PAYLOAD).hexdigest())
        self.assertEqual(receipt["input_bytes"], len(PAYLOAD))
        self.assertEqual(receipt["archive_bytes"], self.args.output.stat().st_size)
        self.assertEqual(receipt["export_returncode"], 0)
        self.assertFalse(process.terminated)
        self.assertEqual(signal.getsignal(signal.SIGTERM), previous)
        self.assertEqual(json.loads(self.args.output.with_name("image.tar.xz.json").read_text()), receipt)
        self.assertFalse(json.loads(self.args.output.with_name("image.tar.xz.intent.json").read_text())["accepted"])

    def test_asset_limit_is_strict_and_failure_retains_partial_archive_and_receipt(self):
        length = len(lzma.compress(PAYLOAD, format=lzma.FORMAT_XZ, check=lzma.CHECK_CRC64, preset=6))
        output, progress = io.BytesIO(), {}
        ARCHIVE.compress_chunks([PAYLOAD], output, progress, length + 1)
        self.assertEqual(len(output.getvalue()), length)
        receipt, _ = self.run_export(limit=length)
        self.assertFalse(receipt["accepted"])
        self.assertFalse(receipt["archive_complete"])
        self.assertIn("asset limit", receipt["error"])
        self.assertEqual(receipt["archive_bytes"], length - 1)
        self.assertEqual(self.args.output.stat().st_size, length - 1)
        self.assertTrue(self.args.output.with_name("image.tar.xz.stderr.log").exists())
        self.assertTrue(self.args.output.with_name("image.tar.xz.json").exists())
        with self.assertRaises(lzma.LZMAError):
            lzma.decompress(self.args.output.read_bytes())

    def test_changed_image_during_export_retains_full_archive_without_accepting_it(self):
        receipt, _ = self.run_export(identities=[IDENTITY, {**IDENTITY, "id": "sha256:" + "c" * 64}])
        self.assertFalse(receipt["accepted"])
        self.assertIn("identity changed", receipt["error"])
        self.assertEqual(lzma.decompress(self.args.output.read_bytes()), PAYLOAD)

    def test_failed_or_stalled_export_retains_prefix_and_closes_only_its_child(self):
        for error in (ValueError("docker image save failed"), TimeoutError("deadline")):
            with self.subTest(error=type(error).__name__):
                self.args.output = self.root / (type(error).__name__ + ".tar.xz")
                process = ExportProcess()
                def failed_export(*args):
                    yield PAYLOAD
                    raise error
                with mock.patch.object(ARCHIVE, "inspect_image", return_value=IDENTITY), \
                        mock.patch.object(ARCHIVE.subprocess, "Popen", return_value=process), \
                        mock.patch.object(ARCHIVE, "export_chunks", side_effect=failed_export):
                    receipt = ARCHIVE.archive_image(self.args)
                self.assertFalse(receipt["accepted"])
                self.assertIn(str(error), receipt["error"])
                self.assertTrue(process.terminated)
                self.assertTrue(process.stdout.closed)
                self.assertGreater(self.args.output.stat().st_size, 0)
                self.assertEqual(receipt["archive_sha256"], hashlib.sha256(self.args.output.read_bytes()).hexdigest())

    def test_cleanup_escalation_is_bounded(self):
        process = ExportProcess()
        process.wait = mock.Mock(side_effect=[subprocess.TimeoutExpired("docker", 10), 0])
        ARCHIVE.close_export(process)
        self.assertTrue(process.terminated)
        self.assertTrue(process.killed)
        self.assertEqual(process.wait.call_args_list, [mock.call(timeout=10), mock.call(timeout=10)])

    def test_unresolved_cleanup_cannot_be_accepted(self):
        with mock.patch.object(ARCHIVE, "close_export", side_effect=TimeoutError("unresolved exporter")):
            receipt, _ = self.run_export()
        self.assertFalse(receipt["accepted"])
        self.assertFalse(receipt["archive_complete"])
        self.assertIn("unresolved exporter", receipt["cleanup_error"])

    def test_empty_export_is_not_a_valid_image_archive(self):
        receipt, _ = self.run_export(chunks=[])
        self.assertFalse(receipt["accepted"])
        self.assertIn("no data", receipt["error"])

    def test_export_stream_reads_pipe_to_eof_and_checks_child_status(self):
        for status in (0, 17):
            with self.subTest(status=status):
                reader, writer = os.pipe()
                os.write(writer, PAYLOAD)
                os.close(writer)
                with os.fdopen(reader, "rb") as stream:
                    process = ExportProcess()
                    process.stdout = stream
                    process.returncode = status
                    if status == 0:
                        self.assertEqual(b"".join(ARCHIVE.export_chunks(process, float("inf"), lambda: False)), PAYLOAD)
                    else:
                        with self.assertRaisesRegex(ValueError, "save failed"):
                            list(ARCHIVE.export_chunks(process, float("inf"), lambda: False))

    def test_interrupted_export_does_not_wait_for_pipe_data(self):
        reader, writer = os.pipe()
        try:
            with os.fdopen(reader, "rb") as stream:
                process = ExportProcess()
                process.stdout = stream
                with self.assertRaises(InterruptedError):
                    list(ARCHIVE.export_chunks(process, float("inf"), lambda: True))
        finally:
            os.close(writer)

    def test_existing_output_or_evidence_is_never_overwritten(self):
        for suffix in ("", ".json", ".intent.json", ".stderr.log"):
            with self.subTest(suffix=suffix), mock.patch.object(ARCHIVE.subprocess, "Popen") as spawn:
                path = self.args.output.with_name(self.args.output.name + suffix)
                path.write_bytes(b"retained sentinel")
                with self.assertRaisesRegex(ValueError, "already exists"):
                    ARCHIVE.archive_image(self.args)
                self.assertEqual(path.read_bytes(), b"retained sentinel")
                spawn.assert_not_called()
                path.unlink()

    def test_bad_arguments_fail_before_any_docker_call(self):
        for field, value in (("image", "--unsafe"), ("revision", "HEAD"), ("architecture", "aarch64"),
                             ("timeout", 0), ("timeout", 3601), ("timeout", True), ("output", self.root / "image.gz")):
            with self.subTest(field=field, value=value), mock.patch.object(ARCHIVE.subprocess, "run") as run, \
                    mock.patch.object(ARCHIVE.subprocess, "Popen") as spawn:
                changed = types.SimpleNamespace(**{**vars(self.args), field: value})
                with self.assertRaises(ValueError):
                    ARCHIVE.archive_image(changed)
                run.assert_not_called()
                spawn.assert_not_called()

    def test_inspection_projects_only_identity_and_rejects_wrong_source_or_platform(self):
        for identity, accepted in ((IDENTITY, True), ({**IDENTITY, "revision": "c" * 40}, False),
                                   ({**IDENTITY, "os": "windows"}, False), ({**IDENTITY, "architecture": "arm64"}, False)):
            with self.subTest(identity=identity), mock.patch.object(ARCHIVE.subprocess, "run",
                    return_value=types.SimpleNamespace(stdout=json.dumps(identity))) as run:
                if accepted:
                    self.assertEqual(ARCHIVE.inspect_image(self.args.image, REVISION, "amd64", float("inf")), IDENTITY)
                else:
                    with self.assertRaisesRegex(ValueError, "platform or frozen revision"):
                        ARCHIVE.inspect_image(self.args.image, REVISION, "amd64", float("inf"))
                self.assertEqual(run.call_args.args[0], ["docker", "image", "inspect", "--format", ARCHIVE.INSPECT_FORMAT, self.args.image])
                self.assertNotIn(".Config.Env", ARCHIVE.INSPECT_FORMAT)


if __name__ == "__main__":
    unittest.main()
