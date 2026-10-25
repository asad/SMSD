#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Save a local Linux image as xz, retaining failures and enforcing the asset size limit."""

import argparse
import hashlib
import json
import lzma
import os
from pathlib import Path
import re
import selectors
import signal
import subprocess
import sys
import time


ASSET_LIMIT = 2_000_000_000
INSPECT_FORMAT = ('{"id":{{json .Id}},"os":{{json .Os}},"architecture":{{json .Architecture}},'
                  '"revision":{{json (index .Config.Labels "org.opencontainers.image.revision")}}}')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def remaining(deadline):
    seconds = deadline - time.monotonic()
    if seconds <= 0:
        raise TimeoutError("Image archive deadline exceeded")
    return seconds


def inspect_image(image, revision, architecture, deadline):
    result = subprocess.run(["docker", "image", "inspect", "--format", INSPECT_FORMAT, image],
                            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=min(30, remaining(deadline)), text=True)
    identity = json.loads(result.stdout)
    require(set(identity) == {"id", "os", "architecture", "revision"}
            and re.fullmatch(r"sha256:[0-9a-f]{64}", identity.get("id", ""))
            and identity["os"] == "linux" and identity["architecture"] == architecture
            and identity["revision"] == revision, "Image platform or frozen revision differs")
    return identity


def export_chunks(process, deadline, interrupted):
    with selectors.DefaultSelector() as observer:
        observer.register(process.stdout, selectors.EVENT_READ)
        while True:
            if interrupted():
                raise InterruptedError("Image export interrupted")
            if not observer.select(min(1, remaining(deadline))):
                continue
            chunk = os.read(process.stdout.fileno(), 1024 * 1024)
            if not chunk:
                break
            yield chunk
    require(process.wait(timeout=remaining(deadline)) == 0, "docker image save failed")


def compress_chunks(chunks, output, progress, limit=ASSET_LIMIT):
    compressor = lzma.LZMACompressor(format=lzma.FORMAT_XZ, check=lzma.CHECK_CRC64, preset=6)
    archive_hash, input_hash = hashlib.sha256(), hashlib.sha256()
    progress.update(archive_bytes=0, input_bytes=0)

    def emit(data):
        available = limit - 1 - progress["archive_bytes"]
        prefix = data[:max(0, available)]
        output.write(prefix)
        archive_hash.update(prefix)
        progress["archive_bytes"] += len(prefix)
        require(len(prefix) == len(data), "Compressed archive reaches the 2,000,000,000-byte asset limit")

    try:
        for chunk in chunks:
            input_hash.update(chunk)
            progress["input_bytes"] += len(chunk)
            emit(compressor.compress(chunk))
        require(progress["input_bytes"] > 0, "docker image save produced no data")
        emit(compressor.flush())
    finally:
        progress.update(archive_sha256=archive_hash.hexdigest(), input_sha256=input_hash.hexdigest())


def close_export(process):
    if process.poll() is None:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=10)
    process.stdout.close()


def archive_image(args):
    require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:/@-]*", args.image), "Expected a local image reference")
    require(re.fullmatch(r"[0-9a-f]{40}", args.revision), "Expected the full frozen source commit")
    require(args.output.name.endswith(".tar.xz"), "Archive output must end with .tar.xz")
    require(args.architecture in ("amd64", "arm64"), "Expected a supported Linux architecture")
    require(type(args.timeout) is int and 0 < args.timeout <= 3600, "Timeout must be 1 to 3600 seconds")
    paths = {name: args.output.with_name(args.output.name + suffix) for name, suffix in
             (("intent", ".intent.json"), ("record", ".json"), ("stderr", ".stderr.log"))}
    require(all(not path.exists() and not path.is_symlink() for path in [args.output, *paths.values()]),
            "Archive output or evidence already exists")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    record = {"schema_version": 1, "accepted": False, "archive_complete": False,
               "archive": args.output.name, "image": args.image, "source_commit": args.revision,
               "architecture": args.architecture, "maximum_asset_bytes_exclusive": ASSET_LIMIT,
               "timeout_seconds": args.timeout, "runtime_validation_performed": False,
               "stderr": paths["stderr"].name, "export_returncode": None, "error": None}
    with paths["intent"].open("x") as stream:
        json.dump(record, stream, indent=2)
        stream.write("\n")
    process, interrupted = None, []
    saved = {item: signal.signal(item, lambda number, frame: interrupted.append(number))
             for item in (signal.SIGTERM, signal.SIGINT)}
    deadline = time.monotonic() + args.timeout
    try:
        with args.output.open("xb") as output, paths["stderr"].open("xb") as log:
            record["image_identity"] = inspect_image(args.image, args.revision, args.architecture, deadline)
            require(not interrupted, "Image export interrupted before launch")
            process = subprocess.Popen(["docker", "image", "save", args.image],
                                       stdout=subprocess.PIPE, stderr=log)
            compress_chunks(export_chunks(process, deadline, lambda: bool(interrupted)), output, record, ASSET_LIMIT)
            require(not interrupted, "Image export interrupted")
            require(inspect_image(args.image, args.revision, args.architecture, deadline)
                    == record["image_identity"], "Image identity changed during export")
            record.update(accepted=True, archive_complete=True)
    except Exception as error:
        record.update(accepted=False, archive_complete=False, error=type(error).__name__ + ": " + str(error))
    finally:
        if process is not None:
            try:
                close_export(process)
            except Exception as error:
                record.update(accepted=False, archive_complete=False,
                               cleanup_error=type(error).__name__ + ": " + str(error))
            record["export_returncode"] = process.returncode
        for item, handler in saved.items():
            signal.signal(item, handler)
        if args.output.exists():
            digest = hashlib.sha256()
            with args.output.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            record.update(archive_bytes=args.output.stat().st_size, archive_sha256=digest.hexdigest())
        with paths["record"].open("x") as stream:
            json.dump(record, stream, indent=2)
            stream.write("\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--architecture", choices=("amd64", "arm64"), default="amd64")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=900, help="Export deadline, plus at most 20 seconds of cleanup")
    args = parser.parse_args()
    try:
        record = archive_image(args)
        print(json.dumps(record, sort_keys=True))
        return 0 if record["accepted"] else 1
    except (ValueError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
