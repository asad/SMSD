# SPDX-License-Identifier: Apache-2.0
"""The Metal accelerator loads, matches the installed smsd and, on an Apple GPU, attaches."""

from importlib.metadata import version

import smsd
import smsd_metal

assert version("smsd_metal") == smsd.__version__
reason = smsd_metal.unavailable_reason()
if reason is None:
    assert smsd_metal.is_available() and smsd_metal.device_info(), smsd.gpu_device_info()
    print("Metal device:", smsd_metal.device_info())
else:
    print("No usable Apple GPU here:", reason)
