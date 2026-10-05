# Microsoft C++ and OpenMP runtimes

Windows wheels bundle Microsoft `msvcp140.dll` and `vcomp140.dll`, renamed by
wheel repair to avoid DLL name collisions. These components retain their
embedded Microsoft copyright notices and are governed by Microsoft's runtime
license terms. SMSD's Apache-2.0 license applies to SMSD's own code.

The unmodified Microsoft runtime license documents are included here:

- [Visual C++ Runtime 2015–2022](LICENSE-2022.docx), from
  [Microsoft's license page](https://visualstudio.microsoft.com/license-terms/vs2022-cruntime/).
- [Visual C++ V14 Redistributable and Runtime 2026](LICENSE-2026.docx), from
  [Microsoft's license page](https://visualstudio.microsoft.com/license-terms/vs2026-ga-visualcpp-v14-redist-runtime/).

The applicable terms follow the bundled runtime version. Build logs record
the selected DLL versions and source directories. Runtime selection requires
the same major version and at least the compiler toolset's minor version;
see [Microsoft's redistribution guidance](https://learn.microsoft.com/en-us/cpp/windows/determining-which-dlls-to-redistribute?view=msvc-170).

These notices concern the Microsoft components in Windows wheels. Linux and
macOS wheels use their separately licensed OpenMP runtimes.
