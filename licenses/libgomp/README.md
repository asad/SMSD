# GCC OpenMP runtime

Linux wheel repair can bundle GCC's `libgomp` runtime. Its licence is GNU GPL
version 3 or later with the GCC Runtime Library Exception, version 3.1.
The unmodified licence documents are included here and in wheel metadata.

Linux wheels bundle the `libgomp` runtime package of the manylinux_2_28
(AlmaLinux 8) image they are built on. Its source is the `gcc` source package
of the same version in the [AlmaLinux vault](https://vault.almalinux.org/); for
example, `libgomp-8.5.0-28.el8_10.alma.1.x86_64` is built from
[`gcc-8.5.0-28.el8_10.alma.1.src.rpm`](https://vault.almalinux.org/8.10/BaseOS/Source/Packages/gcc-8.5.0-28.el8_10.alma.1.src.rpm).

Upstream sources and licence notices:

- https://gcc.gnu.org/projects/gomp/
- https://github.com/gcc-mirror/gcc/blob/master/libgomp/libgomp.h
- https://github.com/gcc-mirror/gcc/blob/master/COPYING3
- https://github.com/gcc-mirror/gcc/blob/master/COPYING.RUNTIME

SMSD's own source remains licensed under Apache-2.0.
