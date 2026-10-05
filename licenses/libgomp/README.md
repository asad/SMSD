# GCC OpenMP runtime

Linux wheel repair can bundle GCC's `libgomp` runtime. Its license is GNU GPL
version 3 or later with the GCC Runtime Library Exception, version 3.1.
The unmodified license documents are included here and in wheel metadata.

The validated Linux 7.2.0 wheel bundles the AlmaLinux runtime package
`libgomp-8.5.0-28.el8_10.alma.1.x86_64`. Its matching source package is
[`gcc-8.5.0-28.el8_10.alma.1.src.rpm`](https://vault.almalinux.org/8.10/BaseOS/Source/Packages/gcc-8.5.0-28.el8_10.alma.1.src.rpm).

Upstream sources and license notices:

- https://gcc.gnu.org/projects/gomp/
- https://github.com/gcc-mirror/gcc/blob/master/libgomp/libgomp.h
- https://github.com/gcc-mirror/gcc/blob/master/COPYING3
- https://github.com/gcc-mirror/gcc/blob/master/COPYING.RUNTIME

SMSD's own source remains licensed under Apache-2.0.
