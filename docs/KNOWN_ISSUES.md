# Known issues

These are the known limitations of SMSD Pro 8.2.1, each with a workaround where
one exists. Each one has a GitHub issue you can follow; when a fix ships, the
[changelog](../CHANGELOG.md) records it and the issue is closed.

## SMARTS

**Directional bonds on aromatic rings in Java.** ([issue 28](https://github.com/asad/SMSD/issues/28)) In Java, `/` and `\` match
only single, non-aromatic bonds, so `c/c` does not match benzene. C++ and
Python also match aromatic bonds. In Java, write a bond that may be aromatic
without a direction, such as `cc` or `c-,:c`.

**Stereo double bonds on neighbouring aromatic ring atoms.** ([issue 29](https://github.com/asad/SMSD/issues/29)) A SMARTS taken from
a molecule with stereo double bonds on two neighbouring atoms of an aromatic
ring does not match that molecule. For example, `C/C=C/c1ccccc1/C=C/C` does not
match itself. Leave the directions off one of the two double bonds, as in
`C/C=C/c1ccccc1C=CC`; that double bond is then matched without its E/Z.

## Drawings

**E/Z on aromatic rings.** ([issue 30](https://github.com/asad/SMSD/issues/30)) Stored E/Z on aromatic ring bonds, as in some
annulenes, may not be shown in generated coordinates and drawings.

**Some large macrocycles.** ([issue 31](https://github.com/asad/SMSD/issues/31)) For a few large macrocycles with many stored ring
double-bond configurations, such as abamectin, no layout is found that shows
every stored E/Z. Coordinate generation and drawing then raise an error
(`ValueError` in Python, `std::invalid_argument` in C++) rather than show the
wrong isomer.
