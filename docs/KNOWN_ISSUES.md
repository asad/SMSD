# Known issues

These are the known limitations of SMSD Pro 8.2.1, each with a workaround where
one exists. Each one has a GitHub issue you can follow; when a fix ships, the
[changelog](../CHANGELOG.md) records it and the issue is closed.

## SMARTS

**Directional bonds inside recursive SMARTS in Java.** ([issue 28](https://github.com/asad/SMSD/issues/28)) In Java, `/` and `\`
inside a recursive SMARTS such as `[$(c/c)]` match only single, non-aromatic
bonds, and `?` and negated directions are not supported there. Outside
`$(...)`, and everywhere in C++ and Python, they also match aromatic bonds.
Inside `$(...)` in Java, write a bond that may be aromatic without a direction,
such as `cc` or `c-,:c`.

## Drawings

**E/Z on aromatic rings.** ([issue 30](https://github.com/asad/SMSD/issues/30)) Stored E/Z on aromatic ring bonds, as in some
annulenes, may not be shown in generated coordinates and drawings.

**Some large macrocycles.** ([issue 31](https://github.com/asad/SMSD/issues/31)) For a few large macrocycles with many stored ring
double-bond configurations, such as abamectin, no layout is found that shows
every stored E/Z. Coordinate generation and drawing then raise an error
(`ValueError` in Python, `std::invalid_argument` in C++) rather than show the
wrong isomer.
