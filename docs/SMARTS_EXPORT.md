# SMARTS export

`SmartsWriter` and `SearchResult.Match.toSmarts` are new in 8.1.0. Maven Central
stays at 8.0.3 in this release, so compile this example against a
[Java source build](HOWTO-INSTALL.md#java-source-build).

Java can write a portable SMARTS query from an MCS match or a whole molecular graph:

```java
import com.bioinception.smsd.core.MCS;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchResult;
import com.bioinception.smsd.core.SmartsWriter;
import com.bioinception.smsd.core.Substructure;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class SmartsExportExample {
    public static void main(String[] args) throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        MolGraph query = new MolGraph(parser.parseSmiles("CCO"));
        MolGraph target = new MolGraph(parser.parseSmiles("CCCO"));
        MolGraph libraryMolecule = new MolGraph(parser.parseSmiles("CCOCC"));
        MCS.Options options = new MCS.Options();
        SearchResult common = MCS.find(query).withOptions(options).match(target);
        if (!common.found()) throw new IllegalStateException("No common substructure");
        String pattern = common.matches().get(0).toSmarts(query, options.chemistry);
        SearchResult hit = Substructure.smarts(pattern).match(libraryMolecule);
        if (!hit.found()) throw new IllegalStateException("Exported pattern did not match");
        System.out.println(pattern);
        System.out.println(SmartsWriter.toSmarts(query, options.chemistry));
    }
}
```

`SmartsWriter.toSmarts(graph, chemistry)` writes every atom and bond. A match writes
only its selected atoms and common bonds, including an MCES match that leaves out
a source bond. Ring and aromaticity flags retain their original graph context.
The atom order follows the source indices; the text is deterministic for that
numbering and is not canonical.

The writer supports element wildcards, formal charge including neutral atoms,
isotope wildcards, strict or flexible aromaticity, bond-order modes and ring
membership. As in molecule matching, isotope zero matches any isotope on the
other side. The pattern is a subgraph query; it does not encode the induced MCIS
objective, a search budget or an optimality proof.

Export throws `IllegalArgumentException` for stereo, tautomer-aware matching,
complete-ring matching and strict ring-fusion matching. With flexible aromaticity
and strict bond order, it also rejects a non-aromatic single or double bond
between two aromatic atoms: the conditional resonance rule has no direct portable
bond expression. Use supported options explicitly rather than silently dropping
these constraints.

Portable output assumes molecular targets with consistent aromatic bond flags
and bond orders. Contradictory source annotations supplied through a graph builder
are rejected. Java SMARTS target preparation applies Daylight aromaticity. Strict
export therefore also rejects selected aromatic flags that change under that model,
including an uppercase Kekule benzene graph whose original flags are non-aromatic.
This check uses the full source graph, preserving the context of selected atoms
and bonds. Strict predicate parity with molecule matching requires library targets
whose aromatic flags also agree with that preparation. At most 99 Java ring closures
can be open at once; closure numbers are
reused after they close.
