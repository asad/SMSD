# C++ and Java

The C++ headers and the Java library offer the substructure, MCS and
fingerprint searches of the Python package, with the same default chemistry.
This chapter shows one complete program in each language, how to build and run
it, and where the languages differ. The [C++ guide](../CPP.md) and the
[Java guide](../JAVA.md) are the references for everything else.

## C++

The C++ library is header-only and needs a C++17 compiler: include
`smsd/smsd.hpp` and add `cpp/include` to the include path. An installed CMake
package provides the target `smsd::smsd` through
`find_package(smsd 8.1.0 CONFIG REQUIRED)`; the
[C++ module README](../../cpp/README.md#install-and-use-with-cmake) shows the
install and a minimal consumer project.

This program computes ECFP4 similarities of phenol to five molecules, screens a
prepared `TargetCorpus` at a threshold of 0.5, and runs exact substructure and
MCS searches over the same corpus:

```cpp
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "smsd/smsd.hpp"

int main() {
    const std::vector<std::string> names = {"phenol", "cyclohexanol", "aspirin", "ethanol", "toluene"};
    const std::vector<std::string> smiles = {"Oc1ccccc1", "OC1CCCCC1", "CC(=O)Oc1ccccc1C(=O)O", "CCO",
                                             "Cc1ccccc1"};
    const double expected[] = {1.000, 0.000, 0.020, 0.000, 0.040};
    const auto query = smsd::parseSMILES("Oc1ccccc1");                             // phenol
    const auto queryFp = smsd::batch::fingerprintFromSmiles("Oc1ccccc1", 2, 2048);  // ECFP4, 2,048 bits

    std::cout << std::fixed << std::setprecision(3);
    std::vector<smsd::MolGraph> targets;
    for (size_t i = 0; i < smiles.size(); ++i) {
        targets.push_back(smsd::parseSMILES(smiles[i]));
        const auto fp = smsd::batch::fingerprintFromSmiles(smiles[i], 2, 2048);
        const double similarity = smsd::batch::fingerprintTanimoto(queryFp, fp);
        if (std::abs(similarity - expected[i]) > 0.0005) return 1;
        std::cout << names[i] << " ECFP4 Tanimoto " << similarity << '\n';
    }

    smsd::batch::TargetCorpus corpus;
    corpus.addTargets(targets);
    corpus.prewarm(2);  // computes the fingerprints that screen() needs
    const std::vector<int> similar = corpus.screen(query, 0.5, 2);

    // The targets each corpus search reached; a search that stopped fails the program.
    bool stopped = false;
    auto hits = [&](const std::vector<smsd::CorpusSearchResult>& rows) {
        std::vector<int> found;
        for (const auto& row : rows) {
            if (row.thresholdReached) found.push_back(static_cast<int>(row.targetIndex));
            else if (!row.provenBelowThreshold) stopped = true;
        }
        return found;
    };
    const auto embeddings = hits(smsd::Substructure::find(query).matchCorpus(corpus, false, 2));
    const auto common = hits(smsd::MCS::find(query).withTimeout(1000).matchCorpus(corpus, 6, false, 2));
    if (stopped || similar != std::vector<int>{0} || embeddings != std::vector<int>{0, 1, 2}
        || common != std::vector<int>{0, 1, 2, 4}) return 1;

    auto print = [&](const char* label, const std::vector<int>& indices) {
        std::cout << label;
        for (int index : indices) std::cout << ' ' << names[index];
        std::cout << '\n';
    };
    print("screen at 0.5:", similar);
    print("substructure hits:", embeddings);
    print("MCS of at least 6 atoms:", common);
}
```

The program prints:

```text
phenol ECFP4 Tanimoto 1.000
cyclohexanol ECFP4 Tanimoto 0.000
aspirin ECFP4 Tanimoto 0.020
ethanol ECFP4 Tanimoto 0.000
toluene ECFP4 Tanimoto 0.040
screen at 0.5: phenol
substructure hits: phenol cyclohexanol aspirin
MCS of at least 6 atoms: phenol cyclohexanol aspirin toluene
```

Save it as `screen.cpp` and build it from the repository root:

```bash
c++ -std=c++17 -O2 -pthread -I cpp/include screen.cpp -o screen
./screen
```

The calls used are:

- `smsd::batch::fingerprintFromSmiles(smiles, radius = 2, fpSize = 2048)`: an
  ECFP fingerprint as 64-bit words, and `fingerprintTanimoto(a, b)` to compare
  two of them.
- `TargetCorpus::addTargets(graphs)` and `prewarm(numThreads = 0)`, which
  computes the ECFP4 fingerprints, radius 2 and 2,048 bits, that
  `screen(query, threshold, numThreads = 0)` compares. `screen` returns the
  targets whose similarity is at least the threshold, and nothing before
  `prewarm`.
- `Substructure::matchCorpus(corpus, gpuPrefilter = false, numThreads = 0)` and
  `MCS::matchCorpus(corpus, minimumAtoms, gpuPrefilter = false,
  numThreads = 0)`: one `CorpusSearchResult` per target, in corpus order.
  `thresholdReached` means a match, or a common substructure of at least
  `minimumAtoms` atoms, was found. `provenBelowThreshold` means the target is
  proved not to reach it, either by a completed search or because it was ruled
  out before matching. `result` holds the exact search result of every target
  that was searched. When neither flag is set, the search stopped.

Cyclohexanol is a substructure hit because the default chemistry lets aromatic
bonds match single and double bonds; see
[default chemistry](GETTING_STARTED.md#default-chemistry). Batch operations use
OpenMP when the compiler enables it; the
[installation guide](../HOWTO-INSTALL.md#c-source-build) covers the CMake
build and its options.

## Java

The Java library is `com.bioinceptionlabs:smsd` and the command-line tool
`com.bioinceptionlabs:smsd-cli`. Maven Central has 8.0.3; for 8.1.0, run
`mvn -f java/pom.xml install` from the repository root, which installs both in
the local Maven repository (Gradle then needs `mavenLocal()`). The
[Java module README](../../java/README.md#installation) gives the Maven and
Gradle dependencies, and the code runs on Java 8 or later.

This program computes the same ECFP4 similarities and runs exact substructure
and MCS searches over the list on two threads:

```java
import com.bioinception.smsd.core.FingerprintEngine;
import com.bioinception.smsd.core.MCS;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchResult;
import com.bioinception.smsd.core.Substructure;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/** ECFP4 similarity, then exact substructure and MCS searches, over a small list. */
public class ScreenExample {
    public static void main(String[] args) throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        String[] names = {"phenol", "cyclohexanol", "aspirin", "ethanol", "toluene"};
        String[] smiles = {"Oc1ccccc1", "OC1CCCCC1", "CC(=O)Oc1ccccc1C(=O)O", "CCO", "Cc1ccccc1"};
        double[] expectedSimilarity = {1.000, 0.000, 0.020, 0.000, 0.040};
        boolean[] expectedHit = {true, true, true, false, false};
        int[] expectedAtoms = {7, 7, 7, 3, 6};

        MolGraph query = FingerprintEngine.toMolGraph(parser.parseSmiles("Oc1ccccc1"));  // phenol
        long[] queryFp = FingerprintEngine.ecfp(query, 2, 2048);                           // ECFP4
        List<MolGraph> targets = new ArrayList<>();
        for (String text : smiles) targets.add(FingerprintEngine.toMolGraph(parser.parseSmiles(text)));

        List<SearchResult> embeddings = Substructure.find(query).matchMany(targets, 2);
        List<SearchResult> common = MCS.find(query).withTimeout(1000).matchMany(targets, 2);
        for (int i = 0; i < targets.size(); i++) {
            SearchResult embedding = embeddings.get(i);
            SearchResult mcs = common.get(i);
            double similarity = FingerprintEngine.tanimoto(queryFp, FingerprintEngine.ecfp(targets.get(i), 2, 2048));
            int atoms = mcs.found() ? mcs.matches().get(0).atomCount() : 0;
            if (!embedding.complete() || !mcs.optimal() || embedding.found() != expectedHit[i]
                || atoms != expectedAtoms[i] || Math.abs(similarity - expectedSimilarity[i]) > 0.0005)
                throw new IllegalStateException("Unexpected result for " + names[i]);
            System.out.println(String.format(Locale.ROOT, "%-13s ECFP4 %.3f  substructure %-5s  MCS atoms %d",
                names[i], similarity, embedding.found(), atoms));
        }
    }
}
```

The program prints:

```text
phenol        ECFP4 1.000  substructure true   MCS atoms 7
cyclohexanol  ECFP4 0.000  substructure true   MCS atoms 7
aspirin       ECFP4 0.020  substructure true   MCS atoms 7
ethanol       ECFP4 0.000  substructure false  MCS atoms 3
toluene       ECFP4 0.040  substructure false  MCS atoms 6
```

Save it as `ScreenExample.java` and compile it from the repository root
against the self-contained JAR of the source build; on Windows, separate the
class path entries with a semicolon:

```bash
javac --release 8 -cp java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar ScreenExample.java
java -cp "java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar:." ScreenExample
```

With JDK 25, `javac --release 8` warns that release 8 is obsolete; the class
still compiles and runs on Java 8 or later.

The Java searches take CDK molecules or `MolGraph` graphs, and
`FingerprintEngine.toMolGraph` standardises a CDK molecule into a graph. Java
has no prepared target corpus, so it screens a list with `matchMany` and
fingerprint comparisons, as above.

## Differences between the languages

| Topic | Python | C++ | Java |
| --- | --- | --- | --- |
| Thread count 0 | The OpenMP default, which follows `OMP_NUM_THREADS` | Searches: one worker per hardware thread; fingerprint work: the OpenMP default, or one thread in a build without OpenMP | Every available processor; `matchMany(targets)` without a count runs on one thread |
| MCIS search | The 8.1.0 search | The 8.1.0 search | The earlier search, so a stopped search can return a different result |
| Prepared target corpus and corpus searches | `TargetCorpus`, `match_corpus` | `TargetCorpus`, `matchCorpus` | No |
| Streamed SMILES screens | `stream_screen_smiles` | No | No |
| Packed fingerprint index files, `butina`, `maxmin` | Yes, with NumPy | No | No |
| NVIDIA CUDA screening | The `smsd-cuda` package | A CUDA source build | No |
| Metal fingerprint screening | No | macOS, when requested explicitly | No |
| Named SMARTS predicates (`$name`) | No | Yes | Yes |
| Default chemistry | The defaults of [default chemistry](GETTING_STARTED.md#default-chemistry) | The same | The same |

The [Java guide](../JAVA.md#results) describes the stopped-search difference,
and the [C++ guide](../CPP.md#optional-metal-fingerprint-screening) the Metal
provider.

## See also

- [C++ guide](../CPP.md) and [C++ module README](../../cpp/README.md): the
  full C++ API, CMake options and the installed package.
- [Java guide](../JAVA.md) and [Java module README](../../java/README.md): the
  full Java API, dependencies and the source build.
- [Command line and containers](COMMAND_LINE.md): the Java command-line tool.
- [Screening a library](SCREENING.md): the screening routes in Python.
