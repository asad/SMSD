# Command line and containers

The command-line tool `smsd` runs one substructure, SMARTS or maximum common
substructure (MCS) search between a query and a target, exports the atom
mappings, and searches every record of an SDF file, reporting the results as
JSON and through its exit status. This chapter covers the tool, the Docker
images that carry it, and the local molecular tools server. Fingerprints and
similarity screening are not available on the command line; the
[limitations](#limitations) section says where to find them.

## Get the command-line tool

The 8.1.0 tool comes from a Java source build or from the CPU
[Docker images](#docker-images). Maven Central and the 8.0.3 release carry the
8.0.3 tool, which the
[installation guide](../HOWTO-INSTALL.md#portable-java-cli) describes. To build
8.1.0, use a JDK (Java 25 LTS preferred) and Maven 3.9 or later, and from the
repository root run:

```bash
mvn -f java/pom.xml clean package
export PATH="$PWD/java/cli/src/scripts:$PATH"
smsd --version
```

The build runs the default tests and writes
`java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar`, which holds the
tool and all its dependencies and runs on Java 8 or later, also as
`java -jar smsd-cli-8.1.0-jar-with-dependencies.jar`. The launcher script
`java/cli/src/scripts/smsd` runs the newest such JAR in `java/cli/target`, or
one next to the script; `smsd.bat` and `smsd.ps1` in the same folder do the
same on Windows. The examples in this chapter call the tool as `smsd`.

| Variable | Effect |
| --- | --- |
| `SMSD_JAR` | The JAR to run, instead of the newest one found (`smsd` and `smsd.bat`) |
| `JAVA_HOME` | Run `$JAVA_HOME/bin/java` instead of the `java` on the path (`smsd` and `smsd.bat`) |
| `JAVA_OPTS` | Options for the Java virtual machine, such as `-Xmx2g` (all three launchers) |

## Synopsis

```text
smsd --Q TYPE --q QUERY --T TYPE --t TARGET [--mode sub|mcs] [options]
```

`--Q` and `--T` name the input types, and `--q` and `--t` give a SMILES or
SMARTS string or a file path:

| Type | `--q` or `--t` holds | Query | Target |
| --- | --- | --- | --- |
| `SMI` | One SMILES string | Yes | Yes |
| `SIG` | One SMARTS string | Yes, in substructure mode | No |
| `MOL` | The path of an MDL MOL file | Yes | Yes |
| `ML2` | The path of a Mol2 file | Yes | Yes |
| `PDB` | The path of a PDB file with one model | Yes | Yes |
| `CML` | The path of a CML file with one molecule | Yes | Yes |
| `SDF` | The path of an SDF file of up to 100,000 records | No | Yes, as a batch |

Type names and the mode ignore case. Each molecule is standardised first
(atom types, implicit hydrogens, rings and aromaticity); a molecule that cannot
be standardised is used as read. Without `--json` the tool writes no report
and answers only through its [exit status](#exit-status); `--json -` writes the
report to standard output. The first lines of `smsd --help` list every option:

```bash
smsd --help | head -5
```

The command prints:

```text
Usage: smsd [-ahmrV] [--benchmark] [--fragments] [--induced] [--json-pretty]
            [--map-pretty] [--json=<jsonOut>] [--map-format=<mapFormat>]
            [--map-out=<mapOut>] [--mode=<mode>] [-o=<mcsFragmentOut>]
            [--objective=<objective>] --q=<q> --Q=<qType> --t=<t> --T=<tType>
            [--threads=<numThreads>] [--timeout=<timeoutMs>]
```

## Substructure search

Substructure mode is the default. Every report carries two effort figures,
`nodes_visited` and `elapsed_ms`, that change from run to run, so the examples
in this chapter delete them with `sed`; leave the filter out to see them. This
search looks for benzene in phenol:

```bash
smsd --Q SMI --q "c1ccccc1" --T SMI --t "Oc1ccccc1" --json - \
  | sed -E 's/,"(nodes_visited|elapsed_ms)":[0-9]+//g'
```

The command prints the following, and `smsd` exits with status 0:

```json
{"query":"c1ccccc1","target":"Oc1ccccc1","mode":"sub","timeout_ms":10000,"exists":true,"status":"Completed","optimal":true,"complete":true,"score":6.0,"bonds":[{"query":[0,1],"target":[2,3]},{"query":[0,5],"target":[1,2]},{"query":[1,2],"target":[3,4]},{"query":[2,3],"target":[4,5]},{"query":[3,4],"target":[5,6]},{"query":[4,5],"target":[1,6]}],"fragments":[[0,1,2,3,4,5]]}
```

`exists` answers the question, and `status`, `optimal` and `complete` say
whether the answer is proved, as described under
[reading a result](GETTING_STARTED.md#reading-a-result). `bonds` pairs each
query bond with its target bond, by atom indices in input order.

A SMARTS query uses the type `SIG`. The first search finds a carbonyl group in
acetic acid and exits with status 0; the second finds none in ethanol and exits
with status 1:

```bash
smsd --Q SIG --q "C=O" --T SMI --t "CC(=O)O" --json - --json-pretty \
  | sed -E '/"(nodes_visited|elapsed_ms)"/d'
smsd --Q SIG --q "C=O" --T SMI --t "CCO" --json - \
  | sed -E 's/,"(nodes_visited|elapsed_ms)":[0-9]+//g'
```

The commands print:

```text
{
  "query" : "C=O",
  "target" : "CC(=O)O",
  "mode" : "sub",
  "timeout_ms" : 10000,
  "exists" : true,
  "status" : "Completed",
  "optimal" : true,
  "complete" : true,
  "score" : 2.0,
  "bonds" : [ {
    "query" : [ 0, 1 ],
    "target" : [ 1, 2 ]
  } ],
  "fragments" : [ [ 0, 1 ] ]
}
{"query":"C=O","target":"CCO","mode":"sub","timeout_ms":10000,"exists":false,"status":"Completed","optimal":true,"complete":true,"score":0.0}
```

## MCS

`--mode mcs` finds a connected maximum common induced subgraph (MCIS) by
default. For N-methylacetamide and propanamide:

```bash
smsd --Q SMI --q "CC(=O)NC" --T SMI --t "CCC(=O)N" --mode mcs \
  --json - --json-pretty \
  | sed -E '/"(nodes_visited|elapsed_ms)"/d'
```

The command prints the following, and `smsd` exits with status 0:

```json
{
  "query" : "CC(=O)NC",
  "target" : "CCC(=O)N",
  "mode" : "mcs",
  "timeout_ms" : 10000,
  "status" : "Completed",
  "optimal" : true,
  "complete" : true,
  "score" : 4.0,
  "bonds" : [ {
    "query" : [ 0, 1 ],
    "target" : [ 1, 2 ]
  }, {
    "query" : [ 1, 2 ],
    "target" : [ 2, 3 ]
  }, {
    "query" : [ 1, 3 ],
    "target" : [ 2, 4 ]
  } ],
  "fragments" : [ [ 0, 1, 2, 3 ] ],
  "objective" : "MCIS",
  "mcs_size" : 4,
  "pairs" : [ [ 0, 1 ], [ 1, 2 ], [ 2, 3 ], [ 3, 4 ] ],
  "query_atoms" : [ 0, 1, 2, 3 ],
  "target_atoms" : [ 1, 2, 3, 4 ],
  "tanimoto" : 0.6666666666666666,
  "similarity_upper_bound" : 1.0,
  "mcs_smiles" : "C(C)(N)=O"
}
```

`mcs_size` counts the common atoms and `pairs` lists them as query and target
indices. `tanimoto` is the MCS similarity 4 / (5 + 5 − 4), not a fingerprint
similarity. `--objective MCES` maximises common bonds instead, and
`--fragments` allows a disconnected result, here the two ethyl groups of
diethyl ether and diethylamine:

```bash
smsd --Q SMI --q "CCOCC" --T SMI --t "CCNCC" \
  --mode mcs --objective MCES --fragments --timeout 1000 --json - \
  | sed -E 's/,"(nodes_visited|elapsed_ms)":[0-9]+//g'
```

The command prints:

```json
{"query":"CCOCC","target":"CCNCC","mode":"mcs","timeout_ms":1000,"status":"Completed","optimal":true,"complete":true,"score":2.0,"bonds":[{"query":[0,1],"target":[0,1]},{"query":[3,4],"target":[3,4]}],"fragments":[[0,1],[3,4]],"objective":"MCES","mcs_size":4,"pairs":[[0,0],[1,1],[3,3],[4,4]],"query_atoms":[0,1,3,4],"target_atoms":[0,1,3,4],"tanimoto":0.3333333333333333,"similarity_upper_bound":1.0,"mcs_smiles":"CC.CC"}
```

For MCES, `score` and `tanimoto` count bonds, while `mcs_size` still counts
atoms. `-o` writes the common graph as SMILES, to a file or, with `-`, to
standard output, and writes nothing when there is no common graph:

```bash
smsd --Q SMI --q "CC(=O)NC" --T SMI --t "CCC(=O)N" --mode mcs -o -
```

The command prints:

```text
C(C)(N)=O
```

## SDF batch search

With `--T SDF`, the query is searched against every record of the file. To
follow the examples, write a small file with the Python package:

```bash
python - <<'EOF'
import smsd

library = {"phenol": "Oc1ccccc1", "cyclohexanol": "OC1CCCCC1",
           "aspirin": "CC(=O)Oc1ccccc1C(=O)O", "ethanol": "CCO",
           "toluene": "Cc1ccccc1"}
graphs = []
for name, smiles in library.items():
    graph = smsd.parse_smiles(smiles)
    graph.name = name
    graphs.append(graph)
smsd.write_sdf(graphs, "targets.sdf")
EOF
```

The report has one row per record under `results`, in file order. A row
carries the record's position as `target_index` but not its title, so read the
titles from the file when you need them. These commands run a substructure
batch and an MCS batch, each on two threads, and print one line per record:

```bash
smsd --Q SMI --q "c1ccccc1" --T SDF --t targets.sdf --threads 2 --json benzene.json
smsd --Q SMI --q "Oc1ccccc1" --T SDF --t targets.sdf --mode mcs --threads 2 --json phenol.json
python - <<'EOF'
import json
import smsd

names = [graph.name for graph in smsd.read_sdf("targets.sdf")]
with open("benzene.json", encoding="utf-8") as handle:
    embeddings = json.load(handle)["results"]
with open("phenol.json", encoding="utf-8") as handle:
    common = json.load(handle)["results"]
for row, mcs in zip(embeddings, common):
    print(names[row["target_index"]], row["exists"], mcs["mcs_size"], round(mcs["tanimoto"], 3))
EOF
```

The commands print:

```text
phenol True 7 1.0
cyclohexanol True 7 1.0
aspirin True 7 0.538
ethanol False 3 0.429
toluene True 6 0.75
```

Benzene is found in cyclohexanol because of the fixed default chemistry; see
[limitations](#limitations). An MCS row has `mcs_size`, `pairs` and `tanimoto`
but no `mcs_smiles`. `--threads 0`, the default, uses every processor, so set
`--threads` on a shared machine. The whole file is read into memory before the
search starts.

## Mapping export

In substructure mode, `-m` enumerates up to 10,000 embeddings and reports them
under `mappings` in place of `exists`. Methanol embeds in methanediol in two
ways:

```bash
smsd --Q SMI --q "CO" --T SMI --t "OCO" -m --json - --json-pretty \
  | sed -E '/"(nodes_visited|elapsed_ms)"/d'
```

The command prints:

```json
{
  "query" : "CO",
  "target" : "OCO",
  "mode" : "sub",
  "timeout_ms" : 10000,
  "status" : "Completed",
  "optimal" : true,
  "complete" : true,
  "score" : 2.0,
  "bonds" : [ {
    "query" : [ 0, 1 ],
    "target" : [ 0, 1 ]
  } ],
  "fragments" : [ [ 0, 1 ] ],
  "mappings" : [ {
    "index" : 0,
    "pairs" : [ [ 0, 1 ], [ 1, 0 ] ],
    "query_atoms" : [ 0, 1 ],
    "target_atoms" : [ 1, 0 ]
  }, {
    "index" : 1,
    "pairs" : [ [ 0, 1 ], [ 1, 2 ] ],
    "query_atoms" : [ 0, 1 ],
    "target_atoms" : [ 1, 2 ]
  } ],
  "tanimoto" : 0.6666666666666666
}
```

`--map-out` writes the mappings of a single-pair search in the format chosen
by `--map-format`. In substructure mode it needs `-m` and otherwise writes
nothing; in MCS mode it works on its own. The `smiles` format writes, for each
mapping, the query and the target as SMILES with matching atom-map numbers,
separated by a tab:

```bash
smsd --Q SMI --q "CO" --T SMI --t "OCO" -m --map-out - --map-format smiles
smsd --Q SMI --q "CC(=O)NC" --T SMI --t "CCC(=O)N" --mode mcs \
  --map-out - --map-format smiles
```

The commands print:

```text
[CH3:1][OH:2]	[CH2:1]([OH:2])O
[CH3:1][OH:2]	[CH2:1](O)[OH:2]
[C:2]([CH3:1])([NH:4]C)=[O:3]	[C:2]([CH2:1]C)([NH2:4])=[O:3]
```

The formats are:

- `json` (the default): an object with `query_size`, `target_size` and
  `mappings`, written without a final newline; `--map-pretty` indents it.
- `smiles` or `smi`: atom-mapped SMILES, as above. `smarts` writes the same
  text.
- `mdl`, `mol` or `sdf`: MOL blocks with generated two-dimensional
  coordinates; in MCS mode the common graph of the first match, and in
  substructure mode the matched target atoms of each mapping. The second header
  line of each block carries a time stamp to the minute, so the output changes
  between runs.

When `--map-out -` and `--json -` both write to standard output, the export
comes first; a `json` export has no final newline and runs straight into the
report, so send them to separate files.
An SDF batch ignores `-m`, `--map-out` and `-o`.

## Hydrogens and induced matching

`--induced` requires an induced embedding: query atoms that are not bonded may
not map to bonded target atoms. Propane embeds in cyclopropane as a path, but
not as an induced subgraph:

```bash
smsd --Q SMI --q "CCC" --T SMI --t "C1CC1" --json - \
  | sed -E 's/,"(nodes_visited|elapsed_ms)":[0-9]+//g'
smsd --Q SMI --q "CCC" --T SMI --t "C1CC1" --induced --json - \
  | sed -E 's/,"(nodes_visited|elapsed_ms)":[0-9]+//g'
```

The commands print the following, and `smsd` exits with status 0 and then 1:

```text
{"query":"CCC","target":"C1CC1","mode":"sub","timeout_ms":10000,"exists":true,"status":"Completed","optimal":true,"complete":true,"score":3.0,"bonds":[{"query":[0,1],"target":[0,1]},{"query":[1,2],"target":[0,2]}],"fragments":[[0,1,2]]}
{"query":"CCC","target":"C1CC1","mode":"sub","timeout_ms":10000,"exists":false,"status":"Completed","optimal":true,"complete":true,"score":0.0}
```

Explicit hydrogens are atoms of the graph. A query with an explicit hydrogen,
such as `[H]OC`, finds no match in methanol written without one, unless `-r`
removes the explicit hydrogens first. `-a` adds explicit hydrogens to both
molecules, so methane must then match a carbon with four hydrogens, which
methanol does not have:

```bash
smsd --Q SMI --q "[H]OC" --T SMI --t "CO" --json - \
  | sed -E 's/,"(nodes_visited|elapsed_ms)":[0-9]+//g'
smsd --Q SMI --q "[H]OC" --T SMI --t "CO" -r --json - \
  | sed -E 's/,"(nodes_visited|elapsed_ms)":[0-9]+//g'
smsd --Q SMI --q "C" --T SMI --t "CO" -a --json - \
  | sed -E 's/,"(nodes_visited|elapsed_ms)":[0-9]+//g'
```

The commands print:

```text
{"query":"[H]OC","target":"CO","mode":"sub","timeout_ms":10000,"exists":false,"status":"Completed","optimal":true,"complete":true,"score":0.0}
{"query":"[H]OC","target":"CO","mode":"sub","timeout_ms":10000,"exists":true,"status":"Completed","optimal":true,"complete":true,"score":2.0,"bonds":[{"query":[0,1],"target":[0,1]}],"fragments":[[0,1]]}
{"query":"C","target":"CO","mode":"sub","timeout_ms":10000,"exists":false,"status":"Completed","optimal":true,"complete":true,"score":0.0}
```

## Exit status

| Status | Single pair | SDF batch |
| --- | --- | --- |
| 0 | The search completed with a match | At least one match, and every search completed |
| 1 | The search completed without a match | No match, and every search completed |
| 2 | Invalid input, an unknown mode, an error, or a search stopped as `TimedOut`, `LimitReached` or `Cancelled` | Any target that failed or stopped |

A stopped search whose report says `"exists": false` proves nothing, which is
why it exits with status 2. `--timeout 0` always stops at once, and reaching
the limit of 10,000 mappings with `-m` also counts as stopped. The exit status
is enough for a script, even without `--json`:

```bash
smsd --Q SMI --q "c1ccccc1" --T SMI --t "Oc1ccccc1"; echo "exit status $?"
smsd --Q SMI --q "c1ccccc1" --T SMI --t "CCO"; echo "exit status $?"
smsd --Q SMI --q "c1ccccc1" --T SMI --t "Oc1ccccc1" --timeout 0; echo "exit status $?"
smsd --Q SMI --q "C1CC" --T SMI --t "CCC" 2> /dev/null; echo "exit status $?"
```

The commands print:

```text
exit status 0
exit status 1
exit status 2
exit status 2
```

The third search timed out and the fourth query is not valid SMILES.

## JSON report

A single-pair report has these keys:

| Key | When | Value |
| --- | --- | --- |
| `query`, `target` | Always | The `--q` and `--t` values as given |
| `mode` | Always | The mode as typed |
| `timeout_ms` | Always | The time limit for the pair |
| `status` | Always | `Completed`, `TimedOut`, `Cancelled` or `LimitReached` |
| `optimal`, `complete` | Always | Whether the answer is proved, and whether the search ran to its end |
| `score` | Always | Substructure: the query atoms when found, otherwise 0.0; MCIS: atoms; MCES: bonds |
| `nodes_visited`, `elapsed_ms` | Always | Effort figures that vary between runs |
| `bonds`, `fragments` | A match was found | The bond pairs and the connected query atom sets of the first match |
| `exists` | Substructure mode without `-m` | Whether an embedding was found |
| `mappings` | `-m` | Each embedding: `index`, `pairs`, `query_atoms` and `target_atoms` |
| `objective` | MCS mode | `MCIS` or `MCES` |
| `mcs_size`, `pairs`, `query_atoms`, `target_atoms` | MCS mode | The common atoms of the first match |
| `tanimoto` | MCS mode or `-m` | MCS similarity: common atoms / (query atoms + target atoms − common atoms), or the same count of bonds for MCES |
| `similarity_upper_bound` | MCS mode | An upper bound on that similarity |
| `mcs_smiles` | MCS mode with a match | The common graph as SMILES |

An SDF batch report has `query`, `target_file`, `target_count`,
`threads_used`, `mode` and `results`. Each row has `target_index`, the outcome
keys from `status` to `elapsed_ms`, `bonds` and `fragments` when it matched,
and then `exists` in substructure mode, or `mcs_size`, `pairs` and `tanimoto`
in MCS mode. A row whose target could not be searched has an `error` message
instead.

The `tanimoto` key compares atom or bond counts of the MCS. It is not a
fingerprint Tanimoto coefficient: for phenol and aspirin it is 0.538, whereas
the ECFP4 Tanimoto coefficient of the pair is 0.020.

## Errors

Usage errors, such as a missing option, an unknown mode, a lower-case
objective or a SMARTS query in MCS mode, print a message and the usage text on
standard error:

```bash
smsd --Q SMI --q "CC" --T SMI --t "CCC" --mode foo 2>&1 | head -1
smsd --Q SMI --q "CCO" --T SMI --t "CCN" --mode mcs --objective mces 2>&1 | head -1
smsd --Q SIG --q "C=O" --T SMI --t "CC=O" --mode mcs 2>&1 | head -1
```

The commands print:

```text
Unknown mode: foo
Invalid value for option '--objective': expected one of [MCIS, MCES] (case-sensitive) but was 'mces'
SMARTS queries require substructure mode
```

Other errors, such as SMILES that cannot be parsed, an unsupported file type,
an unknown `--map-format` or a missing file, print a Java stack trace on
standard error. In every case standard output stays empty and the exit status
is 2. The Mol2 reader also prints warning lines on standard error for valid
files, so redirect standard error when a program reads the report from
standard output.

## Options

| Option | Value | Default | Meaning |
| --- | --- | --- | --- |
| `--Q` | Type | Required | Query type: `SMI`, `SIG`, `MOL`, `ML2`, `PDB` or `CML` |
| `--q` | Text or path | Required | Query SMILES or SMARTS, or a file path |
| `--T` | Type | Required | Target type: `SMI`, `MOL`, `ML2`, `PDB`, `CML` or `SDF` |
| `--t` | Text or path | Required | Target SMILES, or a file path |
| `--mode` | `sub` or `mcs` | `sub` | Substructure or MCS search |
| `--objective` | `MCIS` or `MCES` | `MCIS` | MCS objective; upper case only |
| `--fragments` | Flag | Off | Allow a disconnected MCS |
| `--induced` | Flag | Off | Require an induced embedding (substructure mode) |
| `-m`, `--mappings` | Flag | Off | Enumerate up to 10,000 embeddings (substructure mode) |
| `--timeout` | Milliseconds | 10,000 | Time limit for each query and target pair, in both modes; `0` stops at once |
| `-a`, `--add-hydrogens` | Flag | Off | Add explicit hydrogens to both molecules first |
| `-r`, `--remove-hydrogens` | Flag | Off | Remove explicit hydrogens from both molecules first |
| `--json` | Path or `-` | None | Write the JSON report; without it nothing is printed |
| `--json-pretty` | Flag | Off | Indent the JSON report |
| `--map-out` | Path or `-` | None | Write the mapping export of a single pair |
| `--map-format` | `json`, `smiles`, `smi`, `smarts`, `mdl`, `mol` or `sdf` | `json` | Mapping export format |
| `--map-pretty` | Flag | Off | Indent the JSON mapping export |
| `-o`, `--output` | Path or `-` | None | Write the MCS as SMILES (single-pair MCS mode) |
| `--benchmark` | Flag | Off | Run the substructure and MCS searches of one pair repeatedly and print their average times, which depend on the machine |
| `--threads` | Integer | 0 | Worker threads for an SDF batch; `0` uses every processor |
| `-h`, `--help` | | | Print the usage and exit |
| `-V`, `--version` | | | Print the version and exit |

## Limitations

- There are no fingerprint, similarity or index commands. Use the Python, C++
  or Java library for those; see [screening a library](SCREENING.md) and
  [C++ and Java](OTHER_LANGUAGES.md).
- The chemistry is always the default one: aromatic bonds also match single and
  double bonds, and formal charges are not compared. No option changes it, so
  benzene is found in cyclohexanol, as below. The libraries can switch to
  stricter chemistry.
- `SMI` takes one SMILES string, not a file. A dotted SMILES is one
  disconnected molecule. To search many SMILES, convert them to an SDF file
  first.
- An SDF file may hold up to 100,000 records and is read into memory first.
  Report rows give `target_index` but not the record title.
- `PDB` and `CML` files must hold exactly one molecule or model.
- The time limit is 10,000 ms per pair in both modes, whereas the libraries
  give an MCS a limit that grows with the size of the molecules.

```bash
smsd --Q SMI --q "c1ccccc1" --T SMI --t "OC1CCCCC1" --json - \
  | sed -E 's/,"(nodes_visited|elapsed_ms)":[0-9]+//g'
```

The command prints:

```json
{"query":"c1ccccc1","target":"OC1CCCCC1","mode":"sub","timeout_ms":10000,"exists":true,"status":"Completed","optimal":true,"complete":true,"score":6.0,"bonds":[{"query":[0,1],"target":[2,3]},{"query":[0,5],"target":[1,2]},{"query":[1,2],"target":[3,4]},{"query":[2,3],"target":[4,5]},{"query":[3,4],"target":[5,6]},{"query":[4,5],"target":[1,6]}],"fragments":[[0,1,2,3,4,5]]}
```

## Docker images

The release provides three image archives:

| Image | Archive | Platform | Entry point | Contents |
| --- | --- | --- | --- | --- |
| `smsd:8.1.0-cpu-amd64` | `smsd-8.1.0-cpu-linux-amd64.tar.xz` | linux/amd64 | The command-line tool | Java 25 runtime, CPython 3.14.9 and the `smsd` package |
| `smsd:8.1.0-cpu-arm64` | `smsd-8.1.0-cpu-linux-arm64.tar.xz` | linux/arm64 | The command-line tool | As for amd64 |
| `smsd:8.1.0-cuda12.8` | `smsd-8.1.0-cuda-linux-amd64.tar.xz` | linux/amd64 | `python` | CPython 3.14.9, `smsd`, `smsd-cuda` and the CUDA 12 runtime; no Java |

Load a CPU image and pass the tool's options after the image name; use
`--entrypoint python` for Python:

```bash
docker load --input smsd-8.1.0-cpu-linux-amd64.tar.xz
docker run --rm smsd:8.1.0-cpu-amd64 --version
docker run --rm smsd:8.1.0-cpu-amd64 \
  --Q SMI --q "c1ccccc1" --T SMI --t "Oc1ccccc1" --json -
docker run --rm --entrypoint python smsd:8.1.0-cpu-amd64 \
  -c "import smsd; assert smsd.Substructure.find('CC').match('CCC').found"
```

Without options, the CPU image prints the tool's help. The container runs as
the user `smsd` in the working directory `/work`. If `docker load` stops with
`exec: "xz": executable file not found in $PATH`, the engine cannot
decompress the archive; decompress it on the host instead:

```bash
xz --decompress --stdout smsd-8.1.0-cpu-linux-amd64.tar.xz | docker load
```

The tool in the container reads only the files it is given. To search the
file of the [SDF batch search](#sdf-batch-search), mount the current directory
on `/work` and run as yourself, so the report belongs to you:

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" smsd:8.1.0-cpu-amd64 \
  --Q SMI --q "c1ccccc1" --T SDF --t targets.sdf --threads 2 --json benzene.json
```

The Python in each image holds only the SMSD packages, and in the GPU image the
CUDA runtime: there is no NumPy, pandas, RDKit or `mcp`. Packed fingerprint
indexes, `butina` and `maxmin` clustering, the data-frame and array adapters
and the molecular tools server therefore need an image built on top that adds
those packages.

The GPU image runs Python on an NVIDIA host with the NVIDIA Container Toolkit
and a compatible driver:

```bash
docker load --input smsd-8.1.0-cuda-linux-amd64.tar.xz
docker run --rm --gpus all smsd:8.1.0-cuda12.8 \
  -c "import smsd; print(smsd.__version__); print(smsd.gpu_device_info())"
```

As shipped, it serves prepared `TargetCorpus` screens, streamed screens and
`match_corpus(..., gpu_prefilter=True)`; a packed index on the GPU also needs
NumPy. That a backend is available does not prove that a call used the GPU;
the [CPU and GPU screening](CPU_GPU.md#containers) chapter explains how to
check. The [installation guide](../HOWTO-INSTALL.md#docker) describes building
the images from source.

## Molecular tools server

`smsd-mcp`, installed with `python -m pip install "smsd[mcp]==8.1.0"`, serves
eight molecular tools over standard input and output, with no network
listener. Without the extra it exits with status 1 after a Python traceback
that ends with
`ImportError: The MCP server requires the official SDK; install smsd[mcp]`.
The [molecular tools guide](../MCP_SERVER.md) lists the tools and their limits.

## See also

- [Java guide, command-line interface](../JAVA.md#command-line-interface): the
  reference for the tool and the Java library.
- [Installation guide](../HOWTO-INSTALL.md): source builds, the portable 8.0.3
  tool and the Docker images.
- [C++ and Java](OTHER_LANGUAGES.md): fingerprints and screening in the other
  languages.
- [CPU and GPU screening](CPU_GPU.md): where screening work runs.
