# SMSD Pro user guide

This guide shows how to search, compare and screen molecules with SMSD Pro
8.1.0. Most chapters work in Python; the last two cover the command-line tool,
the Docker images and the C++ and Java libraries. Every Python example is a
complete program that runs on a CPU-only machine, and the output shown below
each one comes from running it on such a machine.

## Who the guide is for

The guide is written for chemists, data scientists and software developers who
need to find substructures and common substructures, compute fingerprints or
screen a compound library. It assumes a working knowledge of Python and of
SMILES, and no knowledge of how SMSD Pro searches. Each chapter starts from a
task and shows the calls that carry it out. The
[reference pages](#reference-pages) list every function and option.

## Chapters

| Chapter | What it covers |
| --- | --- |
| [Getting started](GETTING_STARTED.md) | Installation, a first substructure search and a first MCS search, and the conventions the other chapters rely on: reading a result, molecule input, atom numbering, default chemistry, time limits and threads |
| [Substructure search](SUBSTRUCTURE.md) | Prepared queries, one embedding or every embedding, many targets, prepared corpora, induced matching, chemistry options, hydrogen counts, limits, drawing a hit and R-group decomposition |
| [SMARTS queries](SMARTS.md) | Compiling a pattern, atom and bond primitives, recursive SMARTS, components, stereochemistry, the options that apply, SMARTS written from a molecule, limits and unsupported syntax |
| [Maximum common substructure](MCS.md) | MCIS and MCES, connected and fragment searches, proof and partial results, all optimal mappings, chemistry options, atom weights, reserved atoms, the common graph of a series, MCS similarity, export, drawing and tautomers |
| [Fingerprints and similarity](FINGERPRINTS.md) | Fingerprint families, count and configurable fingerprints, substructure keys and MACCS, many molecules at once, similarity coefficients, encodings and NumPy arrays |
| [Screening a library](SCREENING.md) | Choosing a route, preparing a library, similarity screens in memory and streamed from a file, an exact substructure screen, an MCS-similarity screen and saving the hits |
| [Packed fingerprint index](PACKED_INDEX.md) | Building and reopening an index file, threshold and nearest-neighbour searches, several queries in one pass, neighbours and pairs within the index, clustering, diverse selection and the file format |
| [CPU and GPU screening](CPU_GPU.md) | What runs where, checking the backend, installing the CUDA package, requesting or refusing the GPU, fallback and errors, and comparing the two devices |
| [Command line and containers](COMMAND_LINE.md) | The `smsd` command-line tool, SDF batch search, mapping export, exit status, the JSON report, the Docker images and the local molecular tools server |
| [C++ and Java](OTHER_LANGUAGES.md) | The same searches in the C++ and Java libraries, and the differences between the languages |

## How the chapters fit together

Read [Getting started](GETTING_STARTED.md) first: it installs the package and
sets out the conventions that every later chapter assumes. The next three
chapters cover the exact graph searches.
[Substructure search](SUBSTRUCTURE.md) asks whether a query occurs whole inside
a target, [SMARTS queries](SMARTS.md) describes the pattern language for the
same searches, and [Maximum common substructure](MCS.md) finds the largest
graph that two molecules share. [Fingerprints and similarity](FINGERPRINTS.md)
turns molecules into fingerprints and compares them.

[Screening a library](SCREENING.md) brings these together: fingerprint screens
narrow a library, and the exact searches decide which molecules match.
[Packed fingerprint index](PACKED_INDEX.md) keeps a library's fingerprints in a
file that is built once and searched in later sessions, and
[CPU and GPU screening](CPU_GPU.md) explains which screening steps can run on
an NVIDIA GPU. The last two chapters leave Python:
[Command line and containers](COMMAND_LINE.md) runs single searches and SDF
batch searches from a shell or a Docker image, and
[C++ and Java](OTHER_LANGUAGES.md) shows the same searches in the other
libraries. After Getting started, each chapter can be read on its own, and the
chapters link to each other where a topic continues.

## How the reference chapters are laid out

The chapters from Substructure search to Command line and containers follow a
common pattern, so a detail is found in the same place in each:

- what the feature does and when to use it;
- how to use it, in the order the calls are usually made;
- each kind of search or operation, with a runnable example and its output;
- an options section that lists the methods, arguments or settings with their
  defaults;
- notes on behaviour that is easy to miss;
- see also, with links to related chapters and reference pages.

## Reading a result

Every search returns a `SearchResult`, and four of its fields say what the
result proves:

- `status`: why the search ended, one of `COMPLETED`, `TIMED_OUT`,
  `CANCELLED` or `LIMIT_REACHED`;
- `optimal`: whether the answer is proved; for a substructure search, a match
  exists or a completed search proved absence, and for an MCS, the maximum is
  proved;
- `complete`: whether the search, or an enumeration of every match, ran to its
  end;
- `found`: whether any match was returned.

An empty result proves absence only when `status` is `COMPLETED` and `optimal`
is `True`. A search stopped by its time limit, node limit or cancellation proves
nothing about what it did not find.
[Getting started](GETTING_STARTED.md#reading-a-result) describes the other
fields and shows a stopped search.

## Running the examples

Install SMSD Pro 8.1.0 as shown in
[Getting started](GETTING_STARTED.md#install). Each Python block is a complete
program: copy it into a file and run it with `python`. The values each program
prints are also checked with `assert` statements. Programs that write files
keep them in a temporary directory, and none downloads anything.

A few programs need more than the core package:

- programs that use a packed fingerprint index, clustering or NumPy arrays
  need NumPy, from the `numpy` extra;
- two programs take RDKit molecules, from the `rdkit` extra;
- the GPU branches in [CPU and GPU screening](CPU_GPU.md) run only where the
  `smsd-cuda` package and a usable NVIDIA GPU are present, and the same
  programs take their CPU path everywhere else.

The C++ and Java programs build as described in
[C++ and Java](OTHER_LANGUAGES.md). In a source checkout,
`scripts/check_examples.py` runs every Python, C++ and Java block of this guide
together with the other guides and the example notebooks;
[Check the documented examples](../HOWTO-INSTALL.md#check-the-documented-examples)
gives the command.

## Notebooks

| Notebook | Topic |
| --- | --- |
| [Substructure and SMARTS search](../../python/examples/substructure_search.ipynb) | Molecule queries and SMARTS patterns on a small library of drugs and benzodiazepines: hydrogen counts, functional groups, enumeration, chemistry options, limits, corpus decisions and drawing a hit |
| [Maximum common substructure](../../python/examples/mcs_search.ipynb) | The two objectives, stopped searches, all mappings, chemistry rules, MCS similarity, SMARTS export and the common graph of a short series |
| [Fingerprints and similarity](../../python/examples/fingerprints_similarity.ipynb) | Every fingerprint family, the binary and count coefficients, configurable fingerprints and a similarity ranking of benzodiazepines |
| [Screening a library end to end](../../python/examples/library_screening.ipynb) | Five screens of one library, clustering of the packed index and saving the hits to an SD file |
| [CPU and GPU screening](../../python/examples/cpu_gpu_screening.ipynb) | The screening backend, automatic and requested GPU use, packed-index backends, fallback and a CPU and GPU comparison, with every cell able to run without a GPU |
| [MCS to a small library search](../../python/examples/feature_search.ipynb) | A common graph turned into a SMARTS query, a library search, SVG highlighting, NumPy features, packed-index neighbours, Butina clustering and MaxMin selection |

Open a notebook in Jupyter with a kernel that has SMSD 8.1.0 or later
installed. The notebooks are saved without outputs, perform no downloads,
remove any files they write and report when an optional package is missing.
To run one in Colab, copy the setup cell of
[feature_search.ipynb](../../python/examples/feature_search.ipynb) to the top
of the notebook and run it first; it installs an uploaded SMSD wheel.

## Reference pages

- [Python guide](../PYTHON.md): every function, class and option of the `smsd`
  package.
- [API reference](../API_REFERENCE.md): interfaces, options, proof and stopping
  status, input and output, fingerprints, adapters and GPU scope.
- [Search algorithms](../ALGORITHMS.md): search guarantees, result semantics
  and the differences between the languages.
- [Packed fingerprint corpora](../PACKED_CORPORA.md): the full packed-index API,
  its CPU and CUDA backends, ordering, work and memory.
- [Python examples](../EXAMPLES.md): fingerprints, depiction, layout, stereo,
  files and R-groups, beyond searching.
- [C++ guide](../CPP.md): the header-only C++ library.
- [Java guide](../JAVA.md): the Java library and its command-line tool.
- [Installation](../HOWTO-INSTALL.md): platform requirements, downloads and
  source builds.
