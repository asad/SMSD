# Independent MCS oracles

`ilp_mcis_oracle.py` is a bounded mixed-integer model for checking small molecule
pairs independently of the search implementation. It uses optional NumPy and
SciPy dependencies. They are not required by SMSD itself.

The input is graph JSON with an `atoms` list containing `element` and `aromatic`,
and a `bonds` list of `[first, second, order, aromatic, ring]` entries. Atom indices
start at zero. `graph_from_smsd(graph)` copies this representation from a
`MolGraph` without running a search.

```sh
python benchmarks/ilp_mcis_oracle.py query.json target.json --seconds 15 --output result.json
```

The default is connected MCIS. `--disconnected` allows any number of fragments,
including single atoms. Both modes maximise atoms and then bonds. The chemistry
is the default molecule contract: elements match, bond order is strict and
aromaticity is flexible. Charges, isotopes, ring membership, stereo and tautomers
do not constrain the result. Other options are unsupported by this MCIS helper.

The model enforces injective atom correspondences and induced topology, including
the default aromatic bond and resonance rules. A single-root flow on selected
query edges enforces connectivity. The returned integer witness is checked
against the input graphs separately from the solver constraints.

Every call has explicit variable, row, nonzero, pair-check, atom and solver-node
limits. The time budget covers model construction and the solver cooperatively;
it is not a hard process or memory limit. `model_limit`, `limit_reached`,
`input_rejected`, `dependency_unavailable`, `invalid_witness` and `solver_error`
remain distinct from `optimal`. An existing output file is never replaced.

An optimal result means that HiGHS reported an optimum at zero requested MIP gap
and its integer witness and final bound passed the independent checks. It is not
a formally verified certificate. When no element-compatible atom pair exists,
the empty optimum follows directly without calling a solver. Validation on small pairs does not certify any
previously capped benchmark case.

The helper uses the documented SciPy [MILP interface](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.milp.html)
and [linear constraints](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.LinearConstraint.html).

`ilp_mces_oracle.py` is a separate ANY-bond MCES relaxation for fragmented searches
whose original minimum fragment size is at least two. Atom elements match; bond
order, aromaticity, charge, isotope, rings, stereo and tautomers do not constrain
the mapping. Binary selected-edge correspondences require both mapped endpoints,
injectivity and no isolated selected atoms. The objective is bonds first, then
atoms. With `W = min(query atoms, target atoms) + 1`, the scalar objective is
`W * selected bonds + selected atoms`.

```sh
python benchmarks/ilp_mces_oracle.py query.json target.json --min-fragment-atoms 3 --seconds 15 --output result.json
```

The model drops minimum fragment size and maximum fragment count. Its optimum or
solver dual bound is therefore an upper bound for the stated original contract.
The independently checked incumbent includes its selected-edge components and an
`original_contract_feasible` flag. An incumbent with a two-atom component does not
give an original-contract lower bound when the minimum is three; its score is
retained only as a relaxation incumbent. `original_contract_optimal` additionally
requires both solver optimality and original fragment feasibility.

Raw solver status, objective, numerical dual bound, MIP gap, node count and model
caps are retained even for stopped calls. A stopped numerical dual bound is
rounded upwards with the reported tolerance before deriving a bond upper bound;
it is not a formally verified certificate. The trivial graph-count upper bound
remains available if construction stops before the solver. A record labelled
`optimal` refers to the relaxation; inspect the original-contract fields before
using it to assess a previously capped search. `check_any_mces_witness` validates
an external retained atom/bond mapping under this exact chemistry without
calling a search or solver.
The MCES helper defaults to one HiGHS solver thread (`--solver-threads 1`). This
option is forwarded through SciPy; the exact solver options and any forwarding
warning are retained in each result. It does not cap threads used by unrelated
libraries or provide process-wide CPU isolation.
