# SMSD (CDK‑compat) — Pure‑Java CLI & API

Pure‑Java substructure & MCS search using CDK (no external toolkits).  Jar‑first usage is prioritised; launchers are generated at build.

## Build
```bash
mvn -U clean package
```

## Run (fat JAR)
```bash
java -jar target/smsd-cdk-3.0.0-jar-with-dependencies.jar sub   --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json - --json-pretty
```

## Or use the generated launcher
```bash
target/bin/smsd sub --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json - --json-pretty
```

### Example output (no icons; JSON first policy)
```json
{
  "query": "CCN",
  "target": "CCCNC",
  "mode": "substructure",
  "mcs_type": "MCCS",
  "induced": false,
  "mappings": [
    {
      "index": 0,
      "query_atoms": [0, 1, 2],
      "target_atoms": [1, 2, 3],
      "target_bonds": [1, 2],
      "pairs": [[0,1],[1,2],[2,3]],
      "subgraph_smiles": "[NH]C[CH2]"
    }
  ]
}
```

## Repository layout
- `src/main/java/...` — core code (CDK‑only)
- `src/scripts/` — launchers copied to `target/bin/` at build
- `src/test/java/...` — JUnit tests
- `legacy/` — preserved prior code (read‑only)
- `WHITEPAPER.md` — algorithms & design
- Community docs: `LICENSE`, `WHITEPAPER.md`

## Notes
- Output style is deliberately plain
- To apply SPDX/Apache headers to all Java files. 
