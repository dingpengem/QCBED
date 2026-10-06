# Scientific inputs

- `pattern.dat`: raw experimental FCC Al CBED pattern.
- `Al_FCC.cif`: FCC Al crystallographic structure input.

Start with `examples/fcc_al.yml` or a new GUI project. Data Prep regenerates
geometry, masks, variance and prepared descriptors from these inputs; no saved
preparation cache is required. Historical derived numerical references live in
`tests/fixtures/legacy_reference/`, solely for regression and developer comparisons.
