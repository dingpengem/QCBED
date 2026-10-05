# Immutable variance-closure geometry fixture

These are archived geometry and baseline hashes, not newly fitted geometry.
`beam_catalogue.csv` is copied byte-for-byte from the validated integer-center A7
catalogue. Image crop paths are historical display metadata, not runtime inputs.
`geometry.json` extracts A8 origin/physical radius, A10 conversion, and pre-variance
signal/weight manifests from the previous Patch A closure. Original paths and
SHA256 hashes are recorded; large A8 search traces and synthetic objective reports
are deliberately omitted.

The original pattern/sigma stay in read-only `reference_data`; their hashes are
verified. No physical hkl indexing is asserted. Tests and the diagnostic example
use these fixtures so they do not depend on local temporary run directories.
Do not regenerate them to hide a geometry, beam-order or signal regression.
