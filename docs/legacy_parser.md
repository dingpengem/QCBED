# FOR007 parser contract

## Scope and evidence

`parse_for007` implements the input prefix through `READ(7,452) EXPINT`.
It does not read experimental data from unit 12, or validate the later scan,
structure-factor and optimizer/refinement sections of FOR007. Refinement masks
and coefficients within the supported prefix are validated even when their
values are not exposed in `LegacyInput`. No physics is evaluated.

The separate `parse_for007_scan` API now reuses this prefix reader and then
reads initial scan basis/corners plus the NUMP/NUML experimental-file header.
It does not change `parse_for007` or parse intensities/refinement. See
[legacy scan](legacy_scan.md) for its explicit scope and singleton restriction.

The primary evidence is `private predecessor implementation`, PROGRAM MULTIS and
SFMAP. The related-version private predecessor implementation does not define parser string widths,
capacities or limits. Each row below corresponds to separate source READs;
parenthesized lists within a READ share one cursor.

| Source location | Reads and field interpretation |
| --- | --- |
| Lines 72–87 | Formatted title; real KV, integer KVOPT and three zone indices; integer NSLICE/MESHX/MESHY. |
| Lines 108–123 | Six real geometry fields per slice; integer NTYPE; FORMAT(A2) label and integer ZT per species. |
| Lines 147–166 | Real DSAP and DSMAX in separate reads; integer NREGIONS; interleaved real DEPTH/integer TMASK; integer RNSL list. |
| Lines 218–248 | Integer CONCFLAG; real vacancy coefficients alternating with integer masks; real DILATION and integer mask. |
| Lines 250–291 | Formatted ABSORB; conditional real absorption coefficients alternating with integer masks; per-region real DW list, integer mask list, then integer FSTYPE and three real coefficients. |
| Lines 1770–1774, 1844–1847 | Integer atom count and character RTYPE; S records read species, x, x mask, y, y mask, STARTOCC; G records omit masks, which source sets to zero. |
| Lines 339–364 | Formatted BMNAME; integer MBOUT and all integer H/K pairs in one READ; formatted EXPINT. |

## Supported syntax and record handling

- Whitespace and comma separators; blank records can be skipped while a
  list-directed READ still needs values. A trailing separator before continued
  values is accepted.
- Leading signs, integer/decimal real forms, E/e and D/d exponents, and an
  exponent without its letter (for example `1.25-2`). Leading signs and signs
  following E/e or D/d are accepted, for example `+0.5`, `1E+3`, and `2.0d-4`.
  A sign after an otherwise complete numeric literal is rejected: `1+`, `1-`,
  `0.5+`, and `2.0D3-` are invalid in both integer and real fields.
- Positive integer repetition counts with explicit values, such as `3*0`
  or `2*0.8D0`. Expansion is lazy and does not allocate the repeated list.
- Single- or double-quoted list-directed character values, with doubled
  delimiters recognized inside a string. The supported character field is
  RTYPE, whose decoded value must be exactly S or G, without normalization.
- Continuation is local to one READ, including MBOUT and its dependent beam
  list. Once that READ has enough values, remaining ordinary values on the
  last record are discarded, as in advancing Fortran input. They never leak
  into the next READ. The lexer checks the whole consumed record for
  unsupported syntax even if it occurs after the required values.

Formatted character records are not lexed. Full raw records are retained in
`title_raw`, `absorption_raw`, `output_filename_raw`, and
`experimental_intensity_filename_raw`. Existing title/filename convenience
fields remove trailing spaces only; leading spaces and quotes are preserved.
`Species.label_raw` retains its complete record, and `label` is the first two
characters padded with spaces, from the current source's FORMAT(A2).
Thus a label record `'AL'` yields the literal two-character label `'A`, not AL.
`AtomSite.record_type_raw` retains its original token, including quotes or
repetition; `record_type` stores the decoded S/G value.

Absorption is enabled only when `absorption_raw[:2] == "YE"`. Lowercase,
leading spaces and quoted YES do not activate it. No filename, title,
absorption or RTYPE width is imposed solely from the related-version include.

## Validation and intentional restrictions

All real fields in the supported prefix must be finite, including discarded
dilation, vacancy, absorption, Debye-Waller and surface coefficients. Beam
voltage, cell lengths, slice thickness and region depth must be positive.
Cell gamma must lie strictly between 0 and 180 degrees to describe a
nondegenerate direct cell. Fractional shifts and coordinates are not clipped
or wrapped; starting occupancy is checked for finiteness without inventing an
occupancy range.

Integers are parsed from integer syntax, never float-to-int conversion.
Counts and mesh dimensions must be positive, except atom and output-beam
counts may be zero. Species identifiers must be in 1..NTYPE. Atomic numbers
must be positive. RNSL entries must be positive and sum to NSLICE, consistent
with the driver's region-to-slice assignment. Reciprocal and zone indices may
be negative. No legacy fixed array capacities are copied.

Refinement masks and KVOPT accept 0/1. CONCFLAG and FSTYPE accept 0/1/2.
Absorption C's mask also accepts -1, explicitly documented and handled in the
source. These domain checks are intentionally stricter than a Fortran runtime
READ, which can accept integers that the scientific logic does not define.

Null values (leading/consecutive commas or repetitions without values), slash
termination, comment syntax, complex literals, and quoted numeric values are
explicit errors. Multirecord quoted strings and compiler-specific extensions
are unsupported. A leading comma on a continuation record is rejected; put
the separator after the preceding value instead. Decimal-comma numeric mode
is not supported: commas are separators.

Errors carry filename, section, record number, offending token/record and
expected type or constraint. Model validation is wrapped with the source
record of the implicated field. Input file read/decoding failures also receive
parser context. Unit tests exercise reference input and independent synthetic
cases; they do not establish executable Fortran numerical equivalence.
