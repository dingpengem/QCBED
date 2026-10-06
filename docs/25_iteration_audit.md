# Audit of the apparent 25-iteration stop

The two completed runs in `tmp/al/` with IDs
`qcbed-run-20261006T014531571886Z-a81f1b3ba1ab4358b6f26c32b72a6c1f` and
`qcbed-run-20261006T014939320904Z-ddff24eacbee47febbae81fcd788299b`
have identical recorded optimizer results. This audit reads their existing
checkpoints and replays saved objective values only; it runs no new forward model.

**Stop reason: Evaluation budget reached. Not converged.**

| Item | Recorded / verified value |
|---|---|
| Workflow / stage | AMOEBA (`SIMP`), one stage, one workflow loop |
| Active targets | Depth, (020) real component, (200) real component; slice [0,1] |
| Dynamic vector dimension | 3 independent parameters |
| Completed optimizer iterations | 12 |
| Objective evaluations | 25; all accepted, none rejected |
| Configured iteration maximum | No separate iteration maximum in this historical AMOEBA run |
| Configured stage evaluation budget | `maxfsimp = 20`: legacy **soft search** budget |
| Initialization | N+1 = 4 simplex vertices |
| Search evaluations | 21 |
| Exact terminating branch | `amoeba._run`: `if numf > max_evaluations`, `stop_reason = "max_evaluations"` |
| Relative objective-spread tolerance | `ftol = 1e-6` |
| Final relative spread | 0.010221658140033863 |
| Objective threshold | `chimin = 1e-6`; final mean extreme objective = 0.912651412244863 |
| Convergence | False; neither convergence criterion satisfied |
| Global execution budgets | Not stored as readable values in these historical records; neither caused this normal stage completion |

The metric is `2*abs(y_max-y_min)/(abs(y_max)+abs(y_min))` over the final simplex.
The saved journal reproduces the original 12 iterations and all 25 callback values
exactly. Four initialization calls plus 21 search calls explain 25; **there is no
25-iteration limit** in this run. The legacy soft-budget branch checks `>` at cycle
boundaries and can overshoot. It is preserved rather than silently changing optimizer
behavior. A completed workflow does not necessarily mean convergence.

## Source of the effective budget

`modern_yaml.OPTIMIZER_DEFAULTS` defines `maxfsimp=20`; the new-project builder in
`gui/preparation_model.py` uses that profile. Compact YAML exposes this setting as
`optimizer.stage_budget` and restores 20 if it is omitted. The run checkpoint
confirms the effective value 20, although it cannot establish whether the user
explicitly entered 20 or accepted the default. The same small budget occurs in
instructional/preflight examples. The recommended per-target settings control
D1/D2/PLAMDA, not this global stage budget.

The packaged smoke uses an explicit manifest override (23 in the prior acceptance
manifest); that value is not installed as the new-project default. The relevant
production value **25** is checkpoint cadence (`checkpoint_every`), not a stop
condition. Some examples also use snapshot cadence 25. The AMOEBA adapter's
standalone fallback local budget is 40, but the grouped workflow passes the
configured 20 explicitly. Diagnostic ANNEAL profiles and other literal 25 values
in scattering tables or plotting code do not enter this run's stopping condition.
The current global defaults are 1,000 attempted evaluations per invocation and
3,600 seconds, but these are not substituted for the unavailable historical values.
No evidence supports a smoke-test budget of 25 leaking into the GUI. The existing
20 default is unchanged; its method-dependent meaning is now explicit in help.

## Actual recorded result

Latest, best and retained/final are identical for both runs:

| Physical target | Initial / IAM | Latest = best = retained | Hard bounds |
|---|---|---|---|
| Requested sample depth (Å) | 1000 | 1000.6111111111086 | [800,1200] |
| Re(E_200), slice [0,1] | 19.603768265997132 | 19.30969419192306 | [19,20] |
| Re(E_020), slice [0,1] | 19.603768265997132 | 19.660805303034163 | [19,20] |

Both elastic coefficients use **Real only**. Their imaginary components remain
fixed at −1.7602847776579752e−9. The initial atomic model has B = 0.85 Å², so its
IAM values differ from the tutorial's B = 0.86 Å² values. Depth is a specimen input,
not an IAM coefficient. The recorded gradient-signal objective is
**0.9079870068763699**; it must not be relabeled as a statistical chi-square without
applying the documented objective definition.

Optional Real-only equality would reduce depth + this pair from **3 to 2**
optimizer parameters. Complex equality gives depth + two shared components = **3**.
These optional constraints are unrelated to the termination reason.

## Reporting and explicit limits

Results and D now distinguish completed optimizer iterations, actual objective
evaluations, attempted proposals, and stop reason. Stage reports distinguish
convergence, soft evaluation budgets, exact iteration caps, cancellation and global
execution limits. New checkpoints record readable execution limits. Missing
historical counters are shown as unavailable instead of inferred from evaluations.

D → Advanced execution settings includes an optional **Maximum completed iterations
per stage** (`output.max_stage_iterations`). This is an additional exact cap; blank
preserves the existing engine behavior and legacy stage-budget semantics. A capped
stage retains its best recorded valid evaluation and can advance to the next stage.
The same-limit Resume replays the recorded prefix; changing the cap or equality
layout requires a new compatible run rather than silently extending a checkpoint.
Global evaluation/time interruptions remain resumable under their existing rules.
