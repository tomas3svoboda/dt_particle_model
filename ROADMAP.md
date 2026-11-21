# Closing the two PD-reference qualification gaps

Current reference: complete conservation-audited 64-cell L0 history. Finer
birth/time refinement and PD film applicability remain open. This roadmap
does not require a coupled washer/tower or a population balance.

## 1. Establish the PD transfer boundary

Audit density, heat capacity, viscosity, conductivity, **bulk gas** binary
diffusivity, velocity, Re/Pr/Sc and the definitions of Nu/Sh at the declared
source state. Keep pore diffusivity separate. Count the gas exposure area
reduction once and the jacket duty once.

The current Sc interval 0.6--1.8 is an engineering envelope, not proof that
the physical vapor at Sc about 0.35 is impossible. Also examine Reynolds
range and geometry: the existing bed coefficient comes from thin-bed
through-flow evidence, whereas PD has a stirred wet bed on solid disks
exposed to a freeboard. Passing a numerical guard is not physical validation.

Use primary-source evidence to select a PD coefficient prescription in the
actual regime, retaining native finite-Stefan/component/energy accounting.
Compare the existing source-based bed law with an appropriate freeboard or
agitated-bed treatment. A single-sphere correlation is only a scale comparison
unless its bed/exposure mapping is justified. Record length/velocity, property
evaluation state, mixing and geometry assumptions explicitly. Do not tune
coefficients or exposure to recover a preferred outlet.

**Exit:** one source/property/domain table and one declared PD boundary with
a small physically supported sensitivity bracket. If no suitable correlation
exists, identify an independently measured effective coefficient as missing
evidence; do not manufacture a validity interval or silently widen the guard.

## 2. Diagnose and fix birth on saved states

Reuse the saved wet-feed/birth IVPs. First hold exactly the same physical IVP
for the accepted 0.01-s and refused 0.005-s births. Test 0.0025 s only after
the mechanism is isolated. Record component/energy jump residuals, film
energy and temperature, phase bounds, beta, conditioning and work.
Do not repeat the previous solver/seed sweep without new evidence.

Distinguish failed nonlinear search from an incompatible state/storage
handoff. A carries external water at body temperature; B stores it at the
solved surface temperature. A finite relaxation over shrinking birth time is
a specific target, not a proven sole cause. Previous conservative
surface-film controls also refused, so adopting that candidate is not a fix.

If an admissible root exists, correct residual evaluation, scaling or event
continuation using the existing equations and work budget. Preserve rejected
path rollback. If the handoff itself is incompatible, derive a conservative,
thermodynamically compatible A/B film state. A separate film node can coexist
with a uniform particle body; any such model change must be explicit and
reviewed. Never cool the checkpoint, delete film energy, widen beta or call
the coarser step convergence.

**Exit:** admissible, conserved, smoothly refining local births and accepted
B advances, with negative-path/rollback checks intact.

Steps 1 and 2 can be investigated independently in a serial queue. Freeze
their final boundary and any justified handoff correction together before
expensive whole-history qualification.

## 3. One full history, then matched refinement

First run one fresh cold-feed-to-outlet history under the final boundary.
Check admissible inventories, event order, external-hexane depletion before
core recession, source/film/energy ledgers, vapor export and solver work.
Track uptake, boiling, PD removal, core extinction, film depletion and H/W/T.

Then use genuinely nested radial meshes (proposed 32/64/128) with matched
time controls, and time halving at fixed mesh. Begin with base/half-step
64-cell histories and extend only to separate spatial, temporal and event
errors. Hold uniform-body A fixed and state that radial refinement does not
test its approximation; time refinement still covers A and the handoff.
Freeze observables, event alignment and required reporting accuracy before
running. Conservation alone is not a discretization error estimate.

No tolerance, phase band, seed family or work budget is changed to get a pass.
Any correction needs measured justification. The current 132.8-ppm outlet
is not a calibration target; the result is allowed to change.

## 4. Update the paper and deposit together

Report final boundary provenance/applicability and matched numerical errors.
Recalculate sensitivity claims retained for the final bath, or remove them.
Regenerate all article numbers, plots, audits and manifests. Verify a clean
public reproduction package and publish a new numerical release; retain the
current reference as history.

**Resources:** one numerical worker and one heavy command at a time, half
CPU affinity, low priority and single-threaded numerical libraries. On the
owner's Windows workstation preserve the existing 35% and 25-GiB RAM reserves
and cooperative stops. Measure local costs before estimating campaign time.
This roadmap launches no full numerical campaign.

Source leads: Faner (2008), printed pp. 106/109 (existing source intake);
[Coletto et al. (2022)](https://doi.org/10.1016/j.jfoodeng.2021.110870);
PD topology in [US5992050A](https://patents.google.com/patent/US5992050A/en).
[Wakao and Funazkri (1978)](https://www.sciencedirect.com/science/article/pii/0009250978851203)
is a packed-bed comparison lead, not a selected replacement or proof of PD
applicability. Check its domain and geometry before considering its use.
