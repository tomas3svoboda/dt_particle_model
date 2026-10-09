# Item 38: sparge steam in the desolventizer-toaster, what the records say

A technical note for the owner, 2026-09-27. Read-only on everything outside this
folder. No march was run, nothing was fitted and no physics claim is new. The only
computation is property arithmetic on the frozen authorities (IAPWS-95, Span-Wagner
n-hexane, the certified heteroazeotrope solver) in `steam_arithmetic.py`, output in
`steam_arithmetic.json` (source R23). No plant of this project is named or used; the
industrial numbers are from published sources only.

**Tags.** Each statement carries one tag and a source identifier from `SOURCES.md`:
**[MEASURED]** a record with a number from a run or lane; **[MODEL-DERIVED]** follows
from the formulation's closures or property authorities, with the pointer;
**[PUBLISHED]** a cited source; **[OPEN]** not established.

---

## Question 1. How sparge steam acts in the model

**Latent-heat supply.** Steam condenses on meal colder than its dew point and gives up
its latent heat there. Per kilogram, condensing steam at 100 °C releases 2.26 MJ, and
cooling the same steam by 27 K (from 127 °C) releases 0.055 MJ, so condensation
carries about 41 times more heat than the steam's sensible heat [MODEL-DERIVED R23].
The tower energy balance showed this on a lane of 2026-09-03: 6.46 MW of condensation
latent heat and 0.85 MW of steam sensible heat against 7.19 MW of hexane latent heat,
closing within 1.8 per cent, with the tray jackets at about 16 kW per tray, so more
than 98 per cent of the duty came from the sparge steam [MEASURED R8; that lane ran on
the operating basis later superseded (R6), but the tray-level balance stands]. The
published split is similar: direct steam gives about 75 per cent of the DT heat with
predesolventizing trays and about 90 per cent without [PUBLISHED R14 row B13]. One
kilogram of hexane evaporated at the co-boiling point needs 0.145 kg of condensing water
[MODEL-DERIVED R8, R23]; the lanes show 0.137 to 0.144 [MEASURED R9]. Condensation stops
when the meal reaches the boiling point of its sorbed water, 105.84 °C at 19.0 per cent
moisture and 0.993 bar, 6.4 K above pure water [MODEL-DERIVED R17, R21]; beyond that point
steam heats the meal by convection only. For convection alone the particle paper's
experiments are the reference case: in superheated hexane the surface stays at the
hexane boiling point, 68.71 °C, and Whitaker's correlation gives 0.82 to 0.98 of the
measured constant-rate duty [MODEL-DERIVED against PUBLISHED, R17, R22].

**Stripping agent.** Two effects are separate in the model. (i) *Co-boiling:* liquid
water next to liquid hexane lowers the boiling temperature from 68.71 °C to the
heteroazeotrope, 334.48 K (61.33 °C) on the certified solver [MODEL-DERIVED R2, R23];
published values are 60 to 62 °C [PUBLISHED karnofsky1985recovering via R14 row V5].
Sorbed water alone gives an intermediate plateau, 66.5 to 62.9 °C for 0.05 to 0.20 kg/kg
[MODEL-DERIVED R19]. Water then takes part of the heat, so at 136 °C the hexane loss
rate falls to 0.93 to 0.79 of the pure-hexane rate [MODEL-DERIVED R19]. (ii)
*Equilibrium floor:* the retained hexane at equilibrium depends on the hexane partial
pressure alone. At the last-tray mole fractions 0.003 to 0.03 (100 to 110 °C) the floor
is 4.5 to 66.0 ppm; at 0.1 it is 152 to 222 ppm [MODEL-DERIVED R17, R1]; at the
sparged-section gas Cardarelli gives (0.0064 at 105 °C) it is 11.6 ppm, 20.3 ppm with
the oil arm; at the 75 °C exit vapour (0.612) it is 4834 ppm, a gas the dry meal never
meets [MODEL-DERIVED R1]. Water vapour in the pore gas lowers the floor by the same
dilution: at 136 °C, sorbed water of 0.05 kg/kg takes it from 448.6 to 75.1 ppm
[MODEL-DERIVED R19]. The mass Biot number is 1.64e5, so the particle surface is at
equilibrium with the gas [PUBLISHED R3, Cardarelli Table 5.3]. The gas composition sets
the surface value, and diffusion through the dry shell sets the rate
[MODEL-DERIVED R17].

**Against recondensation of hexane.** A vapour saturated in water above the co-boiling
point is under-saturated in hexane: hexane activity 0.97 at 62 °C, 0.81 at 66 °C and
0.51 at 75 °C [MODEL-DERIVED R23]. Liquid hexane could only form on surfaces colder than
the hexane dew point of the gas: 59.5 °C for a 66 °C dome gas and 53.9 °C for the
published 75 °C exit vapour, whose feed enters at 55 °C [MODEL-DERIVED R23; PUBLISHED R3].
In the lanes, hexane never condenses: the gas-to-meal hexane flux is outward in all
7,866 layer-intervals, and only 4 of 650 film re-formations sit on the condensing side,
by at most 0.21 K [MEASURED R9].

**Steam against a hot inert gas at the same temperature.** What the model *can*
compare: an inert gas supplies only sensible heat. Steam's own sensible heat is
2.04 kJ/(kg K) [MODEL-DERIVED R23], so for the same heat at a 27 K approach a
non-condensing gas would need of the order of 40 times the mass flow of the condensing
steam [MODEL-DERIVED R23; air's heat capacity, about half of steam's, would double that;
the repository holds no air property authority]. At the same hexane partial pressure,
temperature and pressure the equilibrium floor is the same whatever the diluent, because
the frozen isotherm has no water argument [MODEL-DERIVED R19, R17]. What it *cannot*
compare: the pore gas is binary water-hexane, so a third gas is outside the formulation
[OPEN R17]. Water displacing sorbed hexane (80 per cent displaced at 15 per cent
moisture and 30 °C) is published but not modelled [PUBLISHED grant1983factors via R14;
OPEN]. So is condensed water sealing pores [PUBLISHED Mattea via R3; OPEN R10]. Meal
drying in an inert gas, and the role of moist heat in toasting, are outside the records
[OPEN].

## Question 2. How much condensed water becomes bound moisture, and does the sparge rate control it

**The routing in the model.** Condensate always goes first to a surface film. The film
soaks into the hexane-free outer shell up to the bound-water ceiling x_bf, instantly in
the nominal case and not at all in the lower limit. Water behind the front is inert
[MODEL-DERIVED R10, owner ruling D9-a]. A mixed layer reports the dry-matter-weighted
front of its packets [MODEL-DERIVED R11, D19]. The soaking rate has no source: it is
bracketed and swept, never fitted [OPEN R10].

**What the lanes measure.** Free water never exceeds 1.32 per cent of dry mass (mean
0.9 per cent). On SP1 layer 1, which takes about half the tower's condensation, the film
drains to zero within 13 plant seconds. The condensate, 107 kg per tonne fed, is 72 per
cent of the 150 kg per tonne of room in the shell, so in the model almost all condensed
water ends up as sorbed (bound) water within seconds [MEASURED R9]. At the loadings the
trays carry, bound water sits at activity 0.56 to 0.58, a factor of 1.75 below a free
film [MEASURED R12]. The tower condenses 0.46 of whatever sparge it is given: 0.462 at
70 kg/t and 0.466 at 487 kg/t. The shares are MN1 0.514, MN2 0.369 and SP1 0.117
[MEASURED R5, R4 item M-COND-1]. At the ruled reference point, 72.7 kg/t condenses from
152.9 kg/t, inside the published 52 to 76 [MEASURED R6]. These are first intervals of
profile starts, a startup transient.

**The particle paper.** Discharge moisture is an input, not a prediction [R17, subsection "The discharge
temperature and the boiling point of sorbed water"]. The paper reads the relation one way, and it holds the other
way too. In pure steam at 0.993 bar, meal at 105.80 °C holds the evidence cap,
0.2357 kg/kg (19.07 per cent wet basis), and meal at 107.92 °C holds 16.0 per cent
[MODEL-DERIVED R18, R20]. At equilibrium, the bound moisture at discharge is therefore
set by the meal's final temperature in steam, not by the steam flow. The steam rate acts
only through that temperature and through the time the meal gets to approach it
[MODEL-DERIVED R17]. Two limits apply. Against every published isotherm held, the water
branch sits at the retentive edge, which widens the point value to an interval of 100.5
to 107.6 °C [MODEL-DERIVED R21]. Above 19.07 per cent the model has no evidence
[R18].

**Does the sparge rate control it?** Physically, the heat demand sets how much steam
condenses. For 220 kg/t of hexane the model's arithmetic gives 31.8 kg/t of condensing
steam at no temperature rise and 75.7 to 82.1 kg/t with a 38 K rise, and any steam above
that should pass to the dome [MODEL-DERIVED R8]. In the model, the condensate grows in
proportion to the supply [MEASURED R5], which contradicts that picture. The cause is
measured: the gas leaves each layer at the layer-mean interface state [R5]. Published
whole-unit balances condense 0.58 of the direct steam [PUBLISHED Cardarelli 1998 via
R3] and 0.61 [PUBLISHED sipos1961dt via R14 row S8]. The "0.8 to 0.95" target of
M-COND-1 traces to the project's own 56 to 88 kg/t derivation, not to a source
[R15; R14 row B14 note; R7]. Whether 0.58 to 0.61 is the right target is
[OPEN]. Published operating trends: without predesolventizing trays the outlet moisture
is 20.6 per cent, with them 18.3, and it rises by 0.42 points per point of feed hexane
[PUBLISHED R14 row M4]. How discharge moisture depends on the sparge rate, and how the
water divides between bound and free at discharge, are [OPEN]. No tower lane reaches a
steady state: the production point marches to 3.2 plant s [MEASURED R16], against a
residence of 17 min [MEASURED R6].

## Question 3. A healthy water content for the dome vapour, and whether it can be near zero

**The model's statement.** The lanes' dome gas sits on the certified heteroazeotrope:
334.48 K, water mole fraction 0.205 at 101,325 Pa, 0.054 kg of water per kg of hexane
(the record rounds it to 0.055) [MEASURED R2; MODEL-DERIVED R23]. That is the coldest and
driest vapour that can leave over both free water and free hexane at one atmosphere.

**The published domes all ride the water dew line.** Published dome temperatures are 66
to 78 °C, typically 68 to 71 °C, with a theoretical floor of 62 °C [PUBLISHED R14 rows
V1, V2]. The published compositions by mass are 94/6 at 62 °C, 93/7 at 66, 92/8 at 68
to 71 and 91/9 at about 71 [PUBLISHED R14 row V3], plus 88.3/11.7 at 75 °C in Cardarelli
1998 [PUBLISHED R3]. Witte gives 0.05, 0.1, 0.2 and 0.46 kg water per kg hexane at 60,
70, 80 and 90 °C [PUBLISHED R3, Fig. 5.2]. Their water dew points come out at 63.8,
66.5, 68.9, 71.0 and 75.5 °C, within 1.8 K of the stated temperatures. The dew line
itself gives 0.093, 0.184 and 0.471 at 70, 80 and 90 °C, against Witte's 0.1, 0.2 and
0.46 [MODEL-DERIVED R23, ideal vapour]. The dome water load is therefore fixed by the
dome temperature: 0.054 at the lock, 0.073 at 66 °C, 0.099 at 71 °C, 0.129 at 75 °C and
0.184 at 80 °C [MODEL-DERIVED R23]. The dome temperature shows how much sparge steam
passed through uncondensed. The plant's dome is 13 K hotter than the model's because its
flash happens on the predesolventizing trays; the model's predesolventizing trays
evaporate nothing [MEASURED R2].

**The cost.** The latent heat carried by the dome water is 0.13 MJ per kg of hexane at
the lock, 0.18 at 66 °C, 0.23 at 71 °C and 0.31 at 75 °C. Hexane's own latent heat is
0.34 MJ/kg [MODEL-DERIVED R23]. One published plant study used 312.5 t/d of direct steam
at a 79 °C dome set-point and 271.1 t/d at 76 °C, 13 per cent less, with residual hexane
between 439 and 496 ppm at all four set-points [PUBLISHED avelar2018adequacao via R14
rows S10, M6]. Other plants lowered the dome from 72 to 75 °C to 68 to 70 °C [PUBLISHED
fu2024energy, zuo2024dtdc via R14 row V2]. Whether that heat is recovered downstream is
outside the records [OPEN].

**Can it be near zero?** Not while the gas leaves over liquid water. The minimum is the
lock value, 0.054 kg/kg, close to the published 62 °C floor at 6 per cent water (0.064 kg/kg, dew point 63.8 °C)
[MODEL-DERIVED R23; PUBLISHED R14 row V3]. Over hexane-wet meal that holds only sorbed
water, the co-boiling vapour carries 0.015 to 0.043 kg/kg for 0.05 to 0.20 kg/kg of
sorbed water (ideal vapour) [MODEL-DERIVED R19, converted in R23]. The feed brings 5 to
12 per cent water [PUBLISHED R14 row B1], so the load does not go to zero while the gas
last touches feed meal. Superheating the gas above its water dew point is cheap: 10 K
of hexane-vapour superheat costs 0.019 MJ/kg [MODEL-DERIVED R23]. It removes no water,
though; it only moves the gas off the dew line. Leaving with less water per kg of hexane
means condensing more of the steam in the bed and taking the dome towards the lock. That
does three things: more condensation heating of the feed, more moisture on the meal
(Question 2), and, below the gas's hexane dew point (59.5 °C for a 66 °C dome gas), a
risk of hexane recondensing [MODEL-DERIVED R23]. An optimum is [OPEN]. To quantify one,
the model would need a tower steady state, with the dome temperature, residual hexane
and outlet moisture as outputs against the sparge rate; the flash moved onto the
predesolventizing trays (the D12 re-issue); a gas that leaves through the colder top
meal (M-COND-1); and a cost function for steam, dryer duty and condenser load
[OPEN R2, R4, R16].

---

## Proposal: which statements go where

**The particle paper (`paper/01_particle_jfpe`).** Only statements that follow from the
particle formulation and its published comparisons:

1. *The floor against gas composition as the equilibrium side of stripping.* The
   statement is already in the floor subsection. It should first take item_34's
   correction: the sparged-section gas gives 11.6 ppm (20.3 ppm with the oil arm) and
   Cardarelli's Fig. 5.3 gas gives 6.8 to 47.1 ppm, replacing the mis-paired 0.612 at
   105 °C and its 1216.4 and 2085.5 ppm. Then add one sentence: the floor depends on the
   hexane partial pressure alone, so at equilibrium steam strips only by dilution, and
   water displacing sorbed hexane (Grant 1983) is the omitted co-sorption.
2. *The heteroazeotrope lock as pore-gas thermodynamics.* Free water beside free hexane
   co-boils at 334.48 K, 7.4 K below hexane's boiling point, with sorbed water giving
   the intermediate 62.9 to 66.5 °C plateau (item_14). This is Karnofsky's 60 to 62 °C
   with nothing fitted.
3. *The sorption-elevated boiling point read the other way.* In pure steam, the meal
   temperature fixes the equilibrium bound moisture: 19.07 per cent at 105.80 °C and
   16.0 per cent at 107.92 °C. This is stated as an equilibrium relation, with the
   moisture still an input and the kinetics not claimed.
4. Optionally, one sentence in the introduction: condensation carries about 41 times
   the steam's sensible heat per kilogram at a 27 K approach.

**The DTDC paper (`paper/03_dtdc_jfe`, discussion, missing-physics and measurement
sections).**

- The reading that published domes ride the water dew line (all within 1.8 K), the
  dome water load against temperature, and the latent heat it carries.
- The published steam trend against dome set-point.
- The model's dome on the lock, 13 K below the published plant, and why.
- The condensed fraction: 0.46 in the model regardless of supply, with its tray shares,
  against the published whole-unit 0.58 to 0.61. The M-COND-1 target should be
  corrected from 0.8 to 0.95.
- The D9-a routing and the measured free-water ceiling of 1.32 per cent.
- The hexane dew points of the gas as the recondensation margin.
- The heat comparison of steam against a hot inert gas.

**Measurement needed first.**

- The sparge-rate dependence of the condensed fraction and of discharge moisture: a
  tower steady state, or published sparge sweeps.
- The bound-to-free partition inside particles at discharge: transient moisture profiles
  in meal particles at DT temperature, already named in the particle paper.
- Water-hexane co-sorption at high water loading and DT temperature.
- A soybean-meal water adsorption isotherm between 50 and 100 °C.
- A ternary pore gas for any inert-gas comparison.
- Pore sealing by condensed water.
- Any optimum dome temperature.

`physically_qualifying: false. plant_predictive: false.`
