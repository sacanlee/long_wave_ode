# Sub-Saharan Africa Fertility-Crash Scenario: Global Working-Age Population and Re-estimation of Waves 6/7

**Author**: sacanlee | **Date**: 2026-08-23. Script: `scripts/ssa_fertility_scenario.py` (all computations);
data: `data/ssa_scenario_results.json`.

## 1. Scenario definition

**Assumption**: Sub-Saharan Africa's (SSA, UN M49 definition, 48 economies) total
fertility rate (TFR) falls **linearly from its current level to below 2.1 over
2025–2045 (20 years; target 2.0)**, then returns to the UN medium path.

Data: UN WPP 2024 (medium) — OWID `fertility-rate-with-projections` (country TFRs)
+ `population-by-age-group-with-projections` (country age structures); SSA
aggregated by country weighting.

**Baseline facts (SSA aggregate)**: 2025 population **1.271 bn, weighted TFR 4.24**;
under the UN medium scenario 2050 is 2.090 bn with TFR 2.87 — **the UN itself
projects a rapid decline**; this scenario merely moves "reaching 2.1" from the
UN's ~2060s to 2045.

## 2. Projection method (transparent simplifications)

1. SSA annual births b(c) ≈ 0–4 population (year c+2)/0.95 (child-survival
   approximation, constant first order);
2. Scenario TFR path: 4.24 in 2025, linear to 2.0 in 2045, then = UN path;
3. Birth-reduction ratio q(c) = 1 − TFR_scenario(c)/TFR_UN(c), applied to each
   2025–2045 birth cohort; the reduction passes at the same ratio into the
   15–64 surviving population (mortality unchanged at the UN path);
4. World 15–64(Y) = UN value − Σ{2025≤c≤2045, c+15≤Y≤c+64} q(c)·b(c)·0.93.

## 3. Results: world working-age population

| Indicator | UN medium | SSA scenario | Diff |
|---|---|---|---|
| Cumulative 15–64 birth gap (2025–2045 cohorts) | — | ≈ **0.70 bn** | — |
| 15–64 population in 2070 | 6.286 bn | **5.639 bn** | **−0.65 bn (−10.3%)** |
| 15–64 population in 2100 | 6.064 bn | 5.543 bn | −0.52 bn |
| **15–64 peak year** | **2070** | **≈2047** | **23 years earlier** |

| Window n (15–64) | UN medium | SSA scenario | Change |
|---|---|---|---|
| 2030–2040 | +0.543% | +0.543% | 0 (cohorts not yet entering) |
| 2030–2045 | +0.496% | +0.463% | −0.033pp |
| **2045–2070** | **+0.162%** | **−0.252%** | **−0.414pp (turns negative!)** |
| 2070–2085 | −0.137% | −0.153% | −0.016pp |
| 2085–2100 | −0.103% | **+0.039%** | +0.142pp (gap cohorts start to exit; echo effect) |
| 2040–2070 | +0.202% | −0.160% | −0.362pp |

**Interpretation**: from the 2050s, the fertility-gap cohorts (born 2025–45)
enter working age in batches, switching the global 15–64 from "slow growth" to
"net contraction" (2045–70: +0.16% → −0.25%); **the working-age peak moves from
2070 to 2047**; after 2085 the reduced cohorts begin exiting at age 64, producing
an "echo effect" (growth recovers as the gap flows out with cohort ageing).

## 4. Re-estimation of Waves 6/7 (vs the UN medium baseline)

Method: the rise/fall segments of each wave use the **corrected n** above; the
long-wave model (book Eqs. 8.15–8.19, book parameters): main measure =
eigenvalue period × nonlinear (sC) rise/fall shares.

| Wave | Phase n (scenario) | Rise | **Peak** | Fall | **Total period** | Trough |
|---|---|---|---|---|---|---|
| **Wave 6 (from 2030)** | rise +0.463% / fall **−0.252%** | 17.7 y | **≈2048** | 15.8 y | **33.5 y** | ≈2063 |
| UN baseline | +0.495% / +0.163% | 17.8 y | ≈2048 | 16.7 y | 34.5 y | ≈2065 |
| **Wave 7 (from 2063)** | rise −0.026% / fall −0.031% | 16.2 y | **≈2079** | 16.3 y | **32.5 y** | ≈2095 |
| UN baseline (from 2070) | −0.137% / −0.103% | 15.8 y | ≈2086–2087 | 16.1 y | 32.0 y | ≈2102 |

*(Note: within the SSA scenario, Wave 7's start is chained internally to the
Wave 6 trough ≈2063 — unlike the UN baseline, which anchors 2070 (the 15–64
peak); the periods themselves are nearly identical (32.5 vs 32.0 y); only the
peak/trough years shift with the start.)*

## 5. Conclusions

1. **Wave 6 barely moves (33.5 vs 34.5 y; same peak ~2048)**: fertility-decline
   cohorts take ≥15 years to enter working age; their impact hits the labour
   force after 2045, not the 2030–45 upswing — instead it knocks the Wave 6
   **downswing** from "still positive (+0.16%)" to "net contraction (−0.25%)",
   slightly compressing the fall (16.7→15.8 y); the total period shortens by
   about 1 year.
2. **The real change: the working-age peak comes 23 years earlier
   (2070→2047)**, and **the global 15–64 population shrinks by 0.65 bn (−10%)
   during 2047–2070** — the Wave 6 peak (2048) lands exactly on the
   working-age peak (2047): the demographic inflection point coincides with the
   wave peak, and the downslope half is driven by "net labour-force
   contraction".
3. **Wave 7 (from 2063) is almost unchanged (32.5 y)**: chained to the Wave 6
   trough of 2063, its rise/fall n are both ≈ −0.03% (echo effect: reduced
   cohorts start exiting at 64, offsetting the new-entry gap) — the scenario's
   net effect decays with cohort flow; **peak ≈2079, trough ≈2095, period within
   <2% of the UN baseline (32.0 y).** The two waves chain seamlessly:
   Wave 6 2030–2063 (33.5 y) → Wave 7 2063–2095 (32.5 y).
4. **Quantitative implications for the 2050 narrative**: SSA fertility reaching
   target 20 years earlier ⇒ the global working-age "peak–turn-negative"
   inflection comes about a quarter-century earlier, but **the long-wave rhythm
   shortens by only ~1 year (Wave 6) / unchanged (Wave 7)**. Long-wave length is
   insensitive to *second-order changes* in the labour-force growth rate; it
   follows mainly the **level** of n and gw — consistent with the book model's
   historical assessment that wave length is dominated by other parameters.
5. **Limitations**: births backcast from 0–4 population (±5% noise), fixed
   survival coefficient 0.93, migration and participation not endogenous, TFR
   returning to the UN path after 2045; the model's printed inconsistencies
   remain (see the validation document). Directions are robust (earlier peak,
   2045–70 turning negative, Wave 6 slightly shorter); precise years are
   order-of-magnitude judgments.

### 5.1 Final chained timetable (SSA scenario, trough-to-trough)

| Long wave | Trough | Peak | Trough | Total period |
|---|---|---|---|---|
| Wave 6 | 2030 | **≈2048** | ≈2063 | **33.5 y** |
| Wave 7 | 2063 | **≈2079** | ≈2095 | **32.5 y** |

*(UN baseline for comparison: Wave 6 2030→2048→2065 (34.5 y); Wave 7
2070→2086→2102 (32.0 y). Net effect of the SSA crash scenario: Wave 6 −1 y,
Wave 7 +0.5 y — nearly no change at the period level; the change is concentrated
in "the working-age peak coming 23 years earlier (2070→2047)" and "Wave 6's
downswing turning negative (−0.252%)".)*

## 6. Reproduction

```bash
python scripts/ssa_fertility_scenario.py --data-dir <dir>  # all tables above + data/ssa_scenario_results.json
# Data: OWID fertility-rate-with-projections.csv / population-by-age-group-with-projections.csv (UN WPP 2024)
# Parameters adjustable at the top of the script (target TFR=2.0, window 2025-2045, survival=0.93)
```
