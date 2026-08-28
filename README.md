# Classical Long-Wave Model: Reproduction, Validation & Demographic Scenarios

**Author**: sacanlee

Reproduction, verification, and scenario analysis of the five-variable long-wave
(Kondratiev) model of

> Chatzarakis, N., Tsaliki, P., & Tsoulfidis, L. (2022).
> *Economic Growth and Long Cycles: A Classical Political Economy Approach*.
> Routledge. — Section **8.3.3 "Third Stage: The Complete Model"** (Eqs. 8.15–8.19).

The project (1) verifies the script against the numbers printed in the book —
uncovering several **internal print inconsistencies** — (2) evaluates the model
against the five historical Kondratiev waves (book Table 5.1), and (3) projects
the long-wave period for 2050–2080 under UN WPP 2024 working-age population
scenarios, including a Sub-Saharan-Africa fertility-crash variant and a
sign-flip robustness check of Eq. (8.15).

## The model

Five key variables: the rate of profit `r`, the shares of surplus value invested
in variable (`sV`) and constant (`sC`) capital, the devaluation rate `δ` and the
rate of technical change `τ`:

```
(8.15)  r'   = −a1·(sC − sV)·r² + a2·(δ + τ − sV)·r     profit-rate dynamics (with counteracting forces)
(8.16)  sV'  = sC·r − δ − τ − n                         employment / industrial reserve army (Goodwin-type)
(8.17)  sC'  = sC·(τ − sV·r)                            constant-capital investment (mechanisation adoption)
(8.18)  δ'   = b0 + δ·(gw − sC·r)                       devaluation
(8.19)  τ'   = b1·Erf(δ − b2)                           pace of technical change
```

Book parameters: a1=0.04, a2=0.01, b0=0.004, b1=0.005, b2=0.05, n=0.015, gw=0.03.

## Key findings

1. **The book's printed equilibrium {0.0947, 0.1188, 0.8708, 0.05, 0.045} does
   not satisfy its own equations** (max residual 0.0294, coming from (8.16)/(8.17)),
   and its printed closed form (8.25) gives {0.1756, 0.1650, 0.6263, 0.05, 0.045}.
   The printed linearisation matrix is internally self-consistent but implies a
   period of **66.3 years**, not the "about 50 years" of the text. This project
   uses the book's own method (closed-form equilibrium (8.25) + analytical
   Jacobian), which yields **41.7 years** at n=0.015 — the correct order of
   magnitude.
2. **Flipping the sign of `sV` in (8.15) does not fix the inconsistencies**
   (see `docs/sv_sign_flip_report.md`): the F1–F3 problems are unaffected by the
   sign; all sign variants remain unstable at their equilibria; and economically
   the minus sign encodes the Marxian reserve-army channel (Marx, *Capital* I,
   ch. 25 → Glombowski 1983, g_e = ε₀ − ε₁sV), so a flip has no theoretical basis.
3. **Historical fit**: with working-age (15–64) population growth n, the model
   matches the modern waves closely (W4 1940–1982: +3%, W5 1982–2022: +7%),
   underestimates the 19th-century waves (−14% to −35%, parameters are
   modern-calibrated), and — nonlinearly — predicts the historical
   "recovery > recession" asymmetry correctly (W2/W3 exact at 55%).
4. **2050–2080 projection**: under UN WPP 2024, global 15–64 peaks in **2070**
   and turns negative only mildly (2050–80 n = +0.048%/yr). The model gives a
   sixth-wave period of **31–39 years** (most likely ~33–35 y), down from ~40 y
   in Wave 5 — faster, slightly downswing-heavy. Deep shortening (25–28 y) is a
   post-2085 phenomenon. The SSA fertility-crash scenario moves the working-age
   peak 23 years earlier (2047) but barely changes wave lengths (−1 y), because
   wave length is insensitive to second-order changes in n.

## Repository layout

```
long-wave-model/
├── README.md
├── requirements.txt
├── LICENSE
├── scripts/                       # all analysis scripts (English)
│   ├── long_wave_model.py         # model implementation, book self-checks, n scenarios, figure
│   ├── sv_sign_flip.py            # sign-flip robustness check of Eq. (8.15) (4 variants × 6 findings)
│   ├── history_compare.py         # model vs the 5 historical long waves
│   ├── forecast_2050_2080.py      # 2050-2080 UN WPP 2024 three-scenario forecast
│   └── ssa_fertility_scenario.py  # SSA fertility-crash scenario, waves 6/7 re-estimation
├── docs/                          # English analysis reports
│   ├── model_validation_n_scenarios.md
│   ├── history_comparison.md
│   ├── forecast_2050_2080.md
│   ├── ssa_fertility_scenario.md
│   └── sv_sign_flip_report.md
├── data/                          # regenerated result JSONs (English)
│   ├── results.json
│   ├── history_results.json
│   ├── forecast_results.json
│   ├── ssa_scenario_results.json
│   └── sv_flip_results.json
└── figures/
    └── cycles_by_n.png            # profit-rate deviation r(t) per n scenario
```

## Usage

```bash
pip install -r requirements.txt

python scripts/long_wave_model.py                        # self-checks + n scenarios
python scripts/long_wave_model.py --ns 0.015 0 -0.008 -0.015
python scripts/sv_sign_flip.py                           # sign-flip check of (8.15)
python scripts/history_compare.py                        # historical comparison
python scripts/forecast_2050_2080.py                     # 2050-2080 forecast
python scripts/ssa_fertility_scenario.py --data-dir <dir>  # SSA scenario (needs OWID CSVs, see below)
```

Results are written to `data/` (JSON) and `figures/` (PNG).

### Input data

The first four scripts run standalone (all model inputs are hard-coded book
parameters and the wave periodisation from the book's Table 5.1 / App. 8.A).
`ssa_fertility_scenario.py` additionally needs two OWID CSVs (UN WPP 2024),
downloadable from Our World in Data:

- `population-by-age-group-with-projections.csv`
- `fertility-rate-with-projections.csv`

(`curl -L https://ourworldindata.org/grapher/<chart>.csv` works; pass the
directory via `--data-dir`.) The `n` inputs for `history_compare.py` /
`forecast_2050_2080.py` were computed from these OWID/World Bank series on
2026-08-23 (working-age 15–64 CAGR; pre-1950 shares imputed, see
`docs/history_comparison.md` §2.1).

## Main results (quick reference)

| Scenario | n (15–64, annual) | Period (eigenvalue) | Nonlinear sC period |
|---|---|---|---|
| Book baseline n = +0.015 | +1.50% | 41.7 y | 49.4 y |
| n = 0 | 0 | 32.6 y | 36.2 y |
| n = −0.008 (mild decline) | −0.80% | 28.4 y | 30.3 y |
| n = −0.015 (deep decline) | −1.50% | 24.9 y | 25.8 y |
| Optimistic 2050–80 (UN High) | +0.464% | 35.3 y | 38.9 y |
| Medium 2050–80 (UN WPP 2024) | +0.048% | 32.9 y | 36.6 y |
| Pessimistic 2050–80 (UN Low) | −0.302% | 31.0 y | 33.8 y |

## License

MIT — see `LICENSE`. Note: the model equations belong to the cited book
(Chatzarakis, Tsaliki & Tsoulfidis 2022, Routledge); this repository contains an
independent numerical reproduction and analysis.
