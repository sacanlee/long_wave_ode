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

On 2026-09-02 the repository was extended with (i) a **technical report on the
book's Section 8.3.3 errata** (`docs/Technical_Report_Ch8_Model_Errata.docx`),
(ii) a **cycles figure of the nonlinear ODE solution at the book parameters**
(`figures/ode_solution_book_params.png`, produced by
`scripts/long_wave_model.py --cycles-fig`), and (iii) a full **wholesale-price
evidence package** (gold-denominated and local-currency WPI for five countries,
1800–2026: scripts, regenerated JSON series and figures). See
["What was added in this update"](#what-was-added-in-this-update).

On 2026-09-06 the repository was extended with a **class-struggle extension of
the model** (`scripts/class_struggle_ode.py` + `docs/class_struggle_report.md`):
it checks the critique that the book's closure g_e ≅ −sV (Eq. 8.15) eliminates
the *autonomous* (above-reserve-army) wage push from the five-equation model,
introduces the class-struggle level CS via g_e = −sV − CS, and solves the ODE
for a grid of CS values, measuring the effect on the **cycle length** and on
the **upswing/downswing lengths** (see
["Class struggle in the long-wave model (CS extension)"](#class-struggle-in-the-long-wave-model-cs-extension)).

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
5. **Nonlinear solution at the book parameters** (initial values = closed-form
   equilibrium (8.25) × 1.05): within a 300-year window sV, sC, δ and τ show
   **median peak-to-peak periods of ~51–54 y over their first 4–5 cycles**,
   while the rate of profit r has no complete cycle (it decays 0.18 → 0.015,
   the tendency of the profit rate to fall). The oscillation never settles
   into a stationary cycle: the closed-form equilibrium is an unstable focus
   (Re = +0.0212), so amplitudes grow ~×3 per cycle and the spacing lengthens
   (sC: 43 → 82 y over 300 y; up to ~140 y by 600 y, after which sV turns
   negative and the model leaves its economic domain) — see
   `figures/ode_solution_book_params.png` (per-panel annotations plus a
   bottom strip from a 600-year integration) and
   `docs/Technical_Report_Ch8_Model_Errata.docx` (Finding 6).

## Repository layout

```
long-wave-model/
├── README.md
├── requirements.txt
├── LICENSE
├── scripts/                       # all analysis scripts (English)
│   ├── long_wave_model.py         # model implementation, book self-checks, n scenarios,
│   │                              #   cycles_by_n figure + optional --cycles-fig ODE-solution figure
│   ├── sv_sign_flip.py            # sign-flip robustness check of Eq. (8.15) (4 variants × 6 findings)
│   ├── class_struggle_ode.py      # class-struggle CS extension of (8.15): cycle length /
│   │                              #   upswing-downswing lengths vs CS + historical-episode pulse
│   ├── history_compare.py         # model vs the 5 historical long waves
│   ├── forecast_2050_2080.py      # 2050-2080 UN WPP 2024 three-scenario forecast
│   ├── ssa_fertility_scenario.py  # SSA fertility-crash scenario, waves 6/7 re-estimation
│   ├── wpi_build_series.py        # WPI: assemble raw downloads -> annual/quarterly
│   │                              #   gold-denominated series + cycles (JSON)
│   ├── wpi_plot_charts.py         # WPI: gold-denominated country + combined figures
│   ├── wpi_plot_local_currency.py # WPI: local-currency figure, pre-1940 peaks/troughs
│   ├── wpi_dl_destatis.py         # WPI: download helper for the slow Destatis server
│   ├── wpi_em_build.py            # EM (10 developing countries): price/PPP series ->
│   │                              #   PPP-corrected gold-denominated index (JSON)
│   └── wpi_em_sync.py             # EM vs US cycle-sync detection + combined JPG figure
├── docs/                          # English analysis reports
│   ├── model_validation_n_scenarios.md
│   ├── history_comparison.md
│   ├── forecast_2050_2080.md
│   ├── ssa_fertility_scenario.md
│   ├── sv_sign_flip_report.md
│   ├── class_struggle_report.md  # is the class struggle really absent from the model?
│   │                             #   + CS extension, cycle/upswing/downswing results,
│   │                             #   historical reading of every CS level
│   └── Technical_Report_Ch8_Model_Errata.docx   # 6 findings on the book's Section 8.3.3
├── data/                          # regenerated result JSONs (English)
│   ├── results.json
│   ├── history_results.json
│   ├── forecast_results.json
│   ├── ssa_scenario_results.json
│   ├── sv_flip_results.json
│   ├── class_struggle_results.json    # CS scenarios: cycle/rise/fall stats + pulse episode
│   ├── wpi_annual_1800_2026.json      # annual WPI / FX / gold-denominated WPI
│   ├── wpi_quarterly_1982_2026.json   # quarterly gold-denominated WPI
│   ├── wpi_cycles_gold.json           # gold WPI peak/trough cycles (9-yr MA)
│   ├── wpi_pre1940_cycles_local.json  # local-currency WPI pre-1940 cycles
│   ├── wpi_em_annual.json             # EM price/PPP/PPP-corrected gold index (annual)
│   └── wpi_em_sync.json               # EM vs US sync-start years + correlations
└── figures/
    ├── cycles_by_n.png            # linearised Δr(t) per n scenario
    ├── ode_solution_book_params.png   # nonlinear ODE solution cycles at book parameters
    ├── class_struggle_trajectories.png  # CS runs: r/sV/sC paths up to the domain exit
    ├── class_struggle_cycle_metrics.png # cycle length / upswing-downswing / rise share /
    │                               #   profitability as functions of CS
    ├── class_struggle_pulse_episode.png # stylised wage-explosion -> offensive episode vs CS=0
    ├── class_struggle_ratchet_recoveries.png # pointless recoveries at fixed CS=+0.02:
    │                               #   the wave keeps cycling while the profit rate ratchets down
    ├── class_struggle_oscillation_death.png # high-CS regime: r -> ~1e-4-1e-2, cycles die out,
    │                               #   no sV<0 purge, premise sC+sV<=1 violated (extrapolated)
    ├── wpi_chart_US.png ... wpi_chart_FR.png   # gold-denominated WPI per country
    ├── wpi_chart_combined_9yMA_log.png         # five countries, 9-yr MAs
    ├── wpi_chart_local_currency.png            # local-currency WPI, 5 panels
    └── wpi_em_gold_denominated_sync.jpg        # 10 EM + US, sync-start years marked
```

## Usage

```bash
pip install -r requirements.txt

python scripts/long_wave_model.py                        # self-checks + n scenarios
python scripts/long_wave_model.py --ns 0.015 0 -0.008 -0.015
python scripts/long_wave_model.py --cycles-fig figures/ode_solution_book_params.png
python scripts/sv_sign_flip.py                           # sign-flip check of (8.15)
python scripts/class_struggle_ode.py                     # CS scenarios + pulse episode
python scripts/class_struggle_ode.py --cs -0.05 0 0.05   # custom CS grid
python scripts/history_compare.py                        # historical comparison
python scripts/forecast_2050_2080.py                     # 2050-2080 forecast
python scripts/ssa_fertility_scenario.py --data-dir <dir>  # SSA scenario (needs OWID CSVs, see below)
python scripts/wpi_build_series.py --raw-dir <dir>       # WPI series (needs downloads, see below)
python scripts/wpi_plot_charts.py                        # gold-denominated WPI figures
python scripts/wpi_plot_local_currency.py                # local-currency WPI figure
python scripts/wpi_em_build.py --gold-csv wpi_raw_data/gold_monthly.csv   # EM series (needs network)
python scripts/wpi_em_sync.py                            # EM vs US sync + JPG figure
```

Results are written to `data/` (JSON) and `figures/` (PNG/JPG).

### Input data

The first five scripts run standalone (all model inputs are hard-coded book
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

The WPI scripts need the externally downloaded source files described in
["WPI input data"](#wpi-input-data); they are **not** stored in the repository
(raw data ~120 MB). Put them in one directory (default `<repo>/wpi_raw_data/`,
git-ignored) and pass it via `--raw-dir`.

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

---

# What was added in this update

Merged on 2026-09-02 from the working project *price_index_in_gold* (task
folder) plus the redo material of the Section 8.3.3 technical report:

1. `docs/Technical_Report_Ch8_Model_Errata.docx` — the technical report
   (6 findings) on the internal inconsistencies of the book's Section 8.3.3
   (equilibrium / closed form / linearisation / "about 50 years" claim).
   Finding 6 confirms that the nonlinear ODE solution nevertheless shows
   ~50-year cycles with a long-run fall of the profit rate.
2. `figures/ode_solution_book_params.png` — cycles of the **nonlinear solution
   of the ODE** at the book parameters
   (a1=0.04, a2=0.01, b0=0.004, b1=0.005, b2=0.05, gw=0.03, n=0.015),
   generated by the `--cycles-fig` option of `scripts/long_wave_model.py`
   (5 variable panels over 300 y with per-panel annotations of the measured
   peak-to-peak spacings — initial values = closed-form equilibrium (8.25) ×
   1.05 — plus a bottom strip that plots every completed-cycle spacing from a
   600-year integration and shows that the period keeps lengthening, i.e. no
   stable cycle exists; footnote explains the cause, the unstable focus).
3. A wholesale-price evidence package (5 countries, 1800–2026):
   - **scripts** `wpi_build_series.py`, `wpi_plot_charts.py`,
     `wpi_plot_local_currency.py`, `wpi_dl_destatis.py` — see
     ["WPI scripts"](#wpi-scripts);
   - **data** `wpi_annual_1800_2026.json`, `wpi_quarterly_1982_2026.json`,
     `wpi_cycles_gold.json`, `wpi_pre1940_cycles_local.json` — regenerated
     result series (JSON, following the repo's data convention);
   - **figures** `wpi_chart_US/UK/DE/JP/FR.png`,
     `wpi_chart_combined_9yMA_log.png`, `wpi_chart_local_currency.png` — see
     ["WPI figures"](#wpi-figures).

## WPI figures

All WPI figures are English-language PNGs. They can be re-generated with the
plot scripts (they read the JSON series in `data/`; no downloads needed for
plotting).

| Figure | Content | Data used (sources) |
|---|---|---|
| `figures/wpi_chart_US.png` … `wpi_chart_FR.png` | Gold-denominated WPI per country, log scale, 1913=100: annual index (grey), centred 3/6/9-yr moving averages (blue/orange/red), peak/trough years labelled on the 9-yr MA | The five WPI series in local currency spliced over 1800–2026 (see source table below) divided by the USD exchange rate and the USD gold price; rebased to 1913=100 |
| `figures/wpi_chart_combined_9yMA_log.png` | Five countries overlaid (9-yr MAs, log), 39 peak/trough markers with collision-free year labels | Same as above |
| `figures/wpi_chart_local_currency.png` | Local-currency WPI, 3×2 subplot grid: every country rebased to its own 1913=100 (log). Pre-1940 span shaded; peaks/troughs of the 9-yr MA labelled there only. Germany panel caps the y-axis at 1e5 because the 1923 hyperinflation annual average (1.26e11) is off scale (annotated with an arrow + note) | The published local-currency WPI series (sources below); no FX/gold needed |
| `figures/ode_solution_book_params.png` | Nonlinear ODE solution (r, sV, sC, δ, τ) at the book parameters over 300 y: sV/sC/δ/τ median spacings ~51–54 y over the first 4–5 cycles, lengthening to ~80–140 y (bottom strip from a 600-y integration) — no stable cycle; r has none (decays to ~0.015) | Model equations (8.15)–(8.19) of the book; parameters as printed (no empirical data) |

**Cycle tables printed on the charts** (peaks `P` / troughs `T`, 9-yr MA,
min. 3-yr spacing, swing ≥ 12%, alternating):

- Gold-denominated WPI: US P 1814/1869/1921/1967/2001, T 1846/1895/1938/1984/2016;
  UK P 1812/1921/1967/2001, T 1895/1936/1983; DE P 1857/1876/1918/1968/2001,
  T 1866/1895/1953/1984/2014; JP P 1922/1969/1998, T 1896/1938/1983;
  FR P 1920/1967/1998, T 1938/1984/2014 (values in `data/wpi_cycles_gold.json`).
- Local-currency WPI before 1940: US P 1814/1868/1921, T 1846/1895/1935;
  UK P 1810/1870/1920, T 1847/1895/1934; DE P 1871/1919, T 1895/1935;
  JP P 1922, T 1932; FR P 1927, T 1933
  (`data/wpi_pre1940_cycles_local.json`; Germany 1919 peak is off scale on the
  chart; Germany 1945–47 has no data; a first trough within 5 years of the
  series start is discarded as an edge artefact).

### WPI data sources (used by the figures above)

The index input series are a mix of FRED/NBER, national statistical offices and
historical compilations, spliced at overlaps (details of every splicing anchor
are in the script comments of `scripts/wpi_build_series.py`):

**Wholesale / producer prices (local currency):**

| Country | Historical segment | Source |
|---|---|---|
| US | 1800–1849 | Warren–Pearson wholesale price index (HSUS E23–42) |
| US | 1850–1894 / 1890–1914 | FRED NBER M0448AUSM323NNBR / M0448BUSM336NNBR |
| US | 1913–2026 | FRED PPIACO — BLS all-commodity PPI (1982=100) |
| UK | 1800–2016 | Bank of England "Three centuries" millennium dataset, table A47 (wages and prices; 2015=100) |
| UK | 1948–2026 (monthly splice) | ONS GB7S output PPI (2015=100) + OECD MEI (GBRPROINDMISMEI) via FRED |
| DE | 1851–1914 | Hamburg wholesale price index — FRED NBER A04054DE00HAMA314NNBR (annual) / M04054DE00HAMM314NNBR (monthly) |
| DE | 1914–1944 | Statistisches Reichsamt (1913=100; via FRB bulletins and the literature), incl. 1922–23 hyperinflation (1923 annual average 1.26e11) |
| DE | 1948–1957 | literature approximation of the Bundesbank long series (±10%, noted in the charts) |
| DE | 1958–2024 / 2010–2026 | OECD MEI DEUPROINDMISMEI via FRED / Destatis statistical report table 61241-b01 (2021=100) — slow server, fetch with `wpi_dl_destatis.py` |
| JP | 1887–1946 | BOJ IMES historical Tokyo wholesale prices (monthly CSVs, Shift-JIS) |
| JP | 1945–1959 | BOJ/MOF national WPI (1934–36=100; level-spliced to the IMES segment at 1945–46) |
| JP | 1955–2024 / 2024–2026 | OECD MEI JPNPROINDMISMEI via FRED / BOJ CGPI monthly-change tables in the monthly-report PDFs |
| FR | 1913–1939 | FRED NBER M04057FRM360NNBR |
| FR | 1940–1945 / 1946–1952 | FRB bulletin (1938=100) / INSEE (1938=100) |
| FR | 1956–2024 / 2005–2025-10 | OECD MEI FRAPROINDMISMEI via FRED / INSEE IPPI monthly (2015=100; release paused after 2025-10) |

**Gold price (USD/oz, used by the gold-denominated figures only):**
1792–1833 official $19.39; 1834–1861 $20.67; 1862–1878 $20.67 × greenback
premium (Mitchell 1908 / HSUS); 1879–1933 $20.67; 1934–1967 $35; 1950–2026
monthly market price (London fixing / COMEX averages cached from the commodity
price project).

**Exchange rates (local currency per USD, gold figures only):** UK — BOE A33
(USD/GBP 1791–2016, inverted) then FRED EXUSUK; Germany — gold parity 4.198
to 1912, MeasuringWorth historical series 1913–1970 (official hyperinflation
values for 1923–24 and 1941–45), FRED EXGEUS (DM) then DEXUSEU (EUR) 1971+;
Japan — Hitotsubashi LTES table 26 (USD/100 yen, 1874–1941) with official
1942–48 rates, MeasuringWorth ~360 (1949–1970), FRED EXJPUS 1971+; France —
gold parity 5.1826 to 1912, MeasuringWorth 1913–1970, FRED EXFRUS (new francs,
×100 for the 1960 redenomination) then DEXUSEU × 655.957 old francs per euro.

Known caveats (also in the figure notes and JSON metadata): series definitions
change over time (wholesale → producer prices); level breaks of a few percent
possible at the splicing points; Germany 1945–47 has no data (occupation);
Japan 1945–59 prices were controlled and the official FX understates the black
market; France 1953–55 is linearly interpolated.

### WPI scripts

| Script | Function |
|---|---|
| `wpi_build_series.py` | Reads the downloaded source files (`--raw-dir`), assembles the annual and quarterly series (FX and gold denominate the WPI), detects the 9-yr-MA peak/trough cycles, and writes `data/wpi_annual_1800_2026.json`, `data/wpi_quarterly_1982_2026.json`, `data/wpi_cycles_gold.json` |
| `wpi_plot_charts.py` | Reads the JSON series and draws the per-country gold-denominated charts plus the combined 9-yr-MA chart (log axes, labelled peaks/troughs, collision-avoided year labels) into `figures/` |
| `wpi_plot_local_currency.py` | Draws the 5-panel local-currency figure with the shaded pre-1940 window, labelled pre-1940 peaks/troughs and the German capped axis; also writes `data/wpi_pre1940_cycles_local.json` |
| `wpi_dl_destatis.py` | Downloads the Destatis producer-price xlsx with retries + integrity check (that server is too slow for plain curl) |

All WPI figures were regenerated from the committed JSON data and verified to be
pixel-identical to the working-project originals.

# Emerging-economy gold-denominated WPI and sync with the United States

Added 2026-09-02. Ten developing economies (China, India, Brazil, Mexico,
Kenya, Nigeria, Egypt, Indonesia, Pakistan, Morocco) are compared with the
United States on the question: *from which year did their long-wave (price)
cycle start to move in sync with the US?*

## What the scripts do

| Script | Function |
|---|---|
| `wpi_em_build.py` | Downloads price indices and PPP factors, builds the annual local-currency price index, backcasts PPP before 1990, computes the PPP-corrected gold-denominated index and writes `data/wpi_em_annual.json` |
| `wpi_em_sync.py` | Correlates each country with the US gold-denominated WPI (from `data/wpi_annual_1800_2026.json`, rebased to 2010=100), detects the sync-start year and draws `figures/wpi_em_gold_denominated_sync.jpg`; writes `data/wpi_em_sync.json` |

Method and formula (as requested):

```
goldden_c(t) = price_index_c(t) / [PPP_c(t) × gold_usd(t)],   rebased 2010 = 100
```

- The PPP factor plays the role of the "fair" exchange rate, so the series
  survives the episodes of excessive currency depreciation of developing
  economies (raw USD-FX-denominated indices collapse there).
- Sync-start criterion: Pearson correlation over a rolling 20-year window
  between the centred 9-yr MA of log(goldden) and the same transform of the
  US series; the country "syncs from year T" if the correlation is ≥ 0.5 for
  15 consecutive years starting at T. A growth-rate variant (1-yr log
  differences) is stored as a robustness check.

## Data sources (used by the figure)

| Series | Source | Coverage |
|---|---|---|
| Brazil price index | FGV IGP-DI, chained from monthly % changes via the Banco Central do Brasil SGS API (series 190, Feb 1944+); general price index with ~60% wholesale weights — the WPI proxy | 1950–2025 |
| Other 9 countries: price index | consumer price index, World Bank `FP.CPI.TOTL` (2010=100) | 1960–2025 (China 1986–2025) |
| PPP conversion factor (LCU per int. $) | World Bank `PA.NUS.PPP` | 1990–2025, backcast before 1990 |
| US CPI (for the backcast) | FRED `CPIAUCSL` annual averages | 1950+ |
| Gold (USD/oz, annual) | monthly market prices from `wpi_raw_data/gold_monthly.csv` (git-ignored; a static copy of the derived annual gold series is embedded in `data/wpi_em_annual.json`) | 1950–2025 |
| US reference curve | gold-denominated WPI of the US from `data/wpi_annual_1800_2026.json` | rebased 2010=100 |

PPP backcast rule (no official PPP exists before 1990):

```
PPP(t) = PPP(1990) × [pi(t)/pi(1990)] / [CPI_US(t)/CPI_US(1990)]   for t < 1990
```

i.e. a constant 1990 real exchange rate — a documented approximation.

## Results (sync-start year vs the US)

| Country | Sync from | Country | Sync from |
|---|---|---|---|
| Brazil | 1969 (earliest measurable, data from 1950) | Nigeria | 1979 (earliest measurable, data from 1960) |
| India | 1979 (earliest measurable) | Egypt | 1979 (earliest measurable) |
| Mexico | 1979 (earliest measurable) | Indonesia | 1979 (earliest measurable) |
| Kenya | 1979 (earliest measurable) | Pakistan | 1979 (earliest measurable) |
| Morocco | 1979 (earliest measurable) | **China** | **2005** (the only genuine late sync) |

Because the rolling 20-y window and the 15-y persistence rule cannot be
evaluated before the data begin, "earliest measurable" means the country
already satisfied the criterion from the first possible evaluation year —
i.e. it was in sync with the US for the whole observable period. China is the
exception: its correlation with the US only crossed and stayed above 0.5 from
2005 (windows covering ~1990–2019+).

## Caveats

1. The official WPI/PPI portals of most of these countries (India OEA/RBI,
   Mexico INEGI/Banxico, China NBS, Egypt CAPMAS, ...) could not be fetched
   headlessly from this environment (SSL/HTTP blocks). Only Brazil therefore
   uses a wholesale-based index (FGV IGP-DI); the other nine use the CPI.
   Swapping in an official WPI/PPI series later is a drop-in replacement of
   the `pi` block of `data/wpi_em_annual.json` (rerun `wpi_em_sync.py`).
2. For CPI-based countries the ratio CPI/PPP tracks the US price level
   closely (PPP absorbs most of the domestic inflation), so their
   gold-denominated series are dominated by the common "US CPI in gold"
   factor. Their sync results mostly reflect the world real gold cycle, and
   the pre-1990 PPP backcast assumes a constant real exchange rate. The
   Brazil (wholesale-based) series and the post-1990 official-PPP decades are
   the least affected parts of the dataset.
3. Series end in 2025 (full calendar years; latest World Bank vintage).

# Class struggle in the long-wave model (CS extension)

Added 2026-09-06. Question (see `docs/class_struggle_report.md` for the full
analysis): the book's 3rd case (g_e < 0, "profit-squeeze", p. 212) treats the
growth rate of the rate of surplus value as the carrier of the class struggle,
but the closure of the complete model sets g_e ≅ −sV (p. 214), a pure function
of the investment in variable capital — so an *autonomous* wage push (a
struggle aiming higher than the reserve army dictates) has no place in
(8.15)–(8.19). **Verdict of the report: the critique is essentially correct**
(with the reserve-army wage push surviving inside g_e = −sV, and with the
observation that the level of distribution is not a state variable at all).

The extension puts the struggle back as **g_e = −sV − CS**, where CS is the
class-struggle level (CS = 0 reproduces the book; CS > 0 = above-normal wage
push — the "historical and moral element" of the value of labour-power;
CS < 0 = below-normal push — Marx's counteracting cause, the depression of
wages below value). Only Eq. (8.15) changes:

```
(8.15-CS)  r' = −a1·(sC − sV)·r² + a2·(δ + τ − sV − CS)·r
```

Main quantitative results at the book parameters (11 CS values from −0.10 to
+0.10; cycles measured within the economically meaningful window):

| Effect of raising CS (more militant) | Magnitude |
|---|---|
| Linearised (eigenvalue) period | essentially unchanged: 41.7 → 41.8 y |
| Nonlinear cycle length | lengthens ≈0.05 y per +0.01 CS on the first and ≈0.16–0.18 y per +0.01 CS on the second completed cycle (second cycle: 45.0 y at CS=−0.10 → 46.6 y at CS=0 → 48.4 y at CS=+0.10; third ~56 y at CS=+0.10); shortening under CS < 0 |
| Upswing (trough→peak) of sC | lengthens ~0.2–0.3 y per +0.01 CS at moderate/high CS (23.7 → 28.3 y across the grid) |
| Downswing (peak→trough) of sC | flat to slightly shorter (19.9 → 19.4 y) |
| Rise share of the wave | 54% (CS = −0.10) → 57% (CS = 0) → 59% (CS = +0.10) |
| Profit-rate level | the dominant effect: the CS component of d ln r/dt is −a₂·CS (−0.01·CS per year), amplified by state feedback; mean r [0,120 y]: 0.150 → 0.141 → 0.133 across the grid |
| Extreme CS | CS ≥ +0.2: cycles stretch to 50–74 y, r collapses toward 1–4% (permanent squeeze / stagnation); CS ≤ −0.2: cycles ~41–43 y, profits kept high (Roaring-Twenties-like repression) |

Historical reading of the levels (directions, not calibrated magnitudes; units:
CS = 0 means the wage share does not keep pace with productivity under
reserve-army pressure, cf. book p. 216):
CS > 0 ↔ wages at or above the productivity norm — the golden age on the
fourth-wave upswing (1940–1966; moderate positive CS, high profit rate),
intensifying into the 1966–1975 wage explosion and the fourth-wave profit
squeeze (1966–1982); CS ≈ 0 ↔ reserve-army-dominated 19th-century waves;
CS < 0 ↔ the post-1980 offensive (union-density collapse, wage share decline,
profit recovery, the 1982–2007 upswing); CS ≪ 0 ↔ the US 1920s;
regime-scale CS > 0 is historically unsustainable and appears only as bounded
episodes (France 1936–38, UK 1974–79) that crises terminate. A stylised
"wage-explosion → offensive" pulse run (`figures/class_struggle_pulse_episode.png`)
shows that, because the book weights the distribution channel with a₂ = 0.01
against a₁ = 0.04 for the OCC channel, even a strong historical-style episode
moves the profit rate by only ~±1%: the model is OCC-determined by
construction — a quantitative confirmation of the critique. At a *fixed*
CS = +0.02 the model produces "pointless recoveries"
(`figures/class_struggle_ratchet_recoveries.png`): the accumulation wave keeps
cycling (boom peaks at t ≈ 50/97/150/213 y) while every boom peaks at a lower
profit rate (r: 0.146 → 0.115 → 0.083 → 0.051) — short, profitless recoveries
into a permanently lower profitability, ended only by the regime change
(CS < 0) or by the model leaving its domain (sC > 1 at ≈86 y, sV < 0 at
≈133 y). See the report for caveats (no distribution state variable, no true
fixed point, ~2-cycle meaningful window) and the suggested sixth-state
extension.

## License

MIT — see `LICENSE`. Note: the model equations belong to the cited book
(Chatzarakis, Tsaliki & Tsoulfidis 2022, Routledge); this repository contains an
independent numerical reproduction and analysis.
