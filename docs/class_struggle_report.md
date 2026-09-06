# Class Struggle and the Long-Wave Model: A CS-Extension of Eqs. (8.15)–(8.19)

**Author**: sacanlee | **Date**: 2026-09-06 | **Script**: `scripts/class_struggle_ode.py`
**Results**: `data/class_struggle_results.json` | **Figures**:
`figures/class_struggle_trajectories.png`, `figures/class_struggle_cycle_metrics.png`,
`figures/class_struggle_pulse_episode.png`

**Reference**: Chatzarakis, N., Tsaliki, P. & Tsoulfidis, L. (2022), *Economic Growth and
Long Cycles: A Classical Political Economy Approach*, Routledge, Ch. 8 (book page numbers
below refer to this edition; in this repository the model equations (8.15)–(8.19) and the
book's Section 8.3.3 have already been reproduced and audited — see
`docs/model_validation_n_scenarios.md`, `docs/sv_sign_flip_report.md` and
`docs/Technical_Report_Ch8_Model_Errata.docx`, findings F1–F6).

**Question posed** (job_0906.txt): the author treats a negative growth rate of the rate of
surplus value, g_e < 0, as the "profit-squeeze" case, i.e. as an active push for higher wages,
but then closes the model with g_e ≅ −sV, a purely mechanical function of the investment in
variable capital sV. Does this eliminate the *conscious* class struggle (a wage push aiming
higher than what the investment in variable capital dictates) from the five-equation model of
Section 8.3? And what happens if a class-struggle parameter CS is put back into the equation —
in particular to the **length of the cycle** and to the **lengths of its upswing and downswing**?

---

## Part I — Is the critique right?

### I.1 What the book actually says, page by page

**1. The "profit-squeeze" case (g_e < 0) is a case with *autonomous* distributional content.**
The 3rd case of the OCC–profit-rate analysis appears at the end of Section 8.2.1.2 (book
p. 212), immediately before Section 8.2.1.3 (book pp. 213–215):

> "3rd Case: g_e < 0, that is, the surplus value is reduced relative to variable capital
> (the case of profit-squeeze)." (p. 212)

The book is explicit that the *sign* of g_e is where class relations enter. Already on p. 207:

> "The bifurcation of equilibrium point B is, therefore, defined by the sign of g_e, which is
> a variable related not only to productive forces of the system (i.e., technological change,
> capital accumulation, etc.) but also to relations of production (i.e., income distribution,
> class struggle, exploitation of labor, etc.)."

and on p. 212, discussing what the g_e < 0 case means for the law of the tendential fall in
the rate of profit (FROP):

> "The significance of point B is only theoretical and may be interpreted in favor of the
> 'profit squeeze' argument, according to which **the intensity of the class struggle leading
> to wage squeeze of profits is the reason behind the Law of the FROP**. Under these
> circumstances, the Law of FROP cannot be viewed as a general economic law independent of
> people's 'will' and capable of realistically describing the dynamics of a capitalist economy."

So at this stage of the book g_e is effectively an **exogenous shift parameter** whose sign and
size are set, among other things, by the intensity of class struggle.

**2. Section 8.2.1.3 then endogenises g_e mechanically.** Book p. 214:

> "According to Marx (*Capital* I, ch. 25), the growth rate of the rate of surplus value, g_e,
> is inversely proportional to the change in the labor force, or, what is the same, to the
> investment in variable capital, sV. Glombowski (1983), for example, hypothesized that
> g_e = ε₀ − ε₁sV, where ε₀ and ε₁ are small positive constants. In our case, we specify a
> somewhat simpler relation according to which g_e ≅ −sV and any remaining differences are
> captured by two small positive parameters a₁ and a₂..."

The profit-rate equation (8.14) contains g_e additively, g_r = ṙ/r = −(sC − sV)r + δ + τ + g_e,
so the substitution g_e ≅ −sV produces the printed Eq. (8.15):

**ṙ = −a₁·(sC − sV)·r² + a₂·(δ + τ − sV)·r** (8.15)

with a₁ = 0.04, a₂ = 0.01. Equation (8.15) is one of the five equations (8.15)–(8.19) of the
complete model of Section 8.3 (book pp. 221–235). g_e itself appears in **no other equation**
of the system: (8.16)–(8.19) contain only r, sV, sC, δ, τ and the constants n, g_w, b₀, b₁, b₂.

**3. The reserve-army channel that *remains*.** Strictly speaking the model does not abolish
wage pressure altogether. Because sV is high when capital accumulation into employment exceeds
the natural growth of the labour force (Eq. 8.16: ṡV = sC·r − δ − τ − n), g_e ≅ −sV means the
rate of surplus value grows slowly (or falls) exactly in the tight-labour phases of the boom —
the classical reserve-army / Goodwin mechanism, which the book itself describes on p. 216:
"The pressure from an increasing IRAL keeps the wage share constant or even decreasing,
increases the discipline of the labor force, and leads to a higher level of labor intensity."
This is a **wage push, but one fully dictated by the state of the labour market**: workers get
higher wages only when, and exactly as much as, accumulation of variable capital dictates.

**4. The critique is therefore essentially right, with three qualifications:**

- **(i) The case taxonomy (g_e < 0 = profit squeeze) is in Section 8.2.1.2, not 8.2.1.3** —
  the user's section reference is off by one subsection, but the two are contiguous and the
  substantive point is unaffected (the profit-squeeze case immediately precedes the g_e ≅ −sV
  closure).
- **(ii) The tension is real and visible inside the book.** The book repeatedly assigns class
  struggle a decisive role *in words* — the sign of g_e on p. 207, the FROP-via-profit-squeeze
  on p. 212, g_e as "the key explanatory variable of the upswings and downswings of the rate
  of profit" (pp. 237–238), and the rate of surplus value as "a key distributive variable
  reflecting the level of class struggle in the sphere of circulation" (p. 239) — but in the
  operative closure of Section 8.3 g_e is replaced by −sV. A working class that pushes wages
  **more than the investment in variable capital dictates** — more than the reserve army
  would justify at the current sV — cannot be represented anywhere in (8.15)–(8.19): at any
  given state of the system, g_e is fully determined. That is exactly the channel the user
  says is missing, and it is missing.
- **(iii) The author has a textual defence in Marx — but it is only half of Marx.** The
  closure follows Marx's *Capital* I, ch. 25, where the wage is treated as the *dependent*
  variable of accumulation (accumulation determines the demand for labour, the reserve army
  regulates the wage, and a wage rise that threatens accumulation is choked off by the
  slackening of accumulation itself). The book even announces the exclusion of the
  Goodwin-type wage-share equation on p. 216: the dual dynamics of employment and the wage
  share "mainly concern the sphere of distribution … and partially reflects, but does not
  capture, the deeper nature of the capitalist mode of production." But Marx's own theory
  leaves room for an independent distributive force: the value of labour-power contains a
  "historical and moral element" (*Capital* I, ch. 6), and distribution is decided in the
  struggle between the classes (the "antinomy, right against right" of *Capital* I, ch. 10).
  A model that wants to *endogenise the waves as struggles over distribution* (as the
  profit-squeeze reading of the book's own 3rd case suggests) cannot collapse g_e into −sV.
- **(iv) Structurally, the problem is deeper than the g_e substitution.** The complete model
  has no state variable for the level of distribution at all: the rate of surplus value e does
  not appear as a state (only its growth rate g_e enters Eq. 8.15). Hence not even a
  *one-off shift* of the wage share — the typical historical product of a struggle (e.g. a
  5-point rise of the wage share in the late 1960s) — has a representation in the system.
  g_e ≅ −sV is the *second* removal of distributional autonomy: first the level of e is dropped
  (Goodwin's wage-share equation is discarded, p. 216), then its growth rate is made a pure
  function of sV.

**Verdict on question 1**: Yes — the critique is substantially correct. The complete model of
Section 8.3 contains no autonomous class-struggle channel: g_e is a deterministic function of
sV alone, and the distributional level is not a state of the system. Two corrections of detail:
the "profit-squeeze case" is in §8.2.1.2 (p. 212), and a *mechanical* reserve-army wage push
does survive inside g_e ≅ −sV; what is eliminated is the autonomous component — the struggle
that aims higher (or lower) than the labour-market condition dictates.

---

## Part II — Formalisation: putting class struggle back into the equation

### II.1 The modified equation

Glombowski's form g_e = ε₀ − ε₁sV already contains an **intercept** ε₀ (the drift of the rate
of surplus value that is independent of the labour market). The book sets (ε₀, ε₁) = (0, 1).
We restore an intercept driven by the class struggle:

**g_e = −sV − CS** (8.15-CS)

with **CS = 0** reproducing the book exactly. Only Eq. (8.15) changes:

**ṙ = −a₁·(sC − sV)·r² + a₂·(δ + τ − sV − CS)·r** (8.15-CS)

CS is measured per year, in the units of g_e and sV.

- **CS > 0**: the class struggle runs *above* the "normal" level dictated by economic
  conditions — the labour movement presses wages harder than the reserve army alone would
  warrant. At any given sV the rate of surplus value grows more slowly (or falls faster);
  in slumps (sV < 0) the wage fall is resisted (real-wage rigidity / union resistance).
  Formally this is a **negative intercept**: even with a stationary labour demand (sV = 0),
  e does not grow — wages press against surplus value.
- **CS < 0**: the struggle runs *below* the normal level — a demobilised working class, or a
  successful capital offensive. Wages lag the labour-market norm in every phase. Formally the
  intercept is positive: the rate of surplus value grows even when the labour market is not
  tight, i.e. wages are held **below the value of labour-power** — in Marx's own list of the
  *counteracting causes* of the FROP (*Capital* III, ch. 14: "depression of wages below the
  value of labour-power").

Both directions are therefore meaningful within Marxist economics: CS > 0 instantiates the
"historical and moral element" of the value of labour-power (workers push the wage baseline
up), CS < 0 instantiates a counteracting cause of the FROP (capital pushes the wage baseline
down).

### II.2 Parameters, reference state and measurement protocol

Book parameters: a₁ = 0.04, a₂ = 0.01, b₀ = 0.004, b₁ = 0.005, b₂ = 0.05, n = 0.015,
g_w = 0.03. Reference state: the book's own closed-form equilibrium (8.25),
{r*, sV*, sC*, δ*, τ*} = {0.1756, 0.1650, 0.6263, 0.05, 0.045}; all runs start from
1.05 × this state (the repository convention). The printed closed form does not involve g_e,
so the same reference state is used for every CS in order to isolate the pure CS effect.

Known caveats of the model, carried over from the repository's earlier audit (F1–F6): the
reference point is **not a true fixed point** of (8.15)–(8.19) — no positive fixed point
exists for the printed parameters — and the nonlinear solution is locally unstable
(dominant eigenvalue +0.0212 ± 0.1505i), with amplitudes growing ~×3 per cycle. The model
therefore leaves its economically meaningful domain (sV < 0; or a share sV, sC exceeding the
whole surplus value, i.e. the book's premise sC + sV ≤ 1) after roughly 85–130 years, i.e.
after about two cycles. All cycle statistics below are computed from the events that
**complete before the first sV < 0 crossing** (t_dom, the repository convention); the first
sC > 1 crossing (t_dom_shares) is reported as well. Script validity: at CS = 0 the
implementation reproduces the repository's model bit-for-bit (identical RHS; Jacobian
matches the analytical sympy Jacobian to 10⁻¹⁰; eigenvalue period 41.7 y at the reference
state, identical to the repository's published 41.7 y).

---

## Part III — Numerical results: what different CS do to the cycle

### III.1 Scenario table (11 values of CS, book parameters, n = 0.015)

sC/sV statistics refer to the completed cycles of sC(t) / sV(t) ending before t_dom
(first sV < 0); "cycles" lists the individual completed peak-to-peak (or trough-to-trough)
spacings in years; rise = median trough→peak (upswing), fall = median peak→trough
(downswing); r_mean is the mean of r over [0, 120 y] (common window for all scenarios);
r_end = r at t_dom.

| CS | t_dom (sV<0) | t_dom (share>1) | T_lin (eigen.) | sC cycles (y) | sC T_med (y) | sC rise (y) | sC fall (y) | sC rise share | sV T_med (y) | sV rise share | r_mean [0,120] | r_end |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| −0.100 | 126 | 87 | 41.7 | [42.2, 45.0] | 43.6 | 23.7 | 19.9 | 54% | 44.4 | 50% | 0.1498 | 0.1146 |
| −0.050 | 129 | 87 | 41.7 | [42.4, 45.8] | 44.1 | 24.4 | 19.7 | 55% | 45.1 | 50% | 0.1454 | 0.1050 |
| −0.020 | 130 | 87 | 41.7 | [42.6, 46.3] | 44.4 | 24.9 | 19.5 | 56% | 45.5 | 51% | 0.1428 | 0.0993 |
| −0.010 | 131 | 87 | 41.7 | [42.6, 46.4] | 44.5 | 25.1 | 19.5 | 56% | 45.6 | 51% | 0.1419 | 0.0975 |
| −0.005 | 131 | 86 | 41.7 | [42.6, 46.5] | 44.6 | 25.1 | 19.4 | 56% | 45.7 | 51% | 0.1415 | 0.0966 |
| **0** | 131 | 86 | 41.7 | [42.7, 46.6] | 44.6 | 25.2 | 19.4 | 57% | 45.8 | 51% | 0.1411 | 0.0956 |
| +0.005 | 132 | 86 | 41.8 | [42.7, 46.7] | 44.7 | 25.3 | 19.4 | 57% | 45.8 | 51% | 0.1407 | 0.0947 |
| +0.010 | 132 | 86 | 41.8 | [42.7, 46.8] | 44.7 | 25.4 | 19.3 | 57% | 45.9 | 51% | 0.1403 | 0.0938 |
| +0.020 | 133 | 86 | 41.8 | [42.8, 46.9] | 44.8 | 25.6 | 19.3 | 57% | 46.1 | 51% | 0.1394 | 0.0920 |
| +0.050 | 135 | 86 | 41.8 | [42.9, 47.4] | 45.2 | 26.1 | 19.1 | 58% | 46.5 | 51% | 0.1370 | 0.0865 |
| +0.100 | 202 | 85 | 41.8 | [43.2, 48.4, 56.4] | 48.4 | 28.3 | 20.1 | 59% | 50.1 | 52% | 0.1331 | 0.0416 |

(All quantitative statements below are robust to the exact cycle-counting protocol: the same
monotone ordering is obtained when the measurement window is the full 6·T window of the
repository convention instead of the t_dom-restricted window.)

### III.2 Findings on the cycle length (job question 2, first part)

1. **The linearised nominal period is essentially CS-invariant**: T_lin = 41.7–41.8 y for the
   whole grid. Reason: CS shifts only the (1,1) entry of the Jacobian at the reference state
   by −a₂·CS ≈ −10⁻⁴ per 0.01 of CS, which barely moves the dominant complex pair
   (+0.0212 ± 0.1505i). The *local clock* of the model is set by the OCC/mechanisation
   dynamics, not by the wage-push channel.
2. **The nonlinear cycle length rises with militancy, by a small amount in the realistic
   band.** First complete cycle: 42.2 y at CS = −0.10 → 42.7 y at CS = 0 → 43.2 y at
   CS = +0.10 (+0.05 y per +0.01 CS). Second cycle: 45.0 → 46.6 → 48.4 y (+0.16–0.18 y per
   +0.01 CS). A sustained, above-normal wage push of 1%/yr therefore lengthens the completed
   long wave by roughly **0.05–0.2 y per 0.01 of CS** (0.1–0.4% of the period); a push of the
   order of 5–10%/yr (CS = +0.05…+0.10) stretches the second cycle to 47.4–48.4 y and the
   third (completed only there) to 56.4 y.
3. **Wage repression shortens the cycle.** At CS = −0.05…−0.10 the first two cycles shrink to
   42.2–45.8 y (the second cycle loses about 1–1.6 y relative to CS = 0); the local
   eigenvalue period stays 41.7 y. The shortening saturates: even at CS = −0.2…−0.3 the
   cycles do not fall below ~41–44 y within the economically meaningful window.
4. **The extreme-militancy regime is qualitatively different**: at CS = +0.2/+0.3 the sC
   spacings read [43.8, 50.4, 61.8] / [44.4, 52.8, 73.6] y — cycles keep stretching and the
   profit rate crawls towards ~1–4% (r_mean 0.08/0.07, r_min 0.014/0.008) — a permanent
   squeeze in which accumulation and employment stall. (Caveat: at such CS the sC + sV ≤ 1
   premise is already violated from ~48–85 y, so the stretched cycles are recorded partly
   outside the economically meaningful domain.)

### III.3 Findings on the upswing and downswing (job question 2, second part)

The wave is **asymmetric by construction** — at CS = 0 the median sC upswing (trough→peak)
lasts 25.2 y against a downswing (peak→trough) of 19.4 y, a 57/43 split (the repository's
historical validation used exactly this property: measured rise shares of 55% on waves 2–3
match the model).

1. **Higher class struggle lengthens the upswing and shortens the downswing.** Across the
   grid the median sC upswing rises from 23.7 y (CS = −0.10) to 28.3 y (CS = +0.10)
   (+0.2–0.3 y per +0.01 CS), the downswing stays flat or shortens slightly
   (19.9 → 19.4 → 19.1 y up to CS = +0.05), and the rise share rises from 54% to 59%. In
   words: with a militant working class the economy spends a *larger share of each wave in
   the phase of expansion of mechanisation (sC rising)*, but does so on a permanently
   squeezed profit rate.
2. **Wage repression shifts the profile the other way**: rise share 54–56% at CS = −0.05…−0.10
   and, in the extreme runs, down to ~51–52% at CS = −0.2…−0.3, with cycle lengths of only
   ~41–44 y. The extreme-repression regime approaches the book's own qualitative "2nd case"
   of Section 8.2.1.2 (g_e ≥ 0): profitability is sustained or rising during accumulation,
   the boom does not generate its own squeeze, and the wave becomes shorter and more
   symmetric.
3. **sV cycles are more symmetric than sC cycles** (rise share ~50–52% at all CS in the
   grid): the employment/hiring share spends roughly equal time rising and falling even
   though its *amplitude* falls with CS (see III.4).

### III.4 Findings on profitability (the level effect — the strongest response)

The dominant and most robust effect of CS is on the **level and trend of the profit rate**,
not on the clock:

- r(60) falls from 0.1485 (CS = −0.10) to 0.1309 (CS = +0.10), i.e. about **−0.0009 per
  +0.01 CS** at the 60-year endpoint (≈ −0.6% of the level per +0.01 CS); mean r over
  [0,120 y] falls from 0.1498 to 0.1331 (≈ −0.0008 per +0.01 CS, ≈ −0.6% of the level).
- An exact analytical anchor: Eq. (8.15-CS) contains the term −a₂·CS·r, so the
  class-struggle component of the profit-rate growth is exactly
  **d ln r/dt = −a₂·CS = −0.01·CS per year** (−0.01% of the level per year at CS = +0.01).
  The observed endpoints are consistent with this direct channel — the 60-year suppression
  (≈ −0.6% per +0.01 CS) equals the direct cumulative effect, and the 120-year suppression
  (≈ −1.5% per +0.01 CS) exceeds it by at most ~×1.2 from the feedback of sV/sC/δ/τ.
  Concretely: an above-normal wage push of CS = +0.01 (= 1% per year) sustained for a decade
  costs the profit rate ≈ 0.1% of its level (≈ 1% over a century); CS = +0.05…+0.10 costs
  ≈ 0.5–1% of the level per decade — the difference between a mild and a severe profit
  squeeze. These small *rates* compound into the level gaps of the table because they act on
  every year of the trajectory.
- Consistently, the first sV < 0 crossing is delayed by militancy (126 y at CS = −0.10, 131 y
  at CS = 0, 202 y at CS = +0.10): with a squeezed profit rate the accumulation boom is
  weaker, so the employment share never overshoots into negative territory as early — the
  economy instead sinks into the low-profit crawl described in III.2.4.

### III.5 The "historical episode" experiment (pulse)

Conscious struggle typically arrives as an *episode*, not as a permanent regime. Figure 3
shows a stylised pulse: CS = +0.02 during t = 18–30 (a "wage explosion" at the boom peak,
1968–73-style), CS = +0.005 during t = 30–44 (the persistent high-wage-share phase), then
CS = −0.015 for t > 44 (the capital counter-offensive, post-1980-style), against the CS = 0
baseline. Because the model's wage-push channel carries the small weight a₂ = 0.01, even this
episode moves the profit rate by only −0.31% (maximum squeeze, at t = 44) to +1.31%
(maximum excess of the offensive, at the end of the window); sV and sC deviate by at most
0.015–0.017 in absolute terms. This *quantifies* how little room the book's calibration
leaves for distributional
struggle to steer the wave — see Part IV for the interpretation and Part V for the limits
this places on the exercise.

---

## Part IV — What the results mean historically

Model time is not calendar time (the repository does not calibrate the model to the data;
cycle lengths come out at 40–50 y). The mapping below therefore compares *directions and
relative magnitudes* with documented historical episodes, using the wage share and the
profit rate as the bridge variables (both are what the model's CS shifts).

### CS ≈ +0.005…+0.02 — "above-normal" struggle: the late-1960s / 1970s analogy

Model behaviour: profit-rate level squeezed (the CS component of d ln r/dt is −a₂·CS, i.e.
roughly −0.05% to −0.2% of the level per decade at CS = +0.005…+0.02, or −0.25% to −1% over
a typical half-century), waves slightly longer and more upswing-heavy, employment-share
troughs deeper during the squeeze, recovery delayed.

Historical counterpart — the wage explosion of 1966–1975: full employment in the OECD core,
union militancy at its postwar peak (French May 1968; Italian *autunno caldo* 1969; UK miners
1972 and 1974, engineers 1971–72; German wildcat strikes 1969–73), real-wage growth running
above productivity growth, and the documented **profit squeeze** of the late boom
(Glyn & Sutcliffe 1972 for the UK; Boddy & Crotty 1975 for the US; Armstrong, Glyn &
Harrison 1991). The US private wage share rose from the mid-1960s to the mid-1970s and the
measured pre-tax profit rate fell by roughly a third between 1965 and 1982
(Duménil & Lévy 1993), the downswing of the fourth Kondratiev wave (1966–1982 in the
repository's periodisation, Table 5.1 of the book). The model reproduces the *direction*: a
sustained positive CS keeps the profit rate below the CS = 0 path and prolongs the phase of
weak accumulation; it does not reproduce the *magnitude* of the 1970s fall, because in the
book's calibration that fall is driven mainly by the OCC channel (a₁ = 0.04) — the FROP —
with the wage push as a modulator (see the caveats in Part V).

### CS ≈ 0 — the "normal" struggle level: the pure reserve-army cycle

Model behaviour: the book's own wave — upswing 25 y, downswing 19 y (rise share 57%), cycle
length 43–47 y, profit rate falling secularly through the OCC channel.

Historical counterpart: periods in which wage movements were dominated by the state of the
labour market, with organised labour too weak (or too integrated) to shift the wage baseline:
the 19th-century British and US industrial cycles (weak craft unions, free immigration, e.g.
the US waves of 1845–1896), and — on the modern side — the "Fordist" decades 1950–1965 in
the US, where real wages rose with productivity and the wage share was broadly trendless
until the mid-1960s. This is the benchmark the book itself describes on p. 216 (the IRAL
"keeps the wage share constant or even decreasing").

### CS ≈ −0.005…−0.02 — "below-normal" struggle: the post-1980 capital offensive

Model behaviour: profits kept high (r_mean up to +0.008–0.01 above the CS = 0 case at
CS = −0.1), wage-bill growth below the reserve-army norm, waves a little shorter, upswings a
little shorter relative to downswings (rise share down to 54–55%), i.e. accumulation runs
"cheaper" for capital.

Historical counterpart — 1979–2000: the monetarist shock and mass unemployment after
1979–82; the symbolic defeats of organised labour (PATCO, 1981; UK miners, 1984–85; German
and Japanese wage moderation from the mid-1980s); US private-sector union density falling
from ~24% (1973) to ~10% (2010s), UK density from ~50% (1980) to ~23%; real wages flat for
the US production worker from 1973 to the mid-1990s while productivity kept rising; the wage
share of US national income falling from the mid-1970s/early 1980s back to and below its
1965 level. Profits recovered correspondingly: the US profit rate rose from its 1982 trough
through the 1990s (Duménil & Lévy), and the period 1982–2007 is the upswing of the fifth
Kondratiev wave in the repository's periodisation. Marx's *Capital* III, ch. 14
("depression of wages below the value of labour-power") names exactly this mechanism as a
counteracting cause of the FROP — the model's CS < 0 is its formalisation, and the model
confirms the expected effect: the fall of the profit rate is slowed.

### CS ≈ −0.1…−0.3 — extreme repression: the 1920s analogy

Model behaviour: profitability sustained near its initial level for more than a century of
model time (r_mean 0.15–0.17 at CS = −0.2…−0.3, almost no squeeze), cycles shortened to
~41–44 y and nearly symmetric (rise share ~51%).

Historical counterpart — the US 1920s under the "American Plan" open-shop offensive: real
manufacturing wages roughly stagnant between 1923 and 1929 while output per worker rose by
about a quarter, pushing the wage share down and the profit share to a peak by 1929 — the
classic "Roaring Twenties" configuration of high and rising profitability with wage
repression (and, as the decade's end shows, one that stored up the overaccumulation crisis —
which this five-equation model, having no credit/finance layer, cannot itself generate).
The earlier US episode 1898–1907 (post-1896 wage-share trough, strike defeats of the
mid-1890s, then rapid accumulation) is a weaker, partial analogue.

### CS ≳ +0.05 — regime-scale militancy: the squeeze that does not end

Model behaviour: profit rate collapses towards 1–8%, cycles stretch to 50–74 y, and the
economy slides into the low-profit crawl (the model leaves its economically meaningful
domain early, so these runs are boundary illustrations rather than forecasts).

Historical counterpart: not a long historical epoch (no capitalist economy sustained a
5–15%/yr excess wage push for long — the resistance of capital, disinvestment and crisis
have always cut it short), but the *episodes* that point that way: France 1936–1938 (Popular
Front wage rises and the 40-hour week, the profit collapse of 1937 and the "capital strike"
of 1937–38, ended by the wage pause of November 1938), the UK 1974–1979 (wage-push under
full employment, the IMF crisis of 1976, and the squeeze that prepared the 1979–80 reversal),
and the US 1970–1974 (wage controls followed by the 1973–75 recession). These episodes are
the historical evidence that regime-scale CS > 0 is not sustainable *as a constant*: the
system responds with crisis and, ultimately, a capital offensive that resets CS below zero —
exactly the sequence stylised in the pulse experiment of Section III.5 and figure 3.

---

## Part V — Conclusions and limits

1. **The critique is right.** Sections 8.2.1.2–8.2.1.3 (book pp. 207–215) acknowledge that the
   sign and intensity of g_e carry the class struggle, and the 3rd case (g_e < 0,
   profit-squeeze, p. 212) makes the wage push the reason behind the FROP — but the closing
   relation g_e ≅ −sV (p. 214) reduces the wage dynamic to the reserve army, and the complete
   five-equation model of Section 8.3 contains no variable through which a struggle "aiming
   higher than what the investment in variable capital dictates" could act. The book's own
   later claims that the rate of surplus value reflects "class struggle and the like"
   (p. 238) and "the level of class struggle in the sphere of circulation" (p. 239) are
   therefore statements about the *concept* of the model, not about its equations.
2. **The CS-extension** (g_e = −sV − CS; only Eq. 8.15 changes) puts the struggle back in a
   minimal, well-defined way: CS > 0 is an above-normal wage push (the "historical and moral
   element" of the value of labour-power), CS < 0 a below-normal push (Marx's counteracting
   cause, the depression of wages below value). CS = 0 reproduces the book exactly.
3. **What CS does to the wave** (the quantitative answer to the job question):
   - the *local (linearised) period is unchanged* (41.7–41.8 y) — CS shifts only the weak
     (1,1) Jacobian entry; the model's clock is set by the OCC mechanism;
   - the *nonlinear* cycle lengthens mildly with militancy (≈ +0.05 y per +0.01 CS on the
     first and ≈ +0.16–0.18 y per +0.01 CS on the second completed cycle; stretching strongly
     only at regime-scale CS ≥ +0.05), and shortens mildly under repression;
   - the *upswing* lengthens and the *downswing* shortens with CS (rise share of the sC wave:
     54% at CS = −0.10 → 57% at CS = 0 → 59% at CS = +0.10);
   - the *level effect dominates*: the CS component of d ln r/dt is −a₂·CS = −0.01·CS
     exactly, so a sustained above-normal struggle of 1%/yr (CS = +0.01) costs ≈ 0.1% of the
     profit rate per decade (≈ 1% over a century), and regime-scale pushes
     (CS ≥ +0.05…+0.10) drive the profit rate towards zero over the simulation horizon
     (stagnation); wage repression keeps profitability high and slows the FROP.
4. **Historical reading** (directions, not point calibration): CS > 0 ↔ the 1966–1975 wage
   explosion and the profit squeeze of the fourth-wave downswing (1966/73–1982); CS ≈ 0 ↔ the
   reserve-army-dominated 19th-century waves and the 1950–1965 Fordist benchmark; CS < 0 ↔ the
   post-1980 offensive (Thatcher/Reagan, union-density collapse, wage-share decline, profit
   recovery, the fifth-wave upswing of 1982–2007); CS ≪ 0 ↔ the 1920s "open-shop" decade;
   regime-scale CS > 0 is unsustainable and historically appears only as bounded episodes
   (1936–38 France, 1974–79 UK), each terminated by crisis and a subsequent offensive.
5. **Limits of the exercise** (honest statement):
   - the direct quantitative channel of the wage push is small *by the book's own
     calibration*: g_e enters only the "counteracting forces" bracket of Eq. (8.15), which
     the authors deliberately weight with a₂ = 0.01 against a₁ = 0.04 for the OCC channel
     (book p. 215: the smaller effect of the counteracting forces is "the transient nature"
     of these forces). Inserting CS therefore changes the profit-rate *level* first and the
     *clock* second — a finding in its own right, since it shows that even a formal admission
     of the class struggle would leave the book's model largely OCC-determined;
   - the level of distribution is not a state variable of the model (the Goodwin wage-share
     equation was dropped, p. 216), so a *level* shift of the wage share — historically the
     most visible product of class struggle — cannot be represented, only the growth-rate
     premium CS. A faithful extension would add the rate of surplus value e (or the wage
     share) as a sixth state with a struggle-augmented adjustment equation; that is the
     natural next step and would make the historical magnitudes (wage-share swings of ±4–7
     points) commensurable with the model;
   - the model has no positive fixed point at the printed parameters and leaves its
     economically meaningful domain after ~85–130 y (findings F1–F6 of this repository), so
     the quantitative statements above refer to the first two completed cycles;
   - model time is not calibrated to calendar time; the historical mapping is qualitative.

---

## References

- Armstrong, P., Glyn, A. & Harrison, J. (1991). *Capitalism since 1945*. Blackwell.
- Boddy, R. & Crotty, J. (1975). Class conflict and macro-policy: the political business
  cycle. *Review of Radical Political Economics*, 7(1), 1–19.
- Chatzarakis, N., Tsaliki, P. & Tsoulfidis, L. (2022). *Economic Growth and Long Cycles: A
  Classical Political Economy Approach*. Routledge. Sections 8.2.1.2 (pp. 207–212),
  8.2.1.3 (pp. 213–215), 8.2.2 (p. 216), 8.3 (pp. 221–235), 8.4 (pp. 236–242).
- Duménil, G. & Lévy, D. (1993). *The Economics of the Profit Rate*. Edward Elgar.
- Glombowski, J. (1983). A Marxian model of long run capitalist development.
  *Zeitschrift für Nationalökonomie*, 43(4), 363–382.
- Glyn, A. & Sutcliffe, B. (1972). *British Capitalism, Workers and the Profit Squeeze*.
  Penguin.
- Goodwin, R. M. (1967). A growth cycle. In C. H. Feinstein (ed.), *Socialism, Capitalism
  and Economic Growth*. Cambridge University Press.
- Marx, K. *Capital*, Vol. I (ch. 6: the value of labour-power; ch. 10: the struggle over
  the working day; ch. 25: the general law of capitalist accumulation); Vol. III, ch. 14
  (the counteracting causes of the FROP).
- Repository reports: `docs/model_validation_n_scenarios.md`, `docs/sv_sign_flip_report.md`,
  `docs/Technical_Report_Ch8_Model_Errata.docx` (findings F1–F6).
