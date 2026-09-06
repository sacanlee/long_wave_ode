# Class Struggle in the Long-Wave Model: The Closure g_e ≅ −sV and an Extension with an Autonomous Wage-Push Term (CS)

**A research note on Eqs. (8.15)–(8.19) of Chatzarakis, Tsaliki and Tsoulfidis (2022)**

Prepared by: sacanlee | 6 September 2026

Companion material (fully reproducible): `scripts/class_struggle_ode.py`;
`data/class_struggle_results.json`; `figures/class_struggle_trajectories.png`,
`figures/class_struggle_cycle_metrics.png`, `figures/class_struggle_pulse_episode.png`.

**Reference**: Chatzarakis, N., Tsaliki, P. & Tsoulfidis, L. (2022), *Economic Growth and
Long Cycles: A Classical Political Economy Approach*, Routledge, Ch. 8. Book page numbers
below refer to this edition. The model equations (8.15)–(8.19) and Section 8.3.3 have been
reproduced and audited numerically in the companion repository (see
`docs/model_validation_n_scenarios.md`, `docs/sv_sign_flip_report.md` and
`docs/Technical_Report_Ch8_Model_Errata.docx` for the internal-consistency findings F1–F6 of
the printed Chapter 8 material).

---

## Abstract

Chapter 8 of the book builds a five-equation model of the long cycle around the profit rate
and the shares of surplus value invested in constant and variable capital. In the qualitative
analysis of Section 8.2.1.2 the growth rate of the rate of surplus value, g_e, carries the
distributional content of the system: its sign is said to be related to "relations of
production (income distribution, class struggle, exploitation of labor)", and its negative
case is presented as the profit-squeeze case in which "the intensity of the class struggle
leading to wage squeeze of profits is the reason behind the Law of the FROP". In Section
8.2.1.3, however, the model is closed by setting g_e ≅ −sV, so that in the complete model
(8.15)–(8.19) the wage dynamic is fully determined by the investment in variable capital
(the reserve-army mechanism), and no variable remains through which an *autonomous* wage push
— a struggle that aims higher (or lower) than the labour-market condition dictates — could
operate. This note examines that observation, proposes a minimal formalisation of an
above/below-normal class struggle by adding a parameter CS (g_e = −sV − CS, with CS = 0
reproducing the book), and solves the modified system numerically for a grid of CS values,
focusing on the length of the cycle and on the lengths of its upswing and downswing. The
main results are: (i) the linearised period is essentially CS-invariant (41.7–41.8 y), but
the nonlinear cycles lengthen with CS (≈ +0.05 y per +0.01 CS on the first, ≈ +0.16–0.18 y
per +0.01 CS on the second completed cycle); (ii) the upswing (trough→peak of sC) lengthens
and the downswing stays flat or shortens slightly, so the rise share of the wave increases
from 54% (CS = −0.10) to 57% (CS = 0) to 59% (CS = +0.10); (iii) the dominant effect of CS
is on the *level* of the profit rate, since Eq. (8.15-CS) contains the exact term
−a₂·CS·r. An analytical mechanism is derived: with the profit rate drifting secularly
downwards, the share dynamics rotate around the receding balance point (sV̄, sC̄) = (τ/r,
A/r), A = δ + τ + n, with the clock ω = √(A·r); CS depresses r, slows this rotation, and the
rising centre biases every wave towards its upswing. The different CS levels are then
confronted with documented historical episodes (the 1966–75 wage explosion and the
fourth-wave profit squeeze; the post-1980 capital offensive; the 1920s; regime-scale
episodes such as France 1936–38 and the UK 1974–79). The note closes with the limits of the
exercise and a suggestion for a fuller extension in which the level of distribution enters
as a sixth state variable.

**Keywords**: long waves; rate of surplus value; class struggle; profit squeeze; industrial
reserve army; Goodwin-type cycles; Kondratiev waves

---

## 1. Introduction

The long-wave model of Chapter 8 of *Economic Growth and Long Cycles* (Chatzarakis, Tsaliki
& Tsoulfidis 2022) describes a capitalist economy with five state variables — the rate of
profit r, the shares of surplus value invested in variable (sV) and constant (sC) capital,
the capital-devaluation rate δ and the rate of technical change τ — through the equations
(book pp. 221–235)

```
(8.15)  r'  = −a1·(sC − sV)·r² + a2·(δ + τ − sV)·r
(8.16)  sV' = sC·r − δ − τ − n
(8.17)  sC' = sC·(τ − sV·r)
(8.18)  δ'  = b0 + δ·(gw − sC·r)
(8.19)  τ'  = b1·Erf(δ − b2)
```

with a₁ = 0.04, a₂ = 0.01, b₀ = 0.004, b₁ = 0.005, b₂ = 0.05, n = 0.015, g_w = 0.03.
Equation (8.15) is derived in Section 8.2.1.3 from a growth-rate equation that contains the
growth rate of the rate of surplus value, g_e, additively, and the model is closed by
substituting g_e ≅ −sV.

A comparison between the qualitative analysis of Sections 8.2.1.1–8.2.1.2 and the closure
adopted in Sections 8.2.1.3–8.3 gives rise to a natural question: the sign and the intensity
of g_e are assigned, in the qualitative part, to the sphere of class relations — and the
negative case of g_e is presented as the profit-squeeze case, i.e. as an active push for
higher wages — yet in the complete five-equation model g_e is made a deterministic function
of the investment in variable capital. It is then unclear whether a *conscious* class
struggle — a wage push aiming higher than what the investment in variable capital dictates —
can be represented anywhere in the model at all.

The purpose of this note is fivefold: (i) to assess that question against the text of the
book; (ii) to propose a minimal formalisation of an above/below-normal class struggle, by
introducing a parameter CS into the g_e relation; (iii) to report the numerical consequences
of different CS values, with particular attention to the length of the cycle and to the
lengths of its upward and downward phases; (iv) to give the analytical mechanism behind
those consequences; and (v) to confront the different CS levels with documented historical
episodes. A brief discussion of the limits of the exercise and of a possible fuller
extension closes the note.

Throughout, model time is not calibrated to calendar time; cycle lengths of 40–50 years
emerge from the parameters of the book. Historical statements therefore refer to directions
and relative magnitudes, not to point calibrations.

## 2. The role of g_e in the model, and the closure adopted

### 2.1 The qualitative role of g_e (Sections 8.2.1.1–8.2.1.2, pp. 205–212)

In the two-variable (OCC, profit-rate) system analysed in Sections 8.2.1.1–8.2.1.2 the
growth rate of the rate of surplus value, g_e, is treated as an exogenous parameter whose
value selects the qualitative regime. Three cases are distinguished (book pp. 207–212):

- the 1st case, g_e = 0: the distribution of income remains constant;
- the 2nd case, g_e ≥ 0: the rate of surplus value grows faster than variable capital, and
  profitability is sustained or rises;
- the **3rd case, g_e < 0**: "the surplus value is reduced relative to variable capital
  (the case of profit-squeeze)" (p. 212) — an active push for higher wages.

The book is explicit that the *sign* of g_e is where class relations enter the system.
Already on p. 207:

> "The bifurcation of equilibrium point B is, therefore, defined by the sign of g_e, which is
> a variable related not only to productive forces of the system (i.e., technological change,
> capital accumulation, etc.) but also to relations of production (i.e., income distribution,
> class struggle, exploitation of labor, etc.)."

and, discussing the meaning of the g_e < 0 case for the tendential fall in the rate of
profit (FROP), p. 212:

> "The significance of point B is only theoretical and may be interpreted in favor of the
> 'profit squeeze' argument, according to which **the intensity of the class struggle leading
> to wage squeeze of profits is the reason behind the Law of the FROP**. Under these
> circumstances, the Law of FROP cannot be viewed as a general economic law independent of
> people's 'will' and capable of realistically describing the dynamics of a capitalist
> economy."

(For completeness: this third case appears in Section 8.2.1.2, which immediately precedes
Section 8.2.1.3.)

### 2.2 The closure adopted in the complete model (Section 8.2.1.3, pp. 213–215)

Section 8.2.1.3 derives the profit-rate dynamics on the capital stock. From the growth-rate
equation

```
r'/r = −(sC − sV)·r + δ + τ + g_e                                    (8.14)
```

the book needs a theory of g_e. The one it adopts is purely mechanical (p. 214):

> "According to Marx (*Capital* I, ch. 25), the growth rate of the rate of surplus value,
> g_e, is inversely proportional to the change in the labor force, or, what is the same, to
> the investment in variable capital, sV. Glombowski (1983), for example, hypothesized that
> g_e = ε₀ − ε₁sV, where ε₀ and ε₁ are small positive constants. In our case, we specify a
> somewhat simpler relation according to which g_e ≅ −sV and any remaining differences are
> captured by two small positive parameters a₁ and a₂..."

Substitution of g_e ≅ −sV into (8.14), with the calibration parameters a₁ and a₂, gives the
printed equation

```
r' = −a1·(sC − sV)·r² + a2·(δ + τ − sV)·r                            (8.15)
```

Two features of the resulting complete model (8.15)–(8.19) should be noted:

1. g_e itself appears in **no equation** of the complete model: the state equations
   (8.16)–(8.19) contain only r, sV, sC, δ, τ and constants, and in (8.15) the wage-push
   content of g_e survives only as the term −sV inside the "counteracting forces" bracket,
   weighted by the small parameter a₂ = 0.01 against a₁ = 0.04 for the OCC term. The book
   itself describes the second bracket as the locus of "the transient nature of the
   counteracting forces" (p. 215).
2. The level of distribution is not a state variable of the model: the rate of surplus
   value e does not enter the state vector (only its growth rate g_e enters Eq. 8.15
   through the derivation above), and the book explicitly declines to carry over Goodwin's
   wage-share equation into the model (p. 216: the dual dynamics of employment and the wage
   share "mainly concern the sphere of distribution … and partially reflects, but does not
   capture, the deeper nature of the capitalist mode of production").

### 2.3 Assessment

The question raised in the introduction can now be answered precisely.

1. **What the closure retains.** The relation g_e ≅ −sV is not the abolition of wage
   pressure altogether. Because sV is high when capital accumulation into employment
   exceeds the natural growth of the labour force (Eq. 8.16: sV′ = sC·r − δ − τ − n), the
   rate of surplus value grows slowly, or falls, exactly in the tight-labour phases of the
   boom — the reserve-army mechanism in the Goodwin tradition, described by the book on
   p. 216: "The pressure from an increasing IRAL keeps the wage share constant or even
   decreasing, increases the discipline of the labor force, and leads to a higher level of
   labor intensity." This is a wage push, but one fully dictated by the state of the labour
   market: workers obtain higher wages only when, and to the extent that, the accumulation
   of variable capital dictates.
2. **What the closure loses.** At any given state of the system (any given sV), g_e is fully
   determined. A working class that pushes wages *more* than the investment in variable
   capital dictates — more than the reserve army would justify at the current sV — cannot be
   represented anywhere in (8.15)–(8.19). The same holds, symmetrically, for a wage push
   *below* the labour-market norm (a demobilised working class or a successful capital
   offensive). And because the distribution level is not a state variable, even a one-off
   shift of the wage share — the typical historical product of a struggle, e.g. a 4–7 point
   rise of the wage share in the late 1960s and early 1970s — has no representation in the
   system.
3. **The textual basis and its limits.** The closure follows Marx's *Capital* I, ch. 25,
   where the wage is the *dependent* variable of accumulation (accumulation determines the
   demand for labour; the reserve army regulates the wage; a wage rise that threatens
   accumulation is choked off by the slackening of accumulation itself). Marx's own theory,
   however, also contains an independent distributive element: the value of labour-power
   has a "historical and moral element" (*Capital* I, ch. 6), and distribution is decided in
   the struggle between the classes (*Capital* I, ch. 10). A model that wishes to
   *endogenise the long waves as struggles over distribution* — the reading suggested by the
   book's own 3rd case of Section 8.2.1.2 — cannot collapse g_e into −sV.
4. **An internal tension.** Later parts of the chapter continue to assign the rate of
   surplus value a decisive distributive role *in words*: g_e is called "the key explanatory
   variable of the upswings and downswings of the rate of profit" (pp. 237–238), and the
   rate of surplus value "a key distributive variable reflecting the level of class struggle
   in the sphere of circulation" (p. 239). These statements characterise the *concept* of
   the model; in its equations the wage dynamic is a function of sV alone.

The observation examined here is therefore, in our assessment, essentially correct, with
three qualifications of detail: (i) the profit-squeeze taxonomy (g_e < 0) is presented in
Section 8.2.1.2 (p. 212), immediately before the closure of Section 8.2.1.3; (ii) a
*mechanical* reserve-army wage push does survive inside g_e ≅ −sV — it is the *autonomous*
component, the struggle that aims higher (or lower) than the labour-market condition
dictates, that is eliminated; and (iii) structurally, the elimination operates twice: the
level of distribution is dropped (no wage-share state equation, p. 216), and then the growth
rate of the rate of surplus value is made a pure function of sV.

## 3. A minimal extension: the class-struggle level CS

### 3.1 Specification

Glombowski's form g_e = ε₀ − ε₁·sV (Glombowski 1983) contains an **intercept** ε₀: the
drift of the rate of surplus value that is independent of the labour market. The book sets
(ε₀, ε₁) = (0, 1). The extension restores an intercept driven by the class struggle:

```
g_e = −sV − CS                                                    (8.15-CS)

r'  = −a1·(sC − sV)·r² + a2·(δ + τ − sV − CS)·r                   (8.15-CS)
```

with **CS = 0** reproducing the book exactly; only Eq. (8.15) changes. CS is measured per
year, in the units of g_e and of sV.

- **CS > 0**: the class struggle runs *above* the "normal" level dictated by economic
  conditions — the labour movement presses wages harder than the reserve army alone would
  warrant. At any given sV the rate of surplus value grows more slowly (or falls faster),
  and in slumps (sV < 0) the wage fall is resisted (real-wage rigidity / union resistance).
  Formally this is a **negative intercept**: even with a stationary labour demand (sV = 0)
  the rate of surplus value does not grow — wages press against surplus value. This
  instantiates the "historical and moral element" of the value of labour-power.
- **CS < 0**: the struggle runs *below* the normal level — a demobilised working class, or a
  successful capital offensive. Wages lag the labour-market norm in every phase. Formally
  the intercept is positive: the rate of surplus value grows even when the labour market is
  not tight, i.e. wages are held **below the value of labour-power** — one of Marx's own
  *counteracting causes* of the FROP (*Capital* III, ch. 14: "depression of wages below the
  value of labour-power").

Both directions are thus meaningful within the classical-Marxian framework.

### 3.2 Parameters, reference state and measurement protocol

Book parameters (p. 233): a₁ = 0.04, a₂ = 0.01, b₀ = 0.004, b₁ = 0.005, b₂ = 0.05,
n = 0.015, g_w = 0.03. Reference state: the book's closed-form equilibrium (8.25),
{r*, sV*, sC*, δ*, τ*} = {0.1756, 0.1650, 0.6263, 0.05, 0.045}; all runs start from 1.05 ×
this state (the convention of the companion repository). The printed closed form does not
involve g_e, so the same reference state is used for every CS, isolating the pure CS effect.

Caveats carried over from the repository's audit of Chapter 8 (findings F1–F6): the
reference point is **not a true fixed point** of (8.15)–(8.19) — no positive fixed point
exists for the printed parameters — and the nonlinear solution is locally unstable
(dominant eigenvalue +0.0212 ± 0.1505i at the reference state), with amplitudes growing by a
factor of roughly three per cycle. The model therefore leaves its economically meaningful
domain — sV < 0, or a share sV, sC exceeding the whole surplus value (the book's premise
sC + sV ≤ 1) — after roughly 85–130 years, i.e. after about two cycles. All cycle
statistics below are computed from events that **complete before the first sV < 0 crossing**
(t_dom, the repository convention); the first time a share exceeds 1 (t_dom_shares) is
reported as well. Implementation validity: at CS = 0 the present code reproduces the
repository's implementation of the book bit-for-bit (identical right-hand sides; Jacobian
matches the analytical sympy Jacobian to 10⁻¹⁰; eigenvalue period 41.7 y at the reference
state).

## 4. Numerical results: cycle length and the lengths of the upswing and downswing

### 4.1 Scenario table

Eleven values of CS between −0.10 and +0.10 (book parameters, n = 0.015). In Table 1, the
sC/sV statistics refer to the completed cycles of sC(t) / sV(t) ending before t_dom (first
sV < 0); "cycles" lists the individual completed peak-to-peak (trough-to-trough) spacings in
years; rise = median trough→peak duration (upswing), fall = median peak→trough duration
(downswing); r_mean is the mean of r over the common window [0, 120 y]; r_end = r at t_dom.
All statements of Section 4 are robust to the cycle-counting protocol: the same monotone
ordering is obtained when the measurement window is the full 6·T window of the repository
convention.

**Table 1 — Cycle statistics of sC and sV and profitability over the CS grid.**

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

*Source: author's computations (script `class_struggle_ode.py`).*

### 4.2 Effect on the cycle length

1. **The linearised nominal period is essentially CS-invariant**: T_lin = 41.7–41.8 y over
   the whole grid. CS shifts only the (1,1) entry of the Jacobian at the reference state,
   by −a₂·CS ≈ −10⁻⁴ per 0.01 of CS, which barely moves the dominant complex pair
   (+0.0212 ± 0.1505i). The local clock of the model is set by the OCC/mechanisation
   dynamics, not by the wage-push channel.
2. **The nonlinear cycle length rises with militancy, by a small amount in the realistic
   band.** First complete cycle: 42.2 y at CS = −0.10 → 42.7 y at CS = 0 → 43.2 y at
   CS = +0.10 (+0.05 y per +0.01 CS). Second cycle: 45.0 → 46.6 → 48.4 y (+0.16–0.18 y per
   +0.01 CS). A sustained above-normal wage push of 1%/yr therefore lengthens the completed
   long wave by roughly 0.05–0.2 y per 0.01 of CS (0.1–0.4% of the period); a push of the
   order of 5–10%/yr (CS = +0.05…+0.10) stretches the second cycle to 47.4–48.4 y and the
   third (completed only there) to 56.4 y.
3. **Wage repression shortens the cycle.** At CS = −0.05…−0.10 the first two cycles shrink
   to 42.2–45.8 y (the second cycle loses about 1–1.6 y relative to CS = 0), while the
   local eigenvalue period stays 41.7 y. The shortening saturates: even at CS = −0.2…−0.3
   the cycles do not fall below ~41–44 y within the economically meaningful window.
4. **The extreme-militancy regime is qualitatively different**: at CS = +0.2/+0.3 the sC
   spacings read [43.8, 50.4, 61.8] / [44.4, 52.8, 73.6] y — the cycles keep stretching,
   and the profit rate crawls towards 1–4% (r_mean 0.081/0.073; r_min 0.014/0.008) — a
   permanent squeeze in which accumulation and employment stall. (At such CS the premise
   sC + sV ≤ 1 is violated from ~48–85 y onward, so the stretched cycles are recorded
   partly outside the economically meaningful domain.)

### 4.3 Effect on the lengths of the upswing and the downswing

The wave is asymmetric by construction: at CS = 0 the median sC upswing (trough→peak)
lasts 25.2 y against a downswing (peak→trough) of 19.4 y — a 57/43 split. (The repository's
historical validation uses exactly this property: measured rise shares of ~55% on waves 2–3
of the book's Table 5.1 match the model.)

1. **Higher class struggle lengthens the upswing and leaves the downswing flat or slightly
   shorter.** Across the grid the median sC upswing rises from 23.7 y (CS = −0.10) to 28.3 y
   (CS = +0.10) (+0.2–0.3 y per +0.01 CS at moderate and high CS); the downswing stays flat
   or shortens slightly (19.9 → 19.4 → 19.1 y up to CS = +0.05); the rise share rises from
   54% to 57% to 59%. In words: with a militant working class the economy spends a *larger
   share of each wave in the expansion phase of mechanisation (sC rising)*, but does so on
   a permanently squeezed profit rate.
2. **Wage repression shifts the profile the other way**: rise share 54–56% at
   CS = −0.05…−0.10 and, in the extreme runs, down to ~51–52% at CS = −0.2…−0.3, with cycle
   lengths of only ~41–44 y. The extreme-repression regime approaches the qualitative
   "2nd case" of Section 8.2.1.2 (g_e ≥ 0): profitability is sustained or rising during
   accumulation, the boom does not generate its own squeeze, and the wave becomes shorter
   and more symmetric.
3. **sV cycles are more symmetric than sC cycles** (rise share ~50–52% at all CS in the
   grid): the employment/hiring share spends roughly equal time rising and falling even
   though its amplitude responds to CS (Section 4.4).

### 4.4 Effect on profitability (the level effect — the strongest response)

The dominant and most robust effect of CS is on the **level and trend of the profit rate**,
not on the clock:

- r(60) falls from 0.1485 (CS = −0.10) to 0.1309 (CS = +0.10), i.e. about −0.0009 per
  +0.01 CS at the 60-year endpoint (≈ −0.6% of the level per +0.01 CS); the mean of r over
  [0, 120 y] falls from 0.1498 to 0.1331 (≈ −0.0008 per +0.01 CS).
- Exact analytical anchor: Eq. (8.15-CS) contains the term −a₂·CS·r, so the class-struggle
  component of the profit-rate growth is exactly **d ln r/dt = −a₂·CS = −0.01·CS per year**
  (−0.01% of the level per year at CS = +0.01). The observed endpoints are consistent with
  this direct channel — the 60-year suppression (≈ −0.6% per +0.01 CS) equals the direct
  cumulative effect, and the 120-year suppression (≈ −1.5% per +0.01 CS) exceeds it by at
  most ~×1.2 through the feedback of sV/sC/δ/τ. Concretely: an above-normal wage push of
  CS = +0.01 (= 1% per year), sustained for a decade, costs the profit rate ≈ 0.1% of its
  level (≈ 1% over a century); CS = +0.05…+0.10 costs ≈ 0.5–1% of the level per decade —
  the difference between a mild and a severe profit squeeze. These small rates compound
  into the level gaps of Table 1 because they act on every year of the trajectory.
- Consistently, the first sV < 0 crossing is delayed by militancy (126 y at CS = −0.10,
  131 y at CS = 0, 202 y at CS = +0.10): with a squeezed profit rate the accumulation boom
  is weaker, so the employment share never overshoots into negative territory as early —
  the economy instead sinks into the low-profit crawl described in Section 4.2 (item 4).

### 4.5 A "historical episode" experiment (pulse)

Conscious struggle typically arrives as an *episode*, not as a permanent regime. As an
illustration, Figure 3 shows a stylised pulse: CS = +0.02 during t = 18–30 (a "wage
explosion" at the boom peak, 1968–73-style), CS = +0.005 during t = 30–44 (a persistent
high-wage-share phase), then CS = −0.015 for t > 44 (a capital counter-offensive,
post-1980-style), against the CS = 0 baseline. Because the model's wage-push channel carries
the small weight a₂ = 0.01, even this episode moves the profit rate by only −0.31% (maximum
squeeze, at t = 44) to +1.31% (maximum excess of the offensive, at the end of the window);
sV and sC deviate by at most 0.015–0.017 in absolute terms. This quantifies how little room
the calibration of the book leaves for distributional struggle to steer the wave (see
Sections 7 and 8 for the interpretation and limits).

**Figure 1** — `class_struggle_trajectories.png`: profit rate r, sV and sC under the CS grid
(each trajectory truncated at its first sV < 0).
**Figure 2** — `class_struggle_cycle_metrics.png`: cycle length, upswing/downswing lengths,
rise share and profitability as functions of CS.
**Figure 3** — `class_struggle_pulse_episode.png`: the stylised wage-explosion → offensive
episode as deviations from the CS = 0 baseline.
**Figure 4** — `class_struggle_ratchet_recoveries.png`: the "pointless recoveries" at a fixed
CS = +0.02 (Section 7.6): the accumulation wave keeps cycling (boom peaks at t ≈ 50, 97, 150,
213 y) while the profit rate ratchets down at every peak (r: 0.146 → 0.115 → 0.083 → 0.051).

## 5. The mechanism in analytical terms: why CS > 0 lengthens the cycle and tilts it toward the upswing

This section derives the mechanism behind Sections 4.2–4.3. All identities below are exact;
the oscillator statements are linearisations around the moving centre of the share
subsystem, verified against the runs.

### 5.1 The exact phase structure: sC is the integral of a "gap"

From (8.17):

```
d ln sC/dt = τ − sV·r
```

so sC rises ⇔ sV·r < τ and falls ⇔ sV·r > τ, and every extremum of sC is a crossing of the
hiring-pressure variable h := sV·r through the threshold τ. In the simulations this is not
an approximation: every detected peak and trough of sC(t) satisfies |sV·r − τ| ≤ 3 × 10⁻⁵.
The threshold itself is quasi-stationary — τ stays within the band [0.0461, 0.0472] (±1%)
over more than a century — because τ′ = b₁·Erf(δ − b₂) integrates only small δ-oscillations
(b₁ = 0.005), and δ is anchored by (8.18), δ′ = b₀ + δ(g_w − sC·r) ≈ 0, at the value that
keeps sC·r ≈ g_w + b₀/b₂ = 0.11. The whole long wave is therefore the story of h crossing a
fixed threshold.

### 5.2 What drives h, and the exact footprint of CS

Combining (8.16) and (8.15-CS):

```
d ln(sV·r)/dt = sV′/sV + r′/r = (sC·r − A)/sV − a1·(sC − sV)·r + a2·(δ + τ − sV) − a2·CS
```

with A := δ + τ + n ≈ 0.11 (quasi-constant, ±2–3%). The first term is the employment side
(hiring accelerates when accumulation sC·r exceeds the "natural" rate A), the bracket is the
profit side, and **CS enters as an exact proportional drag −a₂·CS on the log-slope of the
crossing variable itself** — h contains r, and r′/r contains −a₂·CS. Nothing else in the
system changes.

### 5.3 The clock: rotation around a drifting centre

For a slowly moving r, the two share dynamics have "balance levels" — sC′ = 0 at
sV̄ := τ/r and sV′ = 0 at sC̄ := A/r — both scaling as 1/r. With log-deviations
p := ln(sC/sC̄), q := ln(sV/sV̄) and g := d ln r/dt < 0:

```
ṗ ≈ −τ·q + g,    q̇ ≈ (A·r/τ)·p + g     ⇒     p̈ + A·r·p = −τ·g
```

i.e. a harmonic oscillator whose natural frequency is set by the product of the
"natural-growth" constant and the profit rate:

```
ω = √(A·r) = √((δ + τ + n)·r),     T = 2π/√(A·r)
```

Table 2 verifies this against the runs (cycle means; the absolute level comes out ~10% high
because of the large-amplitude nonlinearity, but the dependence on r is the point).

**Table 2 — Measured cycle spacings vs the oscillator clock 2π/√(A·r).**

| run | cycle | r̄ over the cycle | actual spacing (y) | 2π/√(Ā·r̄) (y) |
|---|---|---:|---:|---:|
| CS = −0.10 | 2 | 0.1434 | 45.0 | 49.6 |
| CS = 0 | 1 | 0.1626 | 42.7 | 46.7 |
| CS = 0 | 2 | 0.1317 | 46.6 | 51.9 |
| CS = +0.10 | 1 | 0.1574 | 43.2 | 47.5 |
| CS = +0.10 | 2 | 0.1203 | 48.4 | 54.3 |

The cross-CS differentials match almost exactly, which is what matters for the comparison:
√(0.1317/0.1203) = 1.046 predicts 48.7 y for the second cycle at CS = +0.10 versus 46.6 y at
CS = 0 (actual: 48.4 y); √(0.1317/0.1434) = 0.958 predicts 44.7 y at CS = −0.10 (actual:
45.0 y). **To first order, the wave period scales as T ∝ r^(−1/2).**

### 5.4 CS lowers r, so it slows the clock

Eq. (8.15-CS) contains −a₂·CS·r, hence d ln r/dt ∋ −a₂·CS exactly. Two consequences are
visible in the runs: r is lower at every date (measured mean d ln r over the first upswing:
−0.0032/yr at CS = −0.10, −0.0043/yr at CS = 0, −0.0054/yr at CS = +0.10 — the extra
−0.001/yr is −a₂·CS plus a small feedback); and, because the δ-balance keeps sC·r ≈ 0.11,
the lower r rides on a higher sC, i.e. the centre (sC̄, sV̄) = (A/r, τ/r) moves outward.
With ω = √(A·r) the rotation slows, so every wave lengthens — by a little in the realistic
band (first cycle 42.2 → 43.2 y; second cycle 45.0 → 46.6 → 48.4 y) and strongly at regime
scale — and the lengthening compounds from cycle to cycle because the r-deficit (hence the
ω-deficit) grows with time, exactly the spacing pattern of Section 4.2.

### 5.5 Why the upswing lengthens and the downswing does not

The asymmetry is produced by the fact that r declines monotonically (d ln r < 0 on every
segment of every run — the FROP drift) and faster the higher the CS:

1. *A rising target tilts the wave toward its rise.* The centre sC̄ = A/r grows at the rate
   h := d ln sC̄/dt = −d ln r/dt > 0 (about 0.3–0.9%/yr, larger under CS > 0). Since
   sC = sC̄·e^p, the mechanisation share rises whenever the orbit's own motion satisfies
   ṗ > −h, i.e. during more than half of every rotation. For a roughly sinusoidal p of
   amplitude ρ, the upswing share is

   ```
   upswing share ≈ ½ + (1/π)·arcsin(h/(ρ·ω))
   ```

   which is increasing in h and decreasing in ω — both pushed the right way by CS > 0 (a
   steeper r-decline raises h; a lower r lowers ω). Magnitude at cycle 1, CS = 0:
   h ≈ 0.0043/yr, ρ ≈ 0.16, ω ≈ 0.126/yr give arcsin(0.0043/0.020) ≈ 0.21 rad, i.e. a
   predicted share of ~57% against the measured 54% (the difference reflects the relaxation
   shape of the real orbit at large amplitude). The *increase* of the share with CS is the
   robust part: measured per-cycle upswing shares rise from 52.6% (CS = −0.10) to 54.2%
   (CS = 0) to 56.3% (CS = +0.10) on cycle 1 (upswings 23.2 → 24.4 → 25.9 y), and from
   51.1% to 52.4% to 54.5% on cycle 2 (upswings 24.2 → 26.0 → 28.3 y) — consistent with the
   median-based grid numbers 54 → 57 → 59% of Table 1. As the amplitude grows in later
   cycles, ρ·ω grows and the arcsin argument shrinks, so the share drifts back toward 50%
   even at fixed CS — also observed (cycle-2 shares below cycle-1 shares).
2. *The down-leg is trimmed because the overshoot above τ shrinks.* The downswing is the
   time h = sV·r stays above the fixed threshold τ. The measured overshoot above τ
   *decreases* with CS (mean sV·r over the first downswing: 0.064 → 0.061 → 0.059 for
   CS = −0.10/0/+0.10), because the employment term (sC·r − A)/sV acts on a larger sV when
   r is lower and compresses the sV excursion. A smaller excursion above a fixed threshold
   ends sooner: first downswing 20.9 → 20.6 → 20.1 y. (Later downswings still lengthen at
   fixed CS, because ω keeps falling with r, but always less than the upswings.)

Net effect: the lengthening of the cycle is carried almost entirely by the upswing
(≈ +1.5–2.3 y per +0.1 CS per cycle) against a flat-to-slightly-shorter downswing
(≈ −0.5–0 y), which is exactly the rise-share movement of Table 1.

### 5.6 Economic reading of the mechanism

The long wave of this model is a *growth-rate* oscillation: its clock is
ω = √((δ + τ + n)·r), the geometric mean of the natural-growth constant and the profit
rate. A class struggle that pushes wages above the reserve-army norm acts on the clock
through the exact term −a₂·CS in the profit-rate growth: it depresses r, slows the (sC, sV)
rotation around the receding balance point (τ/r, A/r), stretches each wave and — because
the receding centre rises faster under militancy — redistributes each wave's time toward its
mechanisation upswing and away from its high-hiring downswing.

## 6. Historical correspondence of the CS levels

Model time is not calendar time (the model is not calibrated to the data; cycle lengths
come out at 40–50 years). The correspondence below therefore compares *directions and
relative magnitudes* with documented historical episodes, using the wage share and the
profit rate as the bridge variables (both are what CS shifts in the model).

A convention fixes the units of the correspondence, because the model does not carry the
level of the wage share as a state variable (Section 2.3). The observable counterpart of CS
is the growth of the wage share relative to productivity, d ln w/dt ≈ ŵ − τ (real-wage
growth minus productivity growth). The "normal" wage dynamic — what the reserve-army
mechanism alone delivers — does not transmit productivity gains fully into wages: the
pressure of an increasing IRAL "keeps the wage share constant or even decreasing"
(book p. 216), i.e. the rate of surplus value tends to rise with productivity (relative
surplus value) unless labour intervenes. Hence the calibration adopted here: wages that
merely *keep pace* with productivity (ŵ ≈ τ, a flat share) are already an above-normal
outcome of struggle and register as **CS > 0**; wages that outrun productivity (a rising
share, the wage explosion of 1966–75) are a larger CS; and wages that fall behind even the
market norm (a falling share: the post-1980 offensive, the 1920s) are **CS < 0**. CS = 0 is
the benchmark in which the share is constant-or-decreasing under reserve-army pressure —
wages do not systematically keep up with productivity.

### 6.1 CS > 0 — wages at or above the productivity norm: the golden age and the wage explosion (1940–1975)

Model behaviour: profit-rate level squeezed (the CS component of d ln r/dt is −a₂·CS, i.e.
roughly −0.05% to −0.2% of the level per decade at CS = +0.005…+0.02, or −0.25% to −1% over
a typical half-century), waves slightly longer and more upswing-heavy, employment-share
troughs deeper during the squeeze, recovery delayed.

Historical counterpart, in two intensities of the same regime:

- *The golden age (1945–1966) — the postwar, institutionalised phase of the fourth-wave
  upswing: a moderate positive CS.* In the US-centred periodisation of the book's Table 5.1
  the fourth wave runs 1940–1982, so its upswing begins in 1940 (the exit from the interwar
  depression) and crests in the mid-1960s; the golden age proper is the postwar part of
  that upswing, 1945–1966. Organised labour stood at the height of its postwar power (US
  union density at its historical peak in the early 1950s, roughly a third of the private
  workforce), with full employment and institutionalised bargaining that tied wages to
  productivity growth (the Fordist wage formula). Because productivity was rising
  throughout, a wage share that merely *held its ground* at a high level (ŵ ≈ τ) already
  exceeded the reserve-army norm, under which the share is constant-or-decreasing (book
  p. 216): in the units adopted here this registers as a moderate CS > 0. The profit rate
  was nevertheless at its postwar maximum — which is exactly what the model predicts for
  such a moderate push: at CS = +0.005 the drag on the profit rate is only ≈ 0.05% per
  decade (Section 4.4; Table 1: mean r over [0,120 y] = 0.1407 at CS = +0.005 against
  0.1411 at CS = 0), and the strength of the upswing overwhelms it. The golden age is
  therefore not a counter-example to the wage-push mechanism: it is the case of a
  productivity-bound push, and high profits coexist with strong labour precisely because
  the push never outruns productivity.
- *1966–1975 — the wage explosion, at the crest and into the early downswing of the fourth
  wave: a larger CS.* Rank-and-file wage drift above the Fordist formula, strike waves
  (French May 1968; Italian *autunno caldo* 1969; UK miners 1972 and 1974; German wildcat
  strikes 1969–73), real-wage growth running above productivity growth, the US wage share
  rising from the mid-1960s to the mid-1970s: wages now outrun the productivity norm. The
  documented profit squeeze of the late boom followed (Glyn & Sutcliffe 1972 for the UK;
  Boddy & Crotty 1975 for the US; Armstrong, Glyn & Harrison 1991), and the measured pre-tax
  profit rate fell by roughly a third between 1965 and 1982 (Duménil & Lévy 1993) — the
  downswing of the fourth Kondratiev wave (1966–1982 in the periodisation of the book's
  Table 5.1). The model reproduces the *direction* (a sustained positive CS keeps the profit
  rate below the CS = 0 path and prolongs the phase of weak accumulation), though not the
  *magnitude* of the 1970s fall: in the calibration of the book that fall is driven mainly
  by the OCC channel (a₁ = 0.04) — the FROP of the downswing — with the wage push as a
  modulator (Sections 7–8).

### 6.2 CS ≈ 0 — the "normal" struggle level: the pure reserve-army cycle

Model behaviour: the book's own wave — upswing ~25 y, downswing ~19 y (rise share ~57%),
cycle length 43–47 y, profit rate falling secularly through the OCC channel.

Historical counterpart: periods in which wage formation was dominated by the labour market,
with organised labour too weak (or too integrated) to shift the wage baseline — above all
the 19th-century industrial cycles (weak unions, mass immigration, e.g. the US waves of
1845–1896). In such periods wages did not systematically keep pace with productivity
growth: under the pressure of the IRAL the wage share was constant or decreasing, the
mechanism the book itself describes on p. 216, and the wave ran with the book's own
rise/fall profile and no autonomous distributive component.

### 6.3 CS ≈ −0.005…−0.02 — a "below-normal" struggle: the post-1980 capital offensive

Model behaviour: profits kept high (mean r over [0,120 y] up to ~+0.009 above the CS = 0
case at CS = −0.10), wage-bill growth below the reserve-army norm, waves a little shorter,
upswings a little shorter relative to downswings (rise share down to 54–55%): accumulation
runs "cheaper" for capital.

Historical counterpart — 1979–2000: the monetarist shock and mass unemployment of
1979–82; the symbolic defeats of organised labour (PATCO, 1981; the UK miners, 1984–85;
German and Japanese wage moderation from the mid-1980s); US private-sector union density
falling from ~24% (1973) to ~10% (2010s), UK density from ~50% (1980) to ~23%; real wages
flat for the US production worker from 1973 to the mid-1990s while productivity kept
rising; the wage share of US national income falling from the mid-1970s/early 1980s back to
and below its 1965 level. Profits recovered correspondingly: the US profit rate rose from
its 1982 trough through the 1990s (Duménil & Lévy), and the period 1982–2007 is the upswing
of the fifth Kondratiev wave in the book's periodisation. Marx's *Capital* III, ch. 14
("depression of wages below the value of labour-power") names exactly this mechanism as a
counteracting cause of the FROP; the model's CS < 0 is its formalisation, and the model
confirms the expected effect: the fall of the profit rate is slowed.

### 6.4 CS ≲ −0.1 — extreme repression: the 1920s analogy

Model behaviour: profitability sustained near its initial level for more than a century of
model time (mean r 0.15–0.17 at CS = −0.2…−0.3, almost no squeeze), cycles shortened to
~41–44 y and nearly symmetric (rise share ~51%).

Historical counterpart — the US 1920s under the "American Plan" open-shop offensive: real
manufacturing wages roughly stagnant between 1923 and 1929 while output per worker rose by
about a quarter, pushing the wage share down and the profit share to a peak by 1929 — the
classic "Roaring Twenties" configuration of high and rising profitability with wage
repression (which, as the decade's end shows, stored up the overaccumulation crisis — a
crisis this five-equation model, having no credit/finance layer, cannot itself generate).
The earlier US episode of 1898–1907 (post-1896 wage-share trough, strike defeats of the
mid-1890s, then rapid accumulation) is a weaker, partial analogue.

### 6.5 CS ≳ +0.05 — regime-scale militancy: the squeeze that does not end

Model behaviour: the profit rate collapses towards 1–8%, the cycles stretch to 50–74 y, and
the economy slides into the low-profit crawl (the model leaves its economically meaningful
domain early, so these runs are boundary illustrations rather than forecasts).

Historical counterpart: no capitalist economy has sustained a 5–15%/yr excess wage push for
long — the resistance of capital, disinvestment and crisis have always cut it short — but
several *episodes* point that way: France 1936–1938 (Popular Front wage rises and the
40-hour week; the profit collapse of 1937 and the "capital strike" of 1937–38, ended by the
wage pause of November 1938); the UK 1974–1979 (wage push under full employment, the IMF
crisis of 1976, and the squeeze that prepared the 1979–80 reversal); the US 1970–1974 (wage
controls followed by the 1973–75 recession). These episodes are historical evidence that
regime-scale CS > 0 is not sustainable *as a constant*: the system responds with crisis and,
ultimately, a capital offensive that resets CS below zero — exactly the sequence stylised in
the pulse experiment of Section 4.5 and Figure 3.

## 7. Discussion: why a stronger working class lengthens the wave instead of bringing the crisis forward

The numerical results of Section 4 may seem to contradict a widely held Marxian intuition:
if a stronger working class raises the wage share and thereby lowers the rate of surplus
value, the rate of profit is squeezed; and since a low (or falling) rate of profit is the
source of crisis, one would expect militancy to *shorten* the boom and bring the crisis
forward. The model instead shows the opposite: sustained CS > 0 lengthens the waves, raises
the share of the wave spent in the upswing, and *delays* the reserve-army collapse, while
CS < 0 (wage repression) shortens the waves and produces the earlier, deeper collapses.
This section explains why this is not a contradiction — either with the model's own logic or
with Marxian crisis theory — and what the model actually says about the role of the class
struggle in the cycle.

### 7.1 The paradox, stated precisely

What the intuition asserts, and what the model confirms and refutes:

- *Confirmed.* CS > 0 lowers the rate of profit at every date and steepens its decline: the
  mean of r over [0, 120 y] falls from 0.1498 (CS = −0.10) to 0.1411 (CS = 0) to 0.1331
  (CS = +0.10), and d ln r/dt contains the exact drag −a₂·CS. In this narrow sense the
  "squeeze" is real: a powerful working class permanently depresses profitability.
- *Refuted.* The squeeze does not accelerate the crisis. The first sV < 0 crossing (the
  moment the employment share turns into absolute disinvestment — the model's reserve-army
  collapse) is *delayed* from 126 y (CS = −0.10) to 131 y (CS = 0) and 202 y (CS = +0.10),
  and at CS ≥ +0.15 it never occurs within the simulation window; the amplitude growth per
  cycle slows down; the wave lengthens (second cycle 45.0 → 48.4 y) and becomes more
  upswing-heavy (rise share 54 → 59%).

The apparent paradox is therefore: the wage push lowers the *level* of profitability, yet it
*mollifies* the cyclical movement instead of unleashing the crisis.

### 7.2 The resolution: the intuition confuses the turning point with the crisis, and the level with the dynamics

The resolution has three parts, all visible in the equations.

**(a) What actually ends the boom in this model.** The downturn is not triggered by the
level of the rate of surplus value; it is triggered by the state of the labour market. From
Section 5.1, the peak of the mechanisation wave occurs at the crossing sV·r = τ: the boom
ends precisely when the profit-rate-weighted hiring pressure overtakes the (quasi-fixed)
pace of technical change — i.e. when the reserve-army mechanism, through rising employment
and rising wages, chokes further accumulation. This is the profit-squeeze *turning point*,
and it is here that the naive intuition is correct: the wage push generated by full
employment is what terminates the boom. In the book's own language this is the cyclical
expression of the profit-squeeze reading of the FROP (p. 212).

**(b) What sets the size of the boom — and hence the severity of the correction.** The
*depth and violence* of the downswing are not set by the wage share but by the amplitude
that the accumulation process has built up during the preceding boom: the larger the
overshoots of sC and sV around their balance levels (A/r, τ/r), and the faster the
amplitude grows from cycle to cycle, the deeper the subsequent correction. The measured
amplitude dynamics confirm the role of profitability as the fuel of this process: under
wage repression (CS < 0) the oscillations grow faster and the reserve-army collapse comes
earlier (first sV < 0 at 126 y and, at CS ≤ −0.3, at ~78 y), whereas under CS > 0 the
amplitude growth slows down and the hiring-pressure variable sV·r barely overshoots its
threshold τ at the top of the boom (measured mean overshoot 0.064 → 0.059 across the grid,
Section 5.5) — the labour-market excess that ends the boom and powers the correction is
smaller. A powerful working class produces exactly this outcome: by keeping the profit
share lower in every phase of the boom, it deprives the accumulation process of the
excess that generates the crash — the rotation slows (Section 5.3: ω = √(A·r) falls with
r) and the system never develops the internal overaccumulation that produces the violent
collapse. The squeeze therefore acts as a permanent, early brake on the boom rather than
as the puncturing of a boom whose speed is given from outside.

**(c) The rate of accumulation is the independent variable, wages the dependent.** This is
the Marxian content of (a)–(b). In Marx's own treatment of accumulation (*Capital* I,
ch. 25), the wage is not the exogenous cause of the cycle: accumulation determines the
demand for labour, the reserve army regulates the wage, and a wage rise that begins to
encroach on accumulation merely slackens the accumulation that produced it — the rise is
self-limiting and cannot, in this logic, be the *first cause* of a crisis. The systemic
crisis belongs to the sphere of overaccumulation: capital accumulates beyond the point at
which it can be valorised, the organic composition rises and the rate of profit falls —
the mechanism that the book itself locates in the OCC-deviation term of Eq. (8.15)
(weighted by a₁ = 0.04, four times the weight of the wage-push bracket) and in the "moment
of over-accumulation" of its qualitative analysis (Section 8.2.1.2, p. 208, and Section
5.3.1 of the book). Overaccumulation is *fueled* by high profitability — that is, by the
exploitation of labour — not by the wage share. Hence the paradoxical inversion of the
intuition: **it is capital's victory over labour — the restoration of a high rate of
surplus value — that re-arms the boom-bust mechanism; a powerful working class converts the
crisis from an acute overaccumulation collapse into a chronic low-profit regime.**

Table 3 documents the second part of this statement numerically. Over the first 120 years,
a 12% spread of the mean profit rate across the CS grid leaves the *rate of accumulation of
constant capital* (sC·r ≈ 0.11, i.e. g_C = sC·r − δ ≈ 0.06/yr) almost unchanged — the
system accumulates at the same pace on a lower rate of profit, in an increasingly mechanised
form (the mean share sC rises from 0.754 to 0.895 and the mean growth rate of sC turns
positive and increases with CS: capital answers the wage push by substituting dead for
living labour, which is precisely the OCC channel that accelerates the *secular* fall of the
profit rate). The rate of profit falls while the mass and the pace of accumulation are
maintained — Marx's own characterisation of the FROP (*Capital* III, ch. 13–15).

**Table 3 — Means over [0, 120 y]: the profit rate falls with CS while the pace of
accumulation (sC·r) is essentially unchanged.**

| CS | mean r | mean sC | mean sV | mean sC·r | mean growth of sC (per yr) |
|---|---:|---:|---:|---:|---:|
| −0.10 | 0.1498 | 0.754 | 0.326 | 0.1116 | −0.0003 |
| 0 | 0.1411 | 0.821 | 0.339 | 0.1129 | +0.0013 |
| +0.10 | 0.1331 | 0.895 | 0.349 | 0.1140 | +0.0032 |

### 7.3 The book's own phase analysis already contains this logic

The qualitative analysis of Section 8.2.1.1 arrives at exactly the same distinction. There,
the sign of g_e selects the regime of the two-variable system: with g_e < 0 (the permanent
squeeze) the equilibrium point B is a *stable node* that attracts all solutions — the
economy is drawn toward the configuration of zero profitability ("a pure labour economy…
no different from Smith's 'rude and early state of society'", p. 212), i.e. toward
stagnation, not toward an explosive crisis; with g_e > 0 the equilibrium becomes a
repelling saddle and the solutions travel along the accumulation path toward the point of
overaccumulation (pp. 207–212). In other words, in the theory of the book itself the
permanent wage squeeze is the *stabilising, stagnation-generating* case, while the 
restoration of exploitation is the case that generates the dynamic of accumulation and
crisis. The CS results of Section 4 are the five-equation counterpart of that analysis:
high CS behaves like the g_e < 0 regime (the oscillation loses amplitude growth and the
economy sinks into a low-profit crawl), while CS < 0 behaves like the g_e > 0 regime
(expanding, shorter, downswing-heavy cycles ending in the early reserve-army collapses of
Table 1).

### 7.4 The two configurations of crisis in the history of capitalism

Read together with Section 6, the model suggests that Marxian crisis theory should
distinguish two configurations, and that the class struggle selects between them rather
than between "crisis" and "no crisis":

1. *Chronic low-profit stagnation (CS > 0) — the model's counterfactual regime.* The
   squeeze lowers profitability; accumulation continues, but on a low rate of profit, in an
   increasingly mechanised form (Table 3); growth is sluggish, the waves are long, and there
   is no single catastrophic collapse — the system "ages" in the low-profit crawl. Two
   things must be said about the historical correspondence of this regime, and they require
   care with the periodisation.
   First, the model's CS is the *autonomous wage push* — wage growth above the norm that the
   labour market and productivity set — not the institutional strength of the working class
   as such. This distinction matters for the golden age, the postwar, institutionalised
   phase of the fourth-wave upswing (in the periodisation of the book's Table 5.1 the fourth
   wave runs 1940–1982, so its upswing begins in 1940 and the golden age proper covers
   1945–1966) — a phase in which a powerful labour movement coexisted with a profit rate at
   its postwar maximum. The point is not, however, that the era registers as
   CS ≈ 0. Because productivity was rising throughout, a wage share that merely *held its
   ground* at a high level (wages growing with productivity, ŵ ≈ τ) already exceeded the
   reserve-army norm, under which the share is constant-or-decreasing (book p. 216); in the
   units fixed in Section 6 the golden age is therefore a *moderate positive* CS — a
   labour movement strong enough to hold the share at a high plateau against the
   productivity trend. The golden age then ceases to look paradoxical at all: the model
   predicts that such a moderate push costs the profit rate only ≈ 0.05% per decade
   (Section 4.4; Table 1: mean r over [0,120 y] = 0.1407 at CS = +0.005 against 0.1411 at
   CS = 0), a drag that the high-accumulation upswing (with its rapid growth of the labour
   supply) simply overwhelms. What would squeeze profits is not a strong labour movement as
   such but a push that *outruns productivity* — which is what makes the second intensity of
   the same regime, the wage explosion of 1966–1975, the genuine squeeze episode.
   Second, that *above-productivity* push — a larger CS — appears with the wage explosion
   of 1966–1975 (Section 6.1) and, in the periodisation of the book's Table 5.1, sits at
   the crest and in the *downswing* of the fourth wave: the wave crests in 1966 (Table 5.1
   of the book dates the fourth-wave downswing 1966–1982), so by the late 1960s the long
   wave is already in its descending phase, and the low profitability of the mid-1970s to
   early 1980s is a phenomenon of that downswing, not the manifestation of a sustained
   large-CS regime. The
   wage explosion played the profit-squeeze role at the turning point (the boom ends in full
   employment and rising wages), but the depth of the 1974–82 fall of the profit rate
   belongs primarily to the downswing's own dynamics — the OCC/mechanisation channel of Eq.
   (8.15), the FROP — with the wage push as a modulator; Section 6.1 notes that the CS
   channel alone cannot reproduce the magnitude of the 1970s fall. The correct historical
   reading is therefore the pulse template of Section 4.5: a *bounded* CS > 0 episode
   (1966–75) coincided with the depressive phase of the fourth wave; the crisis form of that
   phase — stagflation, repeated recessions (1973–75, 1980–82), a profit rate that recovered
   only after the offensive of 1979–82 — is the historical form of the *profitability crisis
   of a long-wave downswing*, not a 1929-type overaccumulation crash. That the wage-share
   peaks of 1974–75 did not produce a 1929 is, in this framework, exactly what the two
   configurations predict: the deepest crashes are prepared by the *restoration* of
   exploitation (configuration 2 below), not by its squeeze.
2. *Acute overaccumulation crash (CS < 0).* Wage repression raises the profit share; the
   boom runs fast, the overshoot of accumulation becomes large, and the system generates its
   own violent correction (in the model: short cycles, downswing-heavy, early sV < 0; at
   CS ≤ −0.2 the economy leaves its meaningful domain after 78–88 years with the sharpest
   reserve-army collapses of the grid). The historical counterparts are the profit booms of
   the 1920s and of 1982–2007 — both prepared by successful offensives against labour —
   which were broken by the deepest crashes of the twentieth and twenty-first centuries
   (1929–33; 2007–08, with the credit system, absent from this model, as the transmission
   mechanism).

The genuinely counter-intuitive — but, in this framework, unavoidable — conclusion is the
second configuration: the most violent crises of accumulation are prepared not by the wage
push of a strong working class but by the profit booms that follow its defeat. This is not a
denial of the profit-squeeze school: the squeeze remains the mechanism of the *turning
point* (the boom ends in full employment and rising wages, exactly as in the crossing
sV·r = τ of Section 5.1). What the model denies is that the wage share is the *depth*
variable of the crisis: the depth is set by the overaccumulation that high profitability
itself creates.

### 7.5 The reformist illusion in the mathematics of the model

Sections 7.1–7.4 have a political-economic corollary. The *reformist illusion* is the
expectation, strongest in the social-democratic tradition, that the accumulation mechanism
itself — above all the absorption of the industrial reserve army in the boom — will
progressively improve the position of labour within capitalism: full employment, rising
wages, shared prosperity, delivered by the market, so that workers need only be patient
(and moderate in their demands) for accumulation to do the work. The mathematical structure
of the Marxian accumulation cycle defeats that expectation twice, and CS sharpens the second
defeat.

1. *Even at CS = 0 the absorption of the reserve army is the mechanism that ends the boom,
   not one that can continue indefinitely.* Employment is absorbed while the hiring share
   sV > 0, and the boom ends exactly when the hiring-pressure variable reaches the choke
   point sV·r = τ — the sC peak (Section 5.1): the wage push generated by absorption itself
   cuts profitability, accumulation slackens, and mechanisation (the rising centre
   sC̄ = A/r, Table 3) recreates the reserve army. "Full absorption with rising wages" is
   therefore not a state of the system; it is the turning point of the wave. At the book's
   parameters the system has no positive fixed point at all (findings F1–F4), so there is no
   equilibrium in which the reformist promise could rest. This is the Goodwin core of the
   argument: the wage share and employment form a predator–prey pair whose "prosperity"
   configuration is cyclically self-reversing (Marx, *Capital* I, ch. 25).
2. *CS > 0 removes even the transient basis of the expectation.* Reformism counsels
   patience: let the boom run, and the absorption of the reserve army will tighten the
   labour market and raise wages "automatically". A positive CS means wages rise before,
   and independently of, the absorption of the reserve army — the push acts at any size of
   the reserve army. The profitability of the absorption process is then squeezed from the
   outset: the mean profit rate over the first upswing falls from ≈ 0.163 (CS = −0.10) to
   ≈ 0.150 (CS = +0.10), and the hiring-pressure product sV·r overshoots its choke
   threshold τ less (0.064 → 0.059, Section 5.5). The scenario in which capital *rapidly*
   absorbs the reserve army on high and rising profits — the empirical spectacle that
   nourishes the reformist expectation — is precisely the scenario that a larger CS
   suppresses, because the wage push binds before the reserve army is exhausted.

One qualification matters, and it is the counter-intuitive result of Sections 7.1–7.2
again. A larger CS does *not* abolish the absorption phase, and it does not weaken the pace
of accumulation itself: over the first upswing the accumulation rate of constant capital is
exactly invariant across the grid (mean sC·r = 0.1144 at CS = −0.10, 0 and +0.10), and the
hiring share sV stays positive longer (the first sV < 0 crossing is delayed from 131 y at
CS = 0 to 202 y at CS = +0.10 and does not occur within 300 y for CS ≥ +0.15). What CS
removes is the *profitability* of the absorption: the same accumulation proceeds on a lower
and falling rate of profit, in an increasingly mechanised form (Table 3). The militant
regime therefore denies reformism its material basis in a different way than mass
unemployment would — by draining the boom of the profit that sustains the promise of shared
prosperity through the market. Conversely, the historical era in which absorption and high
profitability did coexist — the golden age of 1945–1966, with its moderate, productivity-
bound CS (Section 6.1) — is precisely the era in which the reformist illusion flourished;
and it was the wage explosion of 1966–75, together with the profit squeeze of the
fourth-wave downswing, that broke the postwar social-democratic full-employment
expectation. The formal statement of that break is Przeworski's: wage gains that outrun
accumulation reduce investment, so the reformist strategy cannot deliver its promise within
the accumulation mechanism (Przeworski 1985).

In the terms of the model: the reformist illusion rests on the hope that the wage share can
rise *through* accumulation; the model shows that the wage share can rise only *against*
accumulation (CS > 0) — and that even then the gain is purchased at the price of
profitability, a price that capital makes labour repay in the subsequent offensive
(Sections 6.3 and 7.4).

### 7.6 Why CS > 0 shortens the individual bust but prolongs the low-profit state — and the threshold beyond which recovery fails

The results of Section 4.3 (the measured downswing of sC is flat or slightly *shorter* under
CS > 0) and the stagnation discussion of Sections 7.4–7.5 can seem to contradict each other.
They do not: "the depression" is used in two senses that respond to CS in opposite ways, and
separating them also answers the threshold question — recovery in the strong sense fails
only above a well-defined CS.

**(1) The bust segment vs the low-profit state.** The measured "downswing" of Table 1 is
the time between the peak and the trough of the mechanisation wave sC (≈ 19–24 y): the time
the hiring-pressure variable needs to fall from its overshoot back to the choke threshold τ
(Section 5.1). Under CS > 0 this segment is mildly *shortened* at moderate values (median
fall 19.9 y at CS = −0.10 → 19.4 y at CS = 0 → 19.1 y at CS = +0.05), because the wage push
acts as an early brake on the boom: a boom that runs on lower profitability develops smaller
excesses, so the hiring-pressure overshoot above τ is smaller (0.064 → 0.059, Section 5.5)
and the correction is over sooner. The time is not lost but shifted — the upswing lengthens
by more than the downswing shortens (rise share 54 → 59%). This is the sense in which the
numerical results show a "shorter depression" under CS > 0: the *bust segment* shortens.
The second sense is the *low-profit state* — the level of the profit rate around which the
cycles operate. Here CS > 0 shortens nothing: r lies below the CS = 0 path at every date
(r(120 y) = 0.0997 at CS = 0 against 0.0860 at CS = +0.10), and the deficit widens with time,
because the implied rate of surplus value is never restored while the push persists:
g_e = −sV − CS ≤ 0 as long as the hiring share stays above −CS. Every new upswing therefore
restarts from a lower profit rate — the system *ratchets down* even though each individual
bust is short: short recessions into a permanently lower profitability. That ratchet, not
the bust length, is what configuration 1 of Section 7.4 calls stagnation.

**(2) The threshold.** "Recovery" in the strong sense is the point at which the wage-share
dynamic reverses: g_e turns positive, so that the rate of surplus value — whose sustained
suppression is what drives the profit rate toward zero — stops being eroded. This requires
the reserve-army purge to overpower the push: g_e > 0 requires sV < −CS, i.e. the downswing
must shed employment deeply enough. The discipline the market can deliver is bounded, so the
question "does the purge ever occur?" has a threshold answer at the book's parameters. The
first sV < 0 crossing still exists at CS = +0.10 (t = 202 y) and at CS = +0.11 (t = 205 y,
and only to sV ≈ −0.03), and disappears between CS = +0.11 and +0.12: **for CS ≳ +0.115 the
wage share is never purged within 300 y** — g_e ≤ 0 over the whole horizon and r slides
monotonically towards zero (r(300 y) = 0.0076 at CS = +0.10, 0.0058 at +0.12, 0.0024 at
+0.20). Above a
second, much higher threshold the oscillation itself dies: the completed sC cycles per 300 y
fall from three at CS = +0.30 to one at CS = +0.50, the last extremum of sC moving from
t ≈ 185 y (CS = 0.30) to t ≈ 127 y (CS = 0.35–0.45) to t ≈ 67–69 y (CS = 0.50–0.60), after
which the trajectory is monotone — accumulation in an increasingly mechanised form on a
vanishing profit rate. Recovery in the sense that the boom regenerates itself fails
completely only above this second threshold (roughly CS ≈ 0.35–0.5 here: a permanent excess
wage push of 35–50% per year — an extreme counterfactual; by then the premise sC + sV ≤ 1
is long violated, so these are extrapolated regime statements, not economically meaningful
trajectories).

Between CS = 0 and the first threshold the behaviour is smooth, not bimodal. Two effects
compete: the early-brake effect (smaller overshoot, trims the bust) and the slow-clock
effect (ω = √(A·r) falls with r, lengthening every segment, Section 5.3). The second
gradually wins as the r-deficit accumulates: the measured falls of sC run 20.9 → 20.6 →
20.1 y in cycle 1 (CS = −0.10/0/+0.10), are ≈ 23–24 y in cycle 2 at every CS, and stretch to
27.4 y in cycle 3 at CS = +0.10 (a cycle the CS = 0 economy no longer completes inside its
meaningful domain). There is therefore no single critical CS for the length of individual
busts; there is a critical CS for the *purge* (≈ +0.11–0.12 at these parameters), and a much
higher one beyond which even the cycles cease.

**(3) What this means for the 1970s question.** The 1970s-scale push (CS of the order of
0.01–0.03) lies far below the purge threshold, and the model does not claim that such a push
alone prevents recovery — its direct channel (a₂ = 0.01) is small and the level of the wage
share is not a state variable (Sections 2.3 and 7.8). The model's contribution is structural
and qualitative: (i) the variable that a persistent wage push prolongs is not the individual
bust but the low-profit state — the ratchet, whose exit requires the wage push to be broken;
(ii) if the push were maintained above the purge threshold, no restoration would be possible
at all and the low-profit phase would have no internal end (the extrapolated regime of
configuration 1). The profit depression of 1974–82 is the historical ratchet phase; its exit
required CS to fall below zero — the subject of Section 7.7.

**(4) Fixed CS versus regime change: how the runs are to be read.** CS is a parameter: each
run of the ODE keeps CS constant, and this may seem to make "the exit required CS to fall
below zero" meaningless. The resolution has three parts.

- *Within a fixed-CS run, cycle restarts do occur even at CS > 0.* The upswing of sC does
  not wait for the restoration of the rate of surplus value: the boom restarts when the
  hiring-pressure variable falls back below the choke threshold (the sC trough is the
  crossing sV·r = τ from above, Section 5.1), and the profit rate simply keeps declining
  through successive cycles. A run at fixed CS = +0.02 therefore produces repeated
  "recoveries" — each of them profitless in the sense that it starts from a lower r than the
  previous one (the ratchet of point (1)). This is exactly the historical texture of
  1974–82: output and employment recovered in 1975–79 while the profit rate never returned
  to its pre-1973 level, and each recovery stalled again until the regime changed.
- *What a fixed-CS run cannot do is restore profitability.* At any CS ≥ 0 the implied rate
  of surplus value is not restored (point (2)), so no single trajectory with the 1970s
  parameter value can produce the return of high profits. The historical "exit" is therefore
  not an event inside such a trajectory: it is a *change of the parameter itself* — the
  transition to the CS < 0 regime of Section 7.7, represented by a new run (or, within one
  integration, by a piecewise-constant CS — exactly what the pulse experiment of Section
  4.5 already does, where CS switches from +0.02 to +0.005 to −0.015 along a single
  trajectory).
- *The regime change is exogenous to the model.* Nothing in Eqs. (8.15-CS)–(8.19) turns a
  militant regime into an offensive one: the equations do not explain why the wage push is
  broken in 1979–85 (political defeats, mass unemployment as a policy, the state). The model
  treats the class-struggle level as an exogenous parameter, exactly as it treats the
  population growth n and the "normal" growth g_w; the historical narrative of Sections 6–7
  is a *succession of runs* with different parameters, one per epoch — the same methodology
  with which the companion repository compares the waves of Table 5.1 (one run per wave, with
  that wave's n). What the model *does* contribute to the transition is the demonstration
  that the militant regime is unsustainable in the long run — at fixed CS > 0 the profit rate
  slides toward zero and the shares leave their economically meaningful domain — so that the
  political flip is not arbitrary but is the real-world answer to a degeneracy that the
  equations exhibit.

**(5) How long the pointless recoveries last, and how they end.** Figure 4 shows a run with
CS fixed at +0.02 — the 1970s-scale push held constant — against the CS = 0 baseline. The
ratchet is visible cycle by cycle. The accumulation wave continues almost on schedule: the
boom peaks of sC occur at t ≈ 50, 97, 150 and 213 y (peak-to-peak spacings 46.9, 53.0 and
62.8 y, close to the CS = 0 cadence of 46.6, 52.3 and 61.2 y), because the restart of the
boom needs no restoration of the rate of surplus value — only the return of the hiring
pressure below the choke threshold (the green bands of Figure 4's third panel, the phases
sV·r < τ, are exactly the recoveries). But every boom peaks at a lower profit rate than the
last: r at the successive peaks is 0.146 → 0.115 → 0.083 → 0.051 (against 0.148 → 0.117 →
0.086 → 0.053 at CS = 0 at the same dates), and every recovery starts from a lower trough
(r = 0.163 → 0.130 → 0.096 → 0.062 at t ≈ 26, 71, 121, 178 y). Over the economically
meaningful horizon the profit gap to the CS = 0 path widens with every cycle (peak gaps of
≈ 0.002–0.003 at each of the first three peaks; r(120 y) = 0.0967 against 0.0997) — the
accumulated −a₂·CS·t drag plus the state feedback. (In the later, out-of-domain phase the
CS = 0 baseline's own deep crisis — sV < 0 from t ≈ 131 y against ≈ 133 y at CS = +0.02 —
accelerates its decline and narrows the gap again, which is why the comparison is
meaningful only while the shares respect the book's premise.)

Each pointless recovery — one sC upswing — lasts ≈ 25–28 y of model time, and the
succession has no internal end: as long as the wage push holds, the mechanism that restarts
the boom (sV·r falling back below τ) keeps operating, so the recoveries keep coming
indefinitely at fixed CS. The succession ends only in one of two ways. (i) In the model's
own terms the run eventually leaves its economically meaningful domain: the mechanisation
share exceeds the whole surplus value at t ≈ 86 y (sC > 1, the book's premise sC + sV ≤ 1),
and the reserve-army purge — sV < 0, the event that would finally reverse the wage-share
dynamic — arrives only at t ≈ 133 y, because CS = +0.02 lies below the purge threshold of
≈ +0.115 (point (2) above). By then the shares have long violated the book's premise, so
this late purge is an out-of-domain artefact rather than a meaningful restoration: r(220 y)
= 0.045 and the slide continues. (ii) Historically, the succession was ended not by the
model's internal dynamics but by the regime change — CS switched below zero (Section 7.7);
had the push instead been maintained above the purge threshold, even the late purge would
not occur and r would slide toward ≈ 0.006 by 300 y (point (2) above). The conclusion of
the exercise is the one drawn in points (1)–(4): a militant wage push does not stop the
cycle — it empties the recoveries of profit; and because profit is the regulator of
accumulation, a succession of profitless recoveries is not a resolution but a ratchet,
whose only exits are the domain breakdown of the equations or a change in the class
relation. (The correspondence to the 1970s is of the *pattern*, not of the frequency: the
cyclical recoveries of 1975–79 and 1980–81 inside the long downswing restarted activity but
not profitability; the model's wave is the long wave, so its model-time durations should
not be read as calendar durations.)

### 7.7 A stalemated class struggle and the exit from the low-profit state: the offensive and the external labour supply

The question addressed here: if, during the downswing, capital cannot dismantle the
political and organisational power of labour, so that the class struggle remains stalemated
at CS ≥ 0 (the 1970s situation), how could the centre economies have exited the low-profit
state — and how, had no large external labour pool (East Asia, later China) been available
for industrial relocation, could they have done so at all? Section 7.6 fixes the model's
answer to the first part: a stalemate does not lengthen the individual busts; it keeps the
ratchet engaged — the rate of surplus value is not restored while CS ≥ 0, so profitability
stays depressed for as long as the push persists, and above the purge threshold
(CS ≳ +0.11–0.12 at the book's parameters) it cannot be restored at all. Consistent with
Section 7.6 (4), "the exit" is a regime change — the parameter is switched to CS < 0, and
the post-1979 history is represented by the run of that new regime — not an event inside
the stalemated trajectory. The 1974–82
episode fits this reading qualitatively: the wage share held at its 1974–75 peak through the
late 1970s, the organised power of labour survived the two recessions intact, and the
measured profit rate stayed depressed until the offensive of 1979–82 had begun
(Duménil & Lévy 1993). (The caution of Section 7.6 applies: the model's wage-push channel
carries the small weight a₂ = 0.01 and the level of the wage share is not a state variable,
so the model reproduces the qualitative condition — no purge, no restoration — not the full
quantitative depth of the 1970s profit fall, which the downswing's own FROP dynamics carry.)

The exit required CS to fall below zero — the wage push to be broken — and it happened in
two steps:

1. *The internal offensive of 1979–85* (monetarist disinflation and the deliberate
   unemployment of 1980–82; the defeats of organised labour: PATCO 1981, the UK miners
   1984–85; wage concessions and concession bargaining): in the model this is the parameter
   shift to CS < 0, which re-opens the reserve-army channel — the purge arrives earlier
   (sV < 0 at ≈ 126 y at CS = −0.10, ≈ 78 y at CS = −0.30), the rate of surplus value is
   restored, and profitability recovers (r stays at 0.15–0.17 at CS ≤ −0.10, against 0.10
   at CS = 0), at the price of sharper, more downswing-heavy cycles (configuration 2 of
   Section 7.4). An additional, unplanned purge operated through the inflation of 1974–80,
   which eroded real wages (a nominal channel outside this real model).
2. *The external labour supply.* Relocation and the threat of relocation to lower-wage
   regions — East Asia and Latin America from the late 1960s–70s, China from the 1990s and
   especially after its WTO entry in 2001 — enlarged the reserve army facing the centre's
   workers without a frontal political defeat: the reserve army of Marx's general law of
   accumulation became global. In the units of the model this works like a sustained
   negative CS on the centre: US real wages stopped growing with productivity from the
   1970s and stayed flat even at full employment after 2001 (the "China shock":
   Autor, Dorn & Hanson 2013; Freeman 2006 — the world labour force roughly doubled with
   the integration of China, India and the former Soviet bloc). The upswing of the fifth
   wave (1982–2007 in the book's periodisation) is, in this reading, the wave of global
   labour integration; the companion repository's price evidence dates the synchronisation
   of the emerging economies with the US cycle to ~1979.

**The counterfactual: no external labour pool.** Without the possibility of relocation, the
centre would have had to exit the stalemate through the internal routes alone: either a
longer and deeper slump that devalued capital and purged the wage share by mass
unemployment — the template of 1873–96 and 1929–40, long depressions that ended only when
the internal balance of forces was settled (in 1929–40, only the war and the destruction of
the interwar order accomplished it) — or, if labour held, the indefinite low-profit ratchet
of configuration 1. And even after a successful internal purge, the restored-profit boom
would have run into domestic labour scarcity much sooner: without the global pool, the wage
squeeze would have returned at every approach to full employment, and the centre would have
been condemned to either recurring frontal confrontations over distribution or to keeping
its own reserve army large through repeated deliberate recessions. The global labour supply
was therefore not merely one wage-restraint device among others: it is what allowed high
profitability to coexist with (near-)full employment in the centre for a generation — the
combination that a purely internal class compromise could not deliver.

### 7.8 Limits of this interpretation

1. *Permanent vs episodic struggle.* The scenarios of Section 4 hold CS constant, whereas
   the class struggle in history is episodic and answered by capital (disinvestment,
   restructuring, the state). The pulse experiment of Section 4.5 is the more realistic
   template: a bounded militant episode produces a bounded squeeze, and the subsequent
   offensive restores profitability — after which (historically) the overaccumulation
   mechanism of Section 7.4 resumes. The permanent-CS runs should be read as
   counterfactuals: "what if the wage push had been sustained?" — and the low-profit crawl
   they produce is precisely the outcome that capital's counter-offensive exists to prevent.
2. *The wage share is not bounded in the model.* The rate of surplus value e is not a state
   variable (Section 2.3), so a permanently negative g_e can, in principle, drive the
   implied wage share toward unity. In reality e has a floor set by the reproduction
   requirements of labour power and by the resistance of capital; the regime-scale runs of
   Section 4.2.4 should therefore be treated as boundary illustrations, and the historical
   record confirms that they are not sustainable as constants (Section 6.5).
3. *No credit and finance.* The historical crashes of configuration 2 are transmitted
   through the credit system, which is absent from Eqs. (8.15)–(8.19); the model captures
   the real-side mechanism (profitability → overaccumulation → collapse of accumulation),
   not the financial amplification.
4. *The upswing that lengthens under CS > 0 is the mechanisation upswing of sC* — a boom
   that continues at low profitability precisely because accumulation proceeds in an
   increasingly mechanised form (Table 3). This is a Marxian configuration (accumulation
   against a falling rate of profit, the substitution of dead for living labour), but it is
   not a boom in employment or in the wage bill; the sV share is correspondingly flatter
   (rise share ~50–52% at all CS).

## 8. Conclusions and limits

1. **The observation is correct, with qualifications.** Sections 8.2.1.1–8.2.1.2 (book pp.
   205–212) acknowledge that the sign and intensity of g_e carry the class struggle, and the
   3rd case (g_e < 0, profit-squeeze, p. 212) makes the wage push the reason behind the
   FROP — but the closing relation g_e ≅ −sV (p. 214) reduces the wage dynamic to the
   reserve army, and the complete five-equation model of Section 8.3 contains no variable
   through which a struggle "aiming higher than what the investment in variable capital
   dictates" could act. The book's later statements that the rate of surplus value reflects
   "class struggle and the like" (p. 238) and "the level of class struggle in the sphere of
   circulation" (p. 239) characterise the concept of the model, not its equations.
2. **The CS extension** (g_e = −sV − CS; only Eq. 8.15 changes) restores the struggle in a
   minimal, well-defined way: CS > 0 is an above-normal wage push (the "historical and
   moral element" of the value of labour-power), CS < 0 a below-normal push (Marx's
   counteracting cause, the depression of wages below value). CS = 0 reproduces the book
   exactly.
3. **What CS does to the wave** (the quantitative answer):
   - the *local (linearised) period* is unchanged (41.7–41.8 y) — CS shifts only the weak
     (1,1) Jacobian entry; the local clock is set by the OCC mechanism;
   - the *nonlinear cycle* lengthens mildly with militancy (≈ +0.05 y per +0.01 CS on the
     first and ≈ +0.16–0.18 y per +0.01 CS on the second completed cycle; stretching
     strongly only at regime-scale CS ≥ +0.05), and shortens mildly under repression;
   - the *upswing* lengthens and the *downswing* stays flat or shortens slightly (rise
     share of the sC wave: 54% at CS = −0.10 → 57% at CS = 0 → 59% at CS = +0.10);
   - the *level effect dominates*: the CS component of d ln r/dt is −a₂·CS = −0.01·CS
     exactly, so a sustained above-normal struggle of 1%/yr (CS = +0.01) costs ≈ 0.1% of
     the profit rate per decade (≈ 1% over a century), and regime-scale pushes
     (CS ≥ +0.05…+0.10) drive the profit rate towards zero over the simulation horizon;
     wage repression keeps profitability high and slows the FROP;
   - the *mechanism* behind the cycle-shape effects is the rotation of the share dynamics
     around the receding balance point (τ/r, A/r) with the clock ω = √(A·r): CS depresses
     r and thereby slows and tilts the wave (Section 5).
4. **Historical correspondence** (directions, not point calibration; units as fixed in
   Section 6: CS = 0 means the wage share does not keep pace with productivity under
   reserve-army pressure): CS > 0 ↔ wages at or above the productivity norm — the golden
   age on the fourth-wave upswing of 1940–1966 (its institutionalised postwar phase, the
   golden age proper of 1945–1966: a moderate positive CS with a high profit rate),
   intensifying into the 1966–1975 wage explosion and the profit squeeze of the fourth-wave
   downswing (1966–1982); CS ≈ 0 ↔ the reserve-army-dominated 19th-century waves; CS < 0 ↔
   the post-1980 offensive (decline of union density, falling wage share, profit recovery,
   the fifth-wave upswing of 1982–2007); CS ≪ 0 ↔ the 1920s "open-shop" decade; regime-scale
   CS > 0 is unsustainable and historically appears only as bounded episodes (France
   1936–38; the UK 1974–79), each terminated by crisis and a subsequent offensive.
5. **Limits of the exercise**:
   - the direct quantitative channel of the wage push is small *by the calibration of the
     book*: g_e enters only the "counteracting forces" bracket of Eq. (8.15), weighted with
     a₂ = 0.01 against a₁ = 0.04 for the OCC channel (p. 215). Inserting CS therefore
     changes the profit-rate *level* first and the *clock* second — a finding in its own
     right, since it shows that even a formal admission of the class struggle would leave
     the model largely OCC-determined;
   - the level of distribution is not a state variable of the model (the Goodwin
     wage-share equation is dropped, p. 216), so a *level* shift of the wage share —
     historically the most visible product of class struggle — cannot be represented, only
     the growth-rate premium CS. A fuller extension would add the rate of surplus value e
     (or the wage share) as a sixth state with a struggle-augmented adjustment equation;
     that is the natural next step and would make the historical magnitudes (wage-share
     swings of ±4–7 points) commensurable with the model;
   - the model has no positive fixed point at the printed parameters and leaves its
     economically meaningful domain after ~85–130 y (findings F1–F6 of the companion
     repository), so the quantitative statements refer to the first two completed cycles;
   - model time is not calibrated to calendar time; the historical correspondence is
     qualitative.

## References

- Armstrong, P., Glyn, A. & Harrison, J. (1991). *Capitalism since 1945*. Blackwell.
- Autor, D. H., Dorn, D. & Hanson, G. H. (2013). The China syndrome: local labor market
  effects of import competition in the United States. *American Economic Review*, 103(6),
  2121–2168.
- Freeman, R. B. (2006). The great doubling: the challenge of the new global labor market.
  In *Emerging Economies and the Transformation of International Business*. Edward Elgar.
- Boddy, R. & Crotty, J. (1975). Class conflict and macro-policy: the political business
  cycle. *Review of Radical Political Economics*, 7(1), 1–19.
- Chatzarakis, N., Tsaliki, P. & Tsoulfidis, L. (2022). *Economic Growth and Long Cycles: A
  Classical Political Economy Approach*. Routledge. Sections 8.2.1.1–8.2.1.3
  (pp. 205–215), 8.2.2 (p. 216), 8.3 (pp. 221–235), 8.4 (pp. 236–242).
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
- Przeworski, A. (1985). *Capitalism and Social Democracy*. Cambridge University Press.
- Companion repository reports: `docs/model_validation_n_scenarios.md`,
  `docs/sv_sign_flip_report.md`, `docs/Technical_Report_Ch8_Model_Errata.docx`
  (findings F1–F6 on the printed Chapter 8 material).
