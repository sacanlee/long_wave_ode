# Does Flipping the Sign of sV in Eq. (8.15) Fix the Reported Inconsistencies?

**Author**: sacanlee | **Date**: 2026-08-28 | **Script**: `scripts/sv_sign_flip.py` | **Results**: `data/sv_flip_results.json`
**Subject**: Chatzarakis, Tsaliki & Tsoulfidis (2022), *Economic Growth and Long Cycles*, Routledge, Section 8.3.3
**Question**: If the sign of `sV` in Eq. (8.15) is changed from "−" to "+", are the six internal inconsistencies reported earlier (Technical_Report_Ch8_Model_Errata.docx) resolved?

---

## 1. The printed Eq. (8.15) and the origin of the sV signs

Glyph-level verification against the PDF (book pp. 213–214) shows the printed (8.15) is:

**ṙ = −a₁·(sC − sV)·r² + a₂·(δ + τ − sV)·r** (8.15)

`sV` appears twice, both times with a **minus sign**. Its economic meaning follows from the book's own derivation chain:

1. **(8.14)**: ṙ = −(sC − sV)·r + δ + τ + g_e, where g_e = the growth rate of the rate of surplus value;
2. Book p. 213: "According to Marx (Capital I, ch. 25), the growth rate of the rate of surplus value, g_e, is **inversely proportional** to the change in the labor force, or, what is the same, to the investment in variable capital, sV";
3. **Glombowski (1983)**, "A Marxian model of long run capitalist development", *Zeitschrift für Nationalökonomie* 43(4), 363–382, hypothesizes **g_e = ε₀ − ε₁·sV** (ε₀, ε₁ small positive constants);
4. The book adopts the simpler relation **g_e ≅ −sV**, substituted into (8.14) to obtain (8.15);
5. Parameter values (a₁, a₂, n, δ, τ) come from Sasaki (2013) (note 16 of the chapter).

**Economic meaning**: sV = the share of surplus value invested in variable capital (employment / the wage fund).
- **Second bracket (δ + τ − sV)**: δ (devaluation) and τ (technological change) are the *counteracting tendencies* to the falling rate of profit — they push ṙ up; sV (investment in employment) pulls ṙ down — the minus sign encodes Marx's *Capital* I, ch. 25 ("the growth rate of the rate of surplus value is inversely proportional to the change in the labor force") and Glombowski's g_e = ε₀ − ε₁sV;
- **First bracket (sC − sV)**: measures the *deviation* between the propensities to invest in constant and variable capital (i.e., the movement of the organic composition of capital). Book p. 213 states that the first term "measures the effect that the **deviation** between the propensity to invest in constant and variable capital exerts on profitability."

**The economic mechanism behind the minus sign (confirmed against Glombowski's original paper)**: Glombowski's (1983) model is a Goodwin-tradition Marxian growth model: the rate of surplus value σ = M/V is "Marx's measure of the exploitation" and the indicator of income distribution; the degree of employment is β = B/A (employment relative to labour supply); labour supply grows at the exogenous rate n. The model deliberately retains the **reserve-army channel of wage formation** — his footnote criticises Okishio (real wage given), Roemer (rate of surplus value fixed), Laibman ("class struggle neutrality" = "nothing else than the abstraction from influences of the reserve army of labour on wage formation") and Furth–Heertje–van der Veen (subsistence wages under excess labour supply). The companion Glombowski–Krüger model spells the wage equation out as ŵ = −a₁ + a₂β, justified "by reference to the increase in bargaining power as employment approaches ever higher levels", and explicitly labelled Marxian. Hence: **higher employment → reserve army shrinks → workers' bargaining power rises → real wages rise → the (growth of the) rate of surplus value falls.** Goodwin (1983) himself put it this way: "it is in the region of full employment that the problem of the inverse relation of wages and profits arises since, with labour shortages, the real wage tends to rise strongly. **This is the Marxian concept of the fluctuation of the reserve army of labour**."

**If the sign of sV is changed to "+"**: the economic content becomes g_e ≈ **+sV** — the more one invests in employment, the *faster* the rate of surplus value grows. This is directly opposite to Marx (Capital I, ch. 25), to Glombowski's ε₁ > 0, and to the entire reserve-army/Phillips-curve tradition the equation belongs to; the first bracket would also turn from a "deviation" effect into a "level" effect, losing the OCC-deviation interpretation. **From the classical-Marxian framework the minus sign is correct; flipping it has no theoretical basis.**

---

## 2. Testing the four variants against the six errata findings

Variants (only (8.15) is modified; (8.16)–(8.19) unchanged):

| Variant | Eq. (8.15) | Description |
|---|---|---|
| orig | −a₁(sC−sV)r² + a₂(δ+τ−sV)r | As printed |
| A | −a₁(sC+sV)r² + a₂(δ+τ+sV)r | Both sV signs flipped |
| B | −a₁(sC+sV)r² + a₂(δ+τ−sV)r | First bracket only |
| C | −a₁(sC−sV)r² + a₂(δ+τ+sV)r | Second bracket only |

| Errata finding | orig (printed) | A (both +) | B (first only) | C (second only) |
|---|---|---|---|---|
| **F1** Is the claimed equilibrium {0.0947, 0.1188, 0.8708, 0.05, 0.045} an equilibrium? | residual max **0.02939** ✗ | 0.02939 ✗ | 0.02939 ✗ | 0.02939 ✗ |
| **F2** Does printed closed form (8.25) reproduce the claimed equilibrium? | {0.1756, 0.165, 0.6263} ✗ | ✗ | ✗ | ✗ |
| **F3** Correct closed form (analytical) | {−0.2727, −0.165, −0.4033} ✗ | {0.0857, 0.525, 1.2833} ✗ | {−0.0857, −0.525, −1.2833} ✗ | {0.2727, 0.165, 0.4033} ✗ (only sV* formula correct) |
| **F4** Existence of a positive, economically meaningful equilibrium | none ✗ | exists but **sC\*+sV\* = 1.81 > 1** ✗ | none ✗ | exists, sC\*+sV\* = 0.57 ≤ 1 **✓** |
| **F5** Jacobian at the correct equilibrium: stability / period | −0.130±0.049i, T=127y, **unstable** ✗ | +0.012±0.104i, T=**60.4y**, **unstable** ✗ | −0.092±0.044i, T=142y, unstable ✗ | +0.017±0.176i, T=**35.6y**, **unstable** ✗ |
| **F5b** Is the printed matrix the Jacobian of any variant? | max diff 1.07–1.43, none ✗ | ✗ | ✗ | ✗ |
| **F6** Nonlinear simulation period (book claims "≈50 years") | sV, sC, δ, τ ≈ 51–54y ✓ (r drifts) | 57–60y ✗ | — | 34–35y ✗ |

**Findings per item**:
- **F1 cannot be rescued by the flip**: the claimed equilibrium's maximum residual of 0.02939 comes from (8.17)'s sĊ = sC(τ−sV·r) (and (8.16)'s sV̇ = −0.0275) — equations that do not contain (8.15). No variant that modifies only (8.15) can make it zero. The claimed equilibrium {0.0947, 0.1188, 0.8708} is *not* an equilibrium for **any** variant.
- **F2 cannot be rescued**: the printed closed form (8.25) is independent of the variants and still gives {0.1756, 0.165, 0.6263}.
- **F3 cannot be rescued**: printed (8.25) is algebraically wrong for all four variants; only variant C's sV\* formula happens to be right (printed sV\* = 0.165 equals variant C's correct value).
- **F4 partially improves**: variants A and C admit a positive equilibrium. A's r\* = 0.0857 is close to the claimed 0.0947, but sV\* = 0.525 and sC\* = 1.283 are far from the claimed 0.1188/0.8708 and violate the book's own condition sC+sV ≤ 1 (book p. 222; note 17); variant C yields a fully positive, economically self-consistent equilibrium {0.2727, 0.165, 0.4033} — the positive mirror of the original system's only algebraic solution {−0.2727, −0.165, −0.4033} — but it still contradicts the claimed values.
- **F5 cannot be rescued**: for all variants, the dominant eigenvalue at the correct equilibrium has positive real part → locally **unstable**, still contradicting the book's "stable node-focus"; periods A: 60.4y, C: 35.6y, neither "≈50y"; the printed matrix matches no variant's Jacobian at any candidate point (max diff ≥ 1.0).
- **F6 not improved**: variant A's nonlinear period 57–60y (farther from 50y than the original system's 51–54y), variant C's only 34–35y (too short). The original system (printed equations + closed-form equilibrium initial values) remains the only combination reproducing an ≈50-year cycle nonlinearly.

---

## 3. Historical comparison after the flip (book Table 5.1, five waves, working-age n)

Model output = nonlinear period for each variant at its "correct equilibrium ±5% perturbation" (6T window, median of sC(t) peaks):

| Wave | Actual T | orig (printed) | A (both +) | C (second only) |
|---|---|---|---|---|
| W1 1790–1845 (n=0.537%) | 55y | 49.9y (−9%) | 55.6y (+1%) | 26.7y (−51%) |
| W2 1845–1896 (n=0.493%) | 51y | 49.0y (−4%) | 55.5y (+9%) | 26.4y (−48%) |
| W3 1896–1940 (n=0.907%) | 44y | 48.9y (+11%) | 57.0y (+30%) | 29.7y (−33%) |
| W4 1940–1982 (n=1.691%) | 42y | 145y (degenerate) | 60.6y (+44%) | 36.5y (−13%) |
| W5 1982–2022 (n=1.623%) | 40y | 146y (degenerate) | 60.2y (+51%) | 35.9y (−10%) |

- **Variant A**: systematically too long (19th century +1–9%, 20th century +44–51%), and n rising → period lengthening (55.6→60.6), i.e. the direction of the n–T relation still contradicts history (the original model's defect is not repaired);
- **Variant C**: systematically too short (19th century about half, 20th century −10–13%), with the same inverted n–T direction;
- Neither beats the printed original (which is ±10% in the 19th century and degenerates for W4/W5).

---

## 4. 2050–2080 projections after the flip

| Scenario (n, 15–64 annual growth) | Variant A: T_eig / nonlinear sC | Variant C: T_eig / nonlinear sC | Reference: original model (closed-form point) |
|---|---|---|---|
| Optimistic (UN High, +0.464%) | 56.2y / 55.4y | 26.9y / 26.1y | 33–39y |
| Medium (UN Medium, +0.048%) | 54.8y / 54.2y | 23.6y / 22.8y | 33–39y |
| Pessimistic (UN Low, −0.302%) | 53.8y / 53.2y | 20.7y / 20.0y | 33–39y |

- Variant A lengthens the sixth wave to ~54–56y (peak pushed to ~2057–58); variant C shortens it to ~20–27y (peak ~2043–45);
- All variants are locally unstable at their equilibrium, so period results are initial-condition sensitive; both directions conflict with the historical calibration (W4/W5 at 42/40y);
- The original model (printed equations + closed-form equilibrium + 5% perturbation) remains the only combination that matches the historical 40–42y.

---

## 5. Conclusions

1. **Changing the sign of sV in (8.15) from − to + does NOT resolve the six reported internal inconsistencies.**
   - F1/F2/F3 are entirely unaffected (the residual at the claimed equilibrium comes from (8.16)/(8.17), which contain no (8.15); the printed closed form is unchanged and wrong for every variant);
   - F5 (instability, period ≠ 50y, matrix ≠ Jacobian) persists for all variants;
   - The only "improvement" is that variants A/C have positive equilibria, but their values disagree with the claimed ones (A: r\* = 0.0857 ≈ 0.0947, yet sV\* = 0.525 is 4.4× the claimed 0.1188 and sC\* = 1.283 is 1.5×; C: r\* = 0.2727 is 2.9× the claimed value), A violates sC+sV ≤ 1, and C's nonlinear period is only 34–35y.
2. **Economically the flip has no basis**: the minus sign is part of the book's own derivation chain — Marx (Capital I, ch. 25: g_e inversely proportional to the change in the labor force) → Glombowski (1983) g_e = ε₀ − ε₁sV → the book's g_e ≅ −sV → (8.15). The mechanism behind the sign is the reserve-army/Phillips-curve channel formalised in the Goodwin tradition: more investment in employment tightens the labour market, strengthens workers' bargaining power, raises the real wage, and slows the growth of the rate of surplus value. Changing to "+" would assert "the more you invest in employment, the faster the rate of surplus value grows", which contradicts the classical-Marxian framework directly.
3. **The book's numerical problems have another source**: the claimed equilibrium {0.0947, 0.1188, 0.8708} satisfies no sign variant of (8.15) (even the bare requirement of (8.16), sC·r = 0.11, demands r = 0.1263, not 0.0947); the printed closed form (8.25) satisfies (8.16)/(8.18) but not (8.17) (which requires sV·r = τ, i.e. sV\* = 0.256 instead of 0.165); the printed matrix's eigenvalues are self-consistent (period 66.4y) but the matrix is no variant's Jacobian at any candidate point. The problem is the algebra of the printed closed form (8.25) and the mutual inconsistency of the claimed equilibrium/matrix/parameters — not the sign of sV in (8.15).
4. **History and projections**: the two economically admissible variants (A/C) fit history significantly worse than the printed original (A too long by 30–51%, C too short by 10–51%), and the 2050–2080 projections fall into two extremes (C: 20–27y; A: 54–56y), both unreliable. The previous conclusion stands: use the printed equations + closed-form equilibrium + nonlinear dynamics for the "≈50 year" cycle (51–54y), and under the mild post-2050 labour-force decline the most credible sixth-wave period is 33–39y.

---

## References

- Glombowski, J. (1983). A Marxian model of long run capitalist development. *Zeitschrift für Nationalökonomie*, 43(4), 363–382. DOI: 10.1007/BF01283186.
- Glombowski, J., & Krüger, M. — Goodwin-type wage equation ŵ = −a₁ + a₂β ("increase in bargaining power as employment approaches ever higher levels").
- Goodwin, R. M. (1983). On the inverse relation of wages and profits near full employment: "the Marxian concept of the fluctuation of the reserve army of labour."
- Marx, K. *Capital*, Vol. I, ch. 25 ("The General Law of Capitalist Accumulation"), pp. 653, 655.
- Sasaki, H. (2013). Cyclical growth in a Goodwin-Kalecki-Marx model. *Journal of Economics* (source of the parameter values, note 16).
- Tsoulfidis, L., & Paitaridis, D. (2019). Dominance of the organic composition of capital on the rate of profit (empirical basis of the OCC-deviation term, book p. 213).
