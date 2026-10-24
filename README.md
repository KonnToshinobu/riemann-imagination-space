# The Imagination Space — mapping an amateur's pictures onto the research frontier of the Riemann Hypothesis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23607696.svg)](https://doi.org/10.5281/zenodo.23607696) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**妄想空間 — 素人の像をリーマン予想の最前線に写像する**(Japanese version: [README_ja.md](README_ja.md))

October 15–24, 2026. A record of 22 figures and one correspondence table, produced by the back-and-forth between one person with no mathematical training beyond high school and an AI (Claude) that translated his pictures into mathematics and ran numerical experiments on them.

> **This is not a proof.** Every figure here is a numerical re-confirmation of a known theorem, a known conjecture, or a known open point. No new theorem is contained in this repository.
> Its value is a single map: *which pictures, spoken only in metaphor, correspond exactly to which parts of present-day mathematics — and where the correspondence breaks.*

If you are a student starting out in number theory or mathematical physics, this document is meant for you. It shows how far pure imagination can be carried by translation and computation, and exactly where it stops.

---

## 0. One-minute summary

- **The starting picture.** Wrap a cloth of primes around the surface of a coffee cup (a torus); light it with a candle that rotates around the cup; the shadow on the screen forms an interference pattern. The zeros are the nodes of the pattern. The line where "light and shadow balance" is Re s = 1/2.
- **What was checked.** Piece by piece, this picture corresponds to: Weil's explicit formula (the prime bands), unitarity of the dilation flow (where 1/2 comes from), Montgomery–Odlyzko pair correlation (entanglement), the Berry–Keating horizon (the count of zeros), Berry's diagonal approximation (wormhole = ramp), Hardy–Littlewood / Bogomolny–Keating (the prime-pair rhythm = plateau), Connes' absorption spectrum (knitting the holes), the Selberg trace formula (hyperbolic crochet), orthogonality of characters (white light), and Weil's proof of RH for curves (the light cone = Hodge index theorem).
- **The open core.** Weil positivity: showing, from the prime side alone and without using zeros, that "whatever pattern you weave, the shadow of the cloth laid over itself is never negative." For curves over finite fields, the light cone on the surface X×X (signature (+,−,−,…)) guarantees it. For the integers, the **surface itself** — called here the *imagination space* — has no definition yet.
- **Where this ends.** Every part of Weil's 1948 blueprint has been matched to a picture. What is missing is not assembly; it is the place to assemble on.

---

## 1. Picture ↔ mathematics

| Picture (Konno) | Mathematical object | Status | Fig. |
|---|---|---|---|
| Coffee cup (torus) | Noncommutative torus A_θ; modular surface; E×E | defined | — |
| Bands / loops of primes | Closed-orbit terms Λ(n)/√n of the explicit formula, length log p | **theorem** (Weil) | 1 |
| Rotating light source, always head-on at constant brightness | Dilation flow; unitarity ⇒ exponent 1/2 | **theorem** | 3 |
| The line where light and shadow balance | Functional equation ξ(s)=ξ(1−s); fold at Re s=1/2 | **theorem** (Tate, from the product formula) | — |
| False light (a zero off the line) | A term growing like e^{δT} | **theorem** (explicit formula) | 2 |
| Spiral | Eigenmodes x^{iγ}=e^{iγ log x}; Mellin transform | **theorem** | 4 |
| Entanglement (strong nearby, weak far away) | Pair correlation 1−(sin πr/πr)² | **conjecture** (Montgomery), numerically confirmed | 5, 19 |
| Holography / wormhole | Form-factor ramp = diagonal sum over primes | **semiclassical** (Berry 1985), numerically confirmed | 6, 15 |
| Black hole | Counting function of the xp model with a horizon | **semiclassical** (Berry–Keating), numerically confirmed | 14 |
| Density squared; two cups laid on each other | Autocorrelation → Weil's quadratic form | **theorem** (equivalent to RH) | 7 |
| Van der Waals force (sum over layers) | — | **wrong**: primes are not a sum | 8 |
| Chain stitch (pull one thread, all unravels) | Euler product (multiplicative); but the hole stays local | **half right**: damage set by the removed prime's own Euler factor | 8, 9 |
| Knitting with stitches increased each row | Hyperbolic plane; Selberg trace formula; weight ℓ/(2 sinh(ℓ/2)) | **theorem** (with the opposite sign) | 11 |
| Turning the fabric over / purl stitch / knitting the holes | Functional equation / odd cohomological degree (Deninger) / orthogonal complement (Connes) | theorem / program / theorem | 13 |
| A wormhole through the torus's center; balance | Width of M(x), ψ(x)−x is x^{1/2} | **equivalent to RH** (direction unknowable: Pólya's counterexample) | 18 |
| "Add two primes to the previous prime to get the next" | Goldbach for prime gaps; number of representations = S(h) | **conjecture**, no counterexample to 10⁷ | 18 |
| Large primes sit in the gaps left by the spirals | Sieve of Eratosthenes; Mertens; autocorrelation of the gap pattern = S(h) | **theorem** + **conjecture** (HL) | 17 |
| Colored spirals; superposition of light | Dirichlet characters; orthogonality = white light | **theorem** | 19, 20 |
| Two tori crossing perpendicularly | Fibers H, V on E×E: H²=V²=0, H·V=1 (the two null generators of the light cone) | **theorem** (curves) | 22 |
| The place where the two shadows coincide exactly | The diagonal Δ; Lefschetz fixed-point formula | **theorem** | 21 |
| They meet at points, but the points form a line | Closed orbits of length log p; trace formulas for flows (Guillemin, Connes) | **theorem** | 1, 21 |
| Projection onto the normal of the light cone | Γ'=Γ−pH−V (subtracting the pole terms); Hodge index theorem | **theorem** (curves) / **open** (integers) | 21, 22 |
| 1/2 is the single thread woven through the center | The time axis H+V = the poles of ζ at s=0,1 | theorem (curves) + numerics (Fig. 8) | 8, 22 |
| **The imagination space** | A "surface ℤ×ℤ" (candidates: Connes–Consani's arithmetic / scaling site) | **undefined** | — |

---

## 2. The 22 experiments

Figures are in `figures/`, scripts in `scripts/`, zero data in `data/`. Numbers in parentheses are the accuracy of the check.

1. **explicit_formula_orbits.png** — Weil's explicit formula. From 300 zeros alone, the closed prime orbits (length log p^k, weight Λ(n)/√n) are reconstructed (difference 1e-13).
2. **fake_zero.png** — Inserting a false zero on / off the line. An off-line zero grows like e^{δT} and breaks the books, but the prime-side background grows like e^{T/2}, faster — so a distant false zero is invisible to a naive detector. This *is* the difficulty.
3. **unitary_half.png** — (U_t f)(x)=e^{σt}f(e^t x) is unitary only for σ=1/2: the condition that the light source rotates "head-on at constant brightness."
4. **prime_spiral.png** — The sum of prime vectors Λ(n)n^{-1/2}e^{-iγ log n} resonates only at the frequencies of the zeros.
5. **zero_correlation.png** — Spacing distribution and pair correlation of 1500 zeros follow GUE (Montgomery–Odlyzko).
6. **form_factor.png** — Ramp (slope 0.89) and plateau (1.01) in the spectral form factor.
7. **weil_positivity.png** — The "shadow matrix" of band patterns, built from primes only. All eigenvalues ≥ 0 (finite-window Weil positivity); inserting an off-line false zero makes it negative, more so as the window widens.
8. **vdw_test.png** — Remove one prime's coupling / shift the exponent away from 1/2 / rescale uniformly: positivity collapses in every case. Primes are not a sum; they are one structure.
9. **chain_direction.png** — Widening the window L=2→8, the damage from removing one prime saturates and is set by that prime's own Euler-factor anti-resonance 2 log p/(√p+1) (ratio 0.90–0.96).
10. **quiet_band.png** — Explicit formula decomposed by frequency. In the zero-free band 0<g<14 the archimedean term and the prime sum cancel exactly (2.5e-4); primes fill one another's dips.
11. **hyperbolic_crochet.png** — Increasing stitches each row = hyperbolic plane. Selberg's orbit weight ℓ/(2 sinh(ℓ/2)) coincides with the prime weight Λ(n)/√n at ℓ=log p. The only difference between the two trace formulas is the sign of the orbit term.
12. **torus_curvature.png** — Elliptic-section torus: positive (outer) and negative (inner) curvature total +4π/−4π exactly (Gauss–Bonnet). Curvature sign is not the trace-formula sign.
13. **knit_holes.png** — Knitting the holes: the fabric F=Σf(nx) (pattern copied to all integer multiples), decomposed on the critical line, loses light at the same frequencies for every pattern — the zeros. MF=ζ·Mf (1e-14).
14. **horizon_count.png** — Semiclassical count of the xp model with a horizon cutoff matches the count of 1500 zeros; the remaining fluctuation S(T) has mean 0, |S|≤1.23.
15. **wormhole_from_primes.png** — The ramp built from primes alone (diagonal approximation): correlation 0.989 with the zero-side measurement for τ<1; for τ>1 the diagonal sum keeps growing.
16. **prime_pairs_plateau.png** — Prime pairs counted to 10⁷ agree with Hardy–Littlewood S(h) to 0.99–1.00. Form factor built exactly from primes (0.08%) and split into diagonal / off-diagonal: off-diagonal is zero for τ<1 and negative for τ>1, producing the plateau; replacing Λ(n)Λ(n+h) by S(h) alone reproduces it.
17. **spiral_gaps.png** — Wrap spirals of 2,3,5,7 around a cup of circumference 210: every prime >7 sits in a gap. Gap fraction follows Mertens e^{−γ}/log p. Autocorrelation of the gap pattern = S(h).
18. **center_and_dream.png** — M(x)/√x and (ψ(x)−x)/√x stay within ±0.5–0.8 (RH ⇔ width x^{1/2}). Pólya's conjecture fails at 906,150,257 (the *direction* cannot be read off). Every prime gap is a prime or a sum of two primes (Goldbach for gaps; no counterexample).
19. **family_stats.png** — Zeros of Dirichlet L-functions mod 3,4,5,7 (1292 zeros): each color is GUE; the family and ζ have the same spacing distribution.
20. **white_light.png** — Orthogonality of characters as white light: the uncolored (coherent) intensity shows the zero structure; the average over the four colors erases all interference between residue classes (1e-15). The difference — interference between primes in different classes — is where the zeros' entanglement lives. The limit of family averages (large sieve) is the "3/4 barrier."
21. **sister_world.png** — The elliptic curve y²=x³+x+1 counted over 427 primes: |a_p|≤2√p (Hasse), zeros on the circle of radius √p, intersection form on E×E always of signature (+1,−3), Gram matrix [[−2,−a_p],[−a_p,−2p]] orthogonal to the time axis negative semidefinite ⇔ 4p−a_p²≥0.
22. **light_cone.png** — Two tori crossing perpendicularly, H and V, are the two null lines of the light cone (H²=V²=0, H·V=1); time axis H+V (+2), space axis H−V (−2); Δ', Γ' are orthogonal to time and spacelike.

---

## 3. What is proved, what is open

**Theorems (re-confirmed numerically here)**
The explicit formula; unitarity and the exponent 1/2; the functional equation; orthogonality of characters; the Selberg trace formula; Gauss–Bonnet; Hasse's theorem and the Hodge index theorem for curves; Connes' trace formula (absorption spectrum); Weil's criterion (positivity ⇔ RH).

**Conjectures (strongly supported numerically, unproved)**
Montgomery pair correlation; the Hardy–Littlewood prime-pair conjecture; Goldbach; Sato–Tate (proved for elliptic curves in 2011).

**The open core**
Prove Weil positivity from the prime side without using the zeros. For curves this came from the light cone on X×X (signature (+,−,−,…)). For the integers the corresponding surface "ℤ×ℤ" — the *imagination space* of this record — is undefined. It would need three things:
1. Two tori H, V crossing perpendicularly inside it, with an intersection number H·V=1.
2. A diagonal Δ and prime dilations Γ_p, with Δ·Γ_p — the "length of the fixed closed orbit" — reproducing the prime terms Λ(n)/√n of the explicit formula.
3. An intersection form of signature (+,−,−,…): exactly one timelike direction (Hodge index).

Given 3, RH follows. Candidates: Connes–Consani's arithmetic site / scaling site (2014–).

**External context, October 2026**
OpenAI (2026-09-30, not peer reviewed, Lean-formalized): all Dirichlet L-functions and finite-order Hecke L-functions over ℚ(√−3) have no zeros with Re s>7/8 (an alternate proof: 11/12); Siegel zeros excluded. By the explicit formula this is equivalent to "the wormhole's deviation from center is at most x^{7/8}." The method is said to have a barrier at 3/4, and the cubic family over ℚ(√−3) is said to be essential (Lichtman). Figure 20 reads that barrier as the limit of white light — of averaging over a family.

---

## 4. Method (picture → translation → computation → sieve)

1. State a picture in words. No hand calculation.
2. Translate it into a known mathematical object, and extract one prediction: *if this picture is right, what should the numbers show?*
3. Compute, and separate what fits from what does not.
4. Use the misfit to trim the picture; what survives becomes material for the next one.

The pictures that failed (van der Waals as a sum, chain-stitch cascade, curvature sign as trace-formula sign) are kept on purpose. How they failed decided the next picture.

---

## 5. Reproduction

```bash
pip install mpmath numpy scipy matplotlib
cd riemann_lab
export PYTHONPATH=scripts
python3 scripts/explicit_formula.py 300      # generates data/zeros_300.json (~1 min)
python3 scripts/plot_orbits.py
python3 scripts/fake_zero.py
python3 scripts/unitary_half.py
python3 scripts/spiral.py
python3 scripts/correlation.py               # uses data/zeros_1500.json (included; ~7 min to recompute)
python3 scripts/form_factor.py
python3 scripts/weil_positivity.py
python3 scripts/vdw_test.py
python3 scripts/chain_direction.py
python3 scripts/quiet_band.py
python3 scripts/hyperbolic_crochet.py
python3 scripts/torus_curvature.py
python3 scripts/knit_holes.py                # ~3 min
python3 scripts/horizon_count.py
python3 scripts/wormhole_from_primes.py
python3 scripts/prime_pairs_plateau.py
python3 scripts/spiral_gaps.py
python3 scripts/center_and_dream.py
python3 scripts/lzeros.py 3 700              # Dirichlet L zeros (q=3,4,5,7; 15–30 min each; included)
python3 scripts/family_stats.py
python3 scripts/white_light.py
python3 scripts/sister_world.py
python3 scripts/light_cone.py
```

`CHANGELOG_ja.md` is the day-by-day log (Japanese).

---

## 6. Primary references used

- A. Weil, *Sur les "formules explicites" de la théorie des nombres premiers* (1952) — the positivity criterion
- A. Weil, *Sur les courbes algébriques et les variétés qui s'en déduisent* (1948) — RH for curves via the Hodge index theorem
- H. Montgomery, *The pair correlation of zeros of the zeta function* (1973); A. Odlyzko (1987)
- M. Berry, *Semiclassical formula for the number variance of the Riemann zeros* (1988); M. Berry & J. Keating, *H = xp and the Riemann zeros* (1999)
- E. Bogomolny & J. Keating, *Random matrix theory and the Riemann zeros II* (1996) — from Hardy–Littlewood to the plateau
- A. Connes, *Trace formula in noncommutative geometry and the zeros of the Riemann zeta function* (1999); A. Connes & C. Consani, *Weil positivity and trace formula* (2021)
- C. Deninger, *Some analogies between number theory and dynamical systems on foliated spaces* (1998)
- A. Selberg (1956); D. Taimina, *Crocheting Adventures with Hyperbolic Planes* (2009)
- J. Maldacena & L. Susskind, *Cool horizons for entangled black holes* (ER=EPR, 2013); J. Cotler et al., *Black holes and random matrices* (2017)
- G. Hardy & J. Littlewood, *Some problems of 'Partitio Numerorum' III* (1923)
- N. Katz & P. Sarnak, *Random Matrices, Frobenius Eigenvalues, and Monodromy* (1999)
- OpenAI Math Release preprints (2026-09-30 / 10-01 / 10-05); emergentmind.com/openai-math-explorer/results/003

---

## A note to young mathematicians

Nothing here replaces the training you are getting. Every step that *proved* something in this record was done by Weil, Selberg, Montgomery, Berry, Connes and others, with that training. What this record shows is something smaller and, we think, still worth knowing: a person who cannot compute can still *point*, and if someone translates and tests what he points at, the pointing lands on the right places surprisingly often — including on the exact spot where the proof is stuck. The last unsolved step is a definition, not a calculation. Definitions are made of pictures first.

*Pictures: Toshinobu Konno. Translation, numerical experiments, and text: Claude (Anthropic). 2026-10-15 to 10-24.*

## How to cite

Toshinobu Konno (2026). *riemann-imagination-space: v1.0.0 — first public release* (with numerical experiments and translation by Claude, Anthropic). Zenodo. https://doi.org/10.5281/zenodo.23607696

```bibtex
@software{konno2026imagination,
  author    = {Konno, Toshinobu},
  title     = {riemann-imagination-space: v1.0.0 --- first public release},
  year      = {2026},
  publisher = {Zenodo},
  version   = {v1.0.0},
  doi       = {10.5281/zenodo.23607696},
  url       = {https://doi.org/10.5281/zenodo.23607696}
}
```

Entry-point essays: [Medium (English)](articles/medium_en.md) · [Zenn (Japanese)](articles/zenn_ja.md)
