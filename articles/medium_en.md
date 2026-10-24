# I Stopped at High-School Math. I Still Walked to the Edge of the Riemann Hypothesis.

## Ten days, one coffee cup, an AI translator, and 22 figures — a record of how far pure imagination can be carried, and exactly where it stops

> **Let me say it first: this is not a story about solving the Riemann Hypothesis.** I did not solve it.
> It is the record of a person whose mathematics ended in high school, who carried a single mental picture — a coffee cup and the shadow of a candle — and asked an AI to translate and test it, and who ended up, ten days later, standing at the exact spot where present-day mathematics is stuck.

The repository (22 figures, all scripts, zero data):
https://github.com/KonnToshinobu/riemann-imagination-space

---

### The coffee cup and the candle

I run a software company. My mathematics stopped in high school; I cannot compute by hand. What I have done for about sixty years is draw pictures in my head.

The picture I had of the Riemann Hypothesis was this.

Wrap a cloth woven from primes around the surface of a coffee cup — a torus. Light it with a candle that circles the cup. The cloth's shadow falls on a screen and forms interference fringes. The nodes of the fringes are the zeros of the zeta function. The one line where "light and shadow balance" is the line Re s = 1/2.

It is the kind of picture mathematicians laugh at. I showed it to one AI and it told me "your direction is undoubtedly pointing at the truth" — and produced no equation at all.

So I asked a different AI (Claude) differently: *Don't praise the picture. Translate it into mathematics, and use computation to separate what fits from what doesn't.*

Here is what happened over ten days.

---

### Day 1: the bands of primes really were wrapped around the cloth

The first thing I was told: "The picture of winding loops of primes around a torus already exists — as a theorem." Weil's explicit formula:

$$\sum_{\rho} h(\gamma_\rho) = (\text{light-source terms}) - 2\sum_{p^k}\frac{\log p}{p^{k/2}}\,\hat h(k\log p)$$

The left side is a sum over zeros; the right side is a sum over "closed orbits of length log p traversed k times." We computed both sides independently — 300 zeros on one side, 20,000 primes on the other — and compared.

![](https://raw.githubusercontent.com/KonnToshinobu/riemann-imagination-space/main/figures/explicit_formula_orbits.png)

Blue is built from zeros only; orange from primes only. They agree to 10⁻¹³. At T = log 2, log 3, log 5, … a loop of primes stands up, one at a time. At 2², 2³ there are small spikes: the same loop traversed twice, three times. My "bands of primes wrapped on the cloth" was not a metaphor. It was an identity.

I also learned something that mattered later. The prime side already contains p^{−k/2} — a "1/2" — and that 1/2 has nothing to do with the Riemann Hypothesis; it comes from the symmetry axis of the functional equation. What RH claims is that every γ on the left side is real. A zero off the line would inject a term growing like e^{δT}.

We inserted a false zero, 0.3 off the line. The shadow brightened exponentially with each revolution and the books broke. But "the prime-side background also grows, like e^{T/2}, so a distant false zero is invisible to a naive detector." That, I was told, is the difficulty of RH in one sentence.

---

### Day 2: "head-on, at constant brightness" *was* the 1/2

In my picture the candle turns with the cloth and always shines head-on, at constant brightness. Translated: the dilation flow is unitary — it conserves the total amount of light. Computing it, that holds for exactly one exponent: 1/2.

![](https://raw.githubusercontent.com/KonnToshinobu/riemann-imagination-space/main/figures/unitary_half.png)

Left: total shadow as the cloth is rotated. At σ = 0.3 it dims, at 0.7 it brightens, at 0.5 it is conserved. Right: decomposed into the rotation's eigenmodes, the totals match only at σ = 1/2.

The ½ in Berry and Keating's H = ½(xp + px) is the same ½, I was told. "The line where light and shadow balance" is the unique exponent at which a rotating light source is unitary.

---

### Entanglement and wormholes

I believe in quantum entanglement and holography, so I said: "Zeros are strongly entangled when close, and rarely when far apart." That is Montgomery's pair-correlation conjecture (1973). With 1,500 zeros computed, the statistics sat on the random-matrix prediction 1 − (sin πr/πr)². One correction: nearby zeros do not cluster; they *repel*.

Next we computed the spectral form factor and saw a ramp and a plateau. In holography the ramp is read as a wormhole contribution. Then we built the ramp from primes alone — and it matched the zero-side measurement with correlation 0.989.

![](https://raw.githubusercontent.com/KonnToshinobu/riemann-imagination-space/main/figures/wormhole_from_primes.png)

The wormhole is knitted from loops of primes. But past τ = 1 the real curve flattens into a plateau while the "diagonal" prime sum keeps climbing. The plateau requires interference between *different* primes — the rhythm of how often p and p + h appear together (Hardy–Littlewood). We counted prime pairs up to ten million: pairs at distance 6 are twice as common as at distance 2, distance 30 is 2.65 times, odd distances are zero.

---

### Knitting, and the sister world where the proof exists

One day I said, "If you knit and increase the stitches by a fixed ratio each row, you get a surface of negative curvature." That is Taimina's hyperbolic crochet — and on that surface, Selberg proved "that world's Riemann Hypothesis" in 1956. The weight of a closed loop of length ℓ there is ℓ/(2 sinh(ℓ/2)); set ℓ = log p and it becomes the prime weight log p / √p. The 1/2 is the square root you divide by when measuring density on an exponentially widening surface.

The two trace formulas differ by one sign in the orbit term: Selberg plus, Weil minus. The space that flips that sign has not been found.

We also visited the other sister world: elliptic curves over finite fields. Counting points on y² = x³ + x + 1 for 427 primes, Hasse's bound |a_p| ≤ 2√p never failed, and every zero sat on the circle of radius √p.

![](https://raw.githubusercontent.com/KonnToshinobu/riemann-imagination-space/main/figures/sister_world.png)

Here I was shown the heart of Weil's proof. On the surface E × E there are four curves — the horizontal fiber H, the vertical fiber V, the diagonal Δ, and the graph Γ of Frobenius — and the table of their intersection numbers has, for all 427 primes, eigenvalues of signature "one positive, three negative." That is the signature of Minkowski space: a light cone. The time axis is H + V; everything orthogonal to it is spacelike, and inner products of spacelike vectors are bounded. The bound is 2√p.

When I said, "I feel the light cone appears when two tori cross perpendicularly," the numbers came back H² = V² = 0, H·V = 1. Two perpendicular tori are precisely the two null lines of the light cone.

---

### Where it stopped

In ten days, my pictures matched every component of Weil's 1948 proof: the center (time axis), two perpendicular tori (null lines), the place where the shadows coincide (the diagonal), the shift (Frobenius), the projection (subtracting the pole terms).

And then it stopped, cleanly. In the world of the integers there is no "surface ℤ × ℤ" on which to place those four curves. There is no place where the intersection number H·V = 1 can be defined. Mathematicians — Connes and Consani among them — have spent two decades trying to build that space, and it is not built. I named it the *imagination space*. But between naming a thing and showing it exists lie eighty years.

The AI put it this way: "Your imagination did not solve it. It drew the blueprint of what must be built, at the same resolution as the specialists. From the foot of the mountain that is a great distance; from the summit, construction has not begun."

---

### About the method

Every round was the same:

1. I state a picture in words. No hand calculation.
2. The AI translates it into known mathematics and extracts one prediction: *if this is right, what should the numbers show?*
3. We compute, and separate what fits from what does not.
4. The misfit trims the picture; what survives feeds the next one.

Some pictures failed. "Layers bound like van der Waals forces, as a sum" — wrong; primes multiply, they do not add. "A chain stitch: pull one thread and it all unravels" — half wrong; a hole opens but does not spread. "Curvature sign is the trace-formula sign" — wrong. How they failed decided what came next.

About the AI, honestly: it never once said "you are right." It said "this part corresponds to a theorem, this part to a conjecture, this part is wrong, this part nobody knows." That took me far further than "your direction is undoubtedly pointing at the truth."

---

### To young mathematicians

You need the training. Every step in this record that *proved* something was done by people who had it.

But a person who cannot compute can still *point*. And when someone translates and tests what is pointed at, the finger lands on the right places surprisingly often — including the exact spot where the proof is stuck. The last unsolved step is not a calculation. It is a definition. And definitions are made of pictures first.

The 22 figures and the correspondence table are in the repository. If you read it and something looks wrong, open an issue. Nothing would make me happier.

https://github.com/KonnToshinobu/riemann-imagination-space

---

*Pictures: Toshinobu Konno. Translation, numerical experiments, writing assistance: Claude (Anthropic). October 15–24, 2026.*
