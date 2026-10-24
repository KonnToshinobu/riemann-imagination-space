"""
静かな帯域 0<r<14 での釣り合い：周波数ごとの分解
  試験関数 h(r) = e^{-b(r-g)^2}+e^{-b(r+g)^2}  （周波数 g のまわりの窓、幅 1/√b）
  零点側  D(g) = Σ_{γ>0} [h(γ)+h(−γ)]                … 零点の密度（必ず ≥ 0）
  無限遠  A(g) = h(i/2)+h(−i/2) − ĝ(0) log π + (1/2π)∫ h Re ψ(1/4+ir/2) dr
  素数    P_p(g) = −(2/√(πb)) Σ_k (log p) p^{−k/2} e^{−(k log p)²/4b} cos(g k log p)
  明示公式:  D = A + Σ_p P_p
"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.special import digamma
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

zeros = np.array(json.load(open("data/zeros_1500.json")))
b = 2.0
gs = np.linspace(0, 40, 2001)

D = np.array([np.sum(np.exp(-b*(zeros-g)**2) + np.exp(-b*(zeros+g)**2))*2 for g in gs])

r = np.linspace(-80, 80, 64001); repsi = np.real(digamma(0.25+0.5j*r))
def A(g):
    h = np.exp(-b*(r-g)**2) + np.exp(-b*(r+g)**2)
    pole = 2*np.real(np.exp(-b*(0.5j-g)**2) + np.exp(-b*(0.5j+g)**2))
    g0 = np.trapezoid(h, r)/(2*np.pi)
    return pole - g0*np.log(np.pi) + np.trapezoid(h*repsi, r)/(2*np.pi)
Aa = np.array([A(g) for g in gs])

NMAX = 20000
sieve = np.ones(NMAX+1, bool); sieve[:2] = False
for i in range(2, int(NMAX**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0]
def P_p(p):
    out = np.zeros_like(gs); pk = p; k = 1
    while pk <= NMAX:
        u = k*np.log(p)
        out += -(2/np.sqrt(np.pi*b))*np.log(p)*p**(-k/2)*np.exp(-u*u/(4*b))*np.cos(gs*u)
        pk *= p; k += 1
    return out
Pp = {p: P_p(p) for p in primes}
Ptot = sum(Pp.values())
print("明示公式の検算 max|D − (A+ΣP)| =", np.abs(D - (Aa+Ptot)).max())

fig, axes = plt.subplots(3, 1, figsize=(13, 13), dpi=150, sharex=False)

ax = axes[0]
ax.plot(gs, D, color="#1f5fbf", lw=2.4, label="零点の密度 D(g)（零点だけから）")
ax.plot(gs, Aa, color="#777", lw=1.6, label="無限遠の光源 A(g)（ガンマ項＋極）")
ax.plot(gs, Ptot, color="#d9541e", lw=1.6, label="素数の合計 ΣP_p(g)（素数だけから）")
ax.plot(gs, Aa+Ptot, color="k", lw=1.0, ls="--", label="A + ΣP（= D になるはず）")
ax.axvspan(0, zeros[0], color="#ffe9a8", alpha=0.5)
ax.text(zeros[0]/2, ax.get_ylim()[1]*0.0+2.6, "静かな帯域\n（零点なし）", ha="center", fontsize=9, color="#7a5a00")
ax.axhline(0, color="k", lw=0.6)
ax.set_xlim(0, 40); ax.set_ylim(-3, 3)
ax.set_ylabel("振幅"); ax.set_title("① 明示公式を周波数で分解：零点の密度 ＝ 無限遠の光源 ＋ 素数の合計", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="lower right", ncol=2); ax.grid(alpha=0.25)

ax = axes[1]
cols = ["#1f5fbf", "#2a9d4f", "#d9541e", "#8a5fbf", "#c99a2e", "#555", "#17becf"]
for p, c in zip([2, 3, 5, 7, 11, 13, 17], cols):
    ax.plot(gs, Pp[p], color=c, lw=1.6, label=f"P_{p}(g)")
    for m in [1, 3]:
        x = m*np.pi/np.log(p)
        if x < 20: ax.plot([x], [np.interp(x, gs, Pp[p])], "v", color=c, ms=7)
ax.axvspan(0, zeros[0], color="#ffe9a8", alpha=0.5)
ax.axhline(0, color="k", lw=0.6)
ax.set_xlim(0, 20)
ax.set_ylabel("振幅"); ax.set_title("② 素数一つずつの寄与 P_p(g)：▼＝反共鳴 g = π/log p, 3π/log p（p^{-ig} = −1 の周波数）", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, ncol=4, loc="lower right"); ax.grid(alpha=0.25)

ax = axes[2]
cum = Aa.copy()
shades = plt.cm.Oranges(np.linspace(0.3, 1.0, 7))
ax.plot(gs, Aa, color="#777", lw=1.6, label="無限遠だけ A")
for P, c in zip([2, 3, 5, 7, 13, 100, 20000], shades):
    cum = Aa + sum(Pp[p] for p in primes if p <= P)
    ax.plot(gs, cum, color=c, lw=1.4, label=f"A + 素数 ≤ {P}")
ax.plot(gs, D, color="#1f5fbf", lw=2.4, ls="--", label="零点の密度 D")
ax.axvspan(0, zeros[0], color="#ffe9a8", alpha=0.5)
ax.axhline(0, color="k", lw=0.6)
ax.set_xlim(0, 20); ax.set_ylim(-2.5, 2.5)
ax.set_xlabel("周波数 g  （零点の虚部と同じ目盛）"); ax.set_ylabel("振幅")
ax.set_title("③ 素数を一つずつ足していくと、静かな帯域で無限遠の光源が打ち消されて 0 に沈む", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, ncol=3, loc="lower right"); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/quiet_band.png")

band = gs < zeros[0]-1
print("静かな帯域での |A+ΣP| の最大 =", np.abs((Aa+Ptot)[band]).max(), "  |A| の最大 =", np.abs(Aa[band]).max())
