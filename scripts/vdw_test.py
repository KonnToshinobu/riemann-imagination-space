"""
ファンデルワールス仮説の検証：素数の「結合の手」をいじって、影の行列 W の正値性がどう変わるか
  (a) 手を外す：素数 p の寄与 (p, p^2, p^3, ...) を丸ごと除く
  (b) 手の強さを変える：Λ(n) n^{-1/2} → Λ(n) n^{-σ}  (σ ≠ 1/2)
  (c) 手の強さを一律に α 倍
"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.special import digamma
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

NMAX = 100000
lam = np.zeros(NMAX+1); sieve = np.ones(NMAX+1, bool); sieve[:2] = False
for i in range(2, int(NMAX**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
base_p = np.zeros(NMAX+1, int)
for p in np.nonzero(sieve)[0]:
    pk = p
    while pk <= NMAX: lam[pk] = np.log(p); base_p[pk] = p; pk *= p
ns = np.arange(2, NMAX+1); keep = lam[2:] > 0
L_n = np.log(ns[keep]); lam_n = lam[2:][keep]; p_n = base_p[2:][keep]; n_n = ns[keep]

sigma = 0.08
rgrid = np.linspace(-120, 120, 48001)
repsi = np.real(digamma(0.25 + 0.5j*rgrid))
env = np.exp(-sigma**2*rgrid**2)

def W_matrix(L, weights, spacing=0.08):
    centers = np.arange(-L, L + 1e-9, spacing)
    n = len(centers); W = np.zeros((n, n))
    cache = {}
    for j in range(n):
        for k in range(j, n):
            d = round(centers[j]-centers[k], 6)
            if d not in cache:
                pole = 2*np.pi*sigma**2*np.exp(sigma**2/4)*2*np.cosh(d/2)
                g0 = sigma*np.sqrt(np.pi)*np.exp(-d*d/(4*sigma**2))
                gam = np.trapezoid(2*np.pi*sigma**2*env*np.cos(rgrid*d)*repsi, rgrid)/(2*np.pi)
                prim = np.sum(weights*sigma*np.sqrt(np.pi)*(np.exp(-(L_n-d)**2/(4*sigma**2)) + np.exp(-(-L_n-d)**2/(4*sigma**2))))
                cache[d] = pole - g0*np.log(np.pi) + gam - prim
            W[j, k] = W[k, j] = cache[d]
    return np.linalg.eigvalsh(W)

L = 4.0
w_true = lam_n/np.sqrt(n_n)

results = []
def run(label, w):
    ev = W_matrix(L, w); results.append((label, ev.min(), ev.max()))
    print(f"{label:38s} λ_min = {ev.min():+.3e}   λ_max = {ev.max():.3e}")

run("本物 (Λ(n)/√n)", w_true)
for p in [2, 3, 5, 7]:
    run(f"(a) 素数 {p} の手を外す", np.where(p_n == p, 0.0, w_true))
run("(a) 素数 2 と 3 の手を外す", np.where((p_n == 2)|(p_n == 3), 0.0, w_true))
run("(a) 素数 2,3,5,7 の手を外す", np.where(p_n <= 7, 0.0, w_true))
run("(a) 11 以上の手を全部外す", np.where(p_n > 7, 0.0, w_true))
run("(a) 素数の手を全部外す（無限遠のみ）", 0*w_true)
for s in [0.40, 0.45, 0.48, 0.52, 0.55, 0.60]:
    run(f"(b) 手の強さ Λ(n) n^(-{s:.2f})", lam_n/n_n**s)
for a in [0.9, 1.1]:
    run(f"(c) 手の強さを一律 {a} 倍", a*w_true)

labels = [r[0] for r in results]; mins = np.array([r[1] for r in results])
fig, ax = plt.subplots(figsize=(12, 7), dpi=150)
colors = ["#1f5fbf"] + ["#2a9d4f"]*8 + ["#d9541e"]*6 + ["#8a5fbf"]*2
y = np.arange(len(labels))[::-1]
vals = np.abs(mins)
ax.barh(y, vals, color=colors)
ax.set_xscale("log"); ax.set_xlim(1e-16, 1e3)
ax.axvline(1e-12, color="k", lw=0.8, ls=":"); ax.text(1.3e-12, y[0]+0.6, "←ここより左は数値誤差（＝正値性が保たれている）", fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
for yi, m in zip(y, mins):
    ax.text(abs(m)*1.4, yi, f"{m:+.1e}", va="center", ha="left", fontsize=8)
ax.set_xlabel("最小固有値の大きさ |λ_min|（対数目盛）   本物以外はすべて負＝正値性が壊れる")
ax.set_title(f"素数の「結合の手」をいじったとき、重ねた布の影は正のままか（窓 L={L}, 帯幅 σ={sigma}）", fontsize=11, loc="left")
ax.grid(alpha=0.25, axis="x")
fig.tight_layout()
fig.savefig("figures/vdw_test.png")
