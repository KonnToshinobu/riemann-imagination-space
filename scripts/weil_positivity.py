"""
ヴェイユの正値性（有限窓）
  模様の基底 φ_k(u) = exp(-(u-u_k)^2 / 2σ^2),  u_k ∈ [-L, L]
  相互相関 g_jk(u) = ∫ φ_j(v) φ_k(v-u) dv = σ√π exp(-(u-d)^2/4σ^2),  d = u_j - u_k
  そのフーリエ変換 h_jk(r) = 2πσ^2 e^{-σ^2 r^2} e^{i r d}
  明示公式で g_jk の「影」を素数だけから計算:
     W_jk = h(i/2)+h(-i/2) − g(0) log π + (1/2π)∫ h(r) Re ψ(1/4+ir/2) dr − Σ Λ(n) n^{-1/2}[g(log n)+g(−log n)]
  零点側（検算）: Z_jk = Σ_{γ>0} [h_jk(γ)+h_jk(−γ)]  … 零点が軸上なら Gram 行列なので必ず半正定値
  リーマン予想 ⇔ あらゆる L, あらゆる模様で W ≥ 0
"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.special import digamma
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

zeros = np.array(json.load(open("data/zeros_300.json")))
NMAX = 100000
lam = np.zeros(NMAX+1); sieve = np.ones(NMAX+1, bool); sieve[:2] = False
for i in range(2, int(NMAX**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
for p in np.nonzero(sieve)[0]:
    pk = p
    while pk <= NMAX: lam[pk] = np.log(p); pk *= p
ns = np.arange(2, NMAX+1); L_n = np.log(ns); wn = lam[2:]/np.sqrt(ns); keep = wn > 0
L_n, wn = L_n[keep], wn[keep]

sigma = 0.08
rgrid = np.linspace(-120, 120, 48001)
repsi = np.real(digamma(0.25 + 0.5j*rgrid))

def build(L, spacing=0.08):
    centers = np.arange(-L, L + 1e-9, spacing)
    n = len(centers)
    W = np.zeros((n, n)); Z = np.zeros((n, n))
    env = np.exp(-sigma**2*rgrid**2)
    for j in range(n):
        for k in range(j, n):
            d = centers[j] - centers[k]
            # 零点側
            Z[j, k] = Z[k, j] = 2*np.pi*sigma**2 * np.sum(np.exp(-sigma**2*zeros**2) * 2*np.cos(zeros*d))
            # 素数側（明示公式の右辺）
            pole = 2*np.pi*sigma**2*np.exp(sigma**2/4)*2*np.cosh(d/2)
            g0 = sigma*np.sqrt(np.pi)*np.exp(-d*d/(4*sigma**2))
            gam = np.trapezoid(2*np.pi*sigma**2*env*np.cos(rgrid*d)*repsi, rgrid)/(2*np.pi)
            prim = np.sum(wn*sigma*np.sqrt(np.pi)*(np.exp(-(L_n-d)**2/(4*sigma**2)) + np.exp(-(-L_n-d)**2/(4*sigma**2))))
            W[j, k] = W[k, j] = pole - g0*np.log(np.pi) + gam - prim
    return centers, W, Z

def fake_quadruplet(centers, g0=20.0, delta=0.3):
    n = len(centers); F = np.zeros((n, n))
    ts = [g0-1j*delta, g0+1j*delta, -g0-1j*delta, -g0+1j*delta]
    for j in range(n):
        for k in range(n):
            d = centers[j]-centers[k]
            F[j, k] = np.real(sum(2*np.pi*sigma**2*np.exp(-sigma**2*t*t)*np.exp(1j*t*d) for t in ts))
    return F

Ls = [1.0, 2.0, 3.0, 4.0, 5.0]
mins_W, mins_Z, mins_F, maxs = [], [], [], []
spectra = {}
for L in Ls:
    c, W, Z = build(L)
    ew = np.linalg.eigvalsh(W); ez = np.linalg.eigvalsh(Z)
    F = fake_quadruplet(c); ef_ = np.linalg.eigvalsh(Z + F)
    mins_W.append(ew.min()); mins_Z.append(ez.min()); mins_F.append(ef_.min()); maxs.append(ew.max())
    spectra[L] = ew
    print(f"L={L:3.1f}  基底 {len(c):3d}  |W−Z|max={np.abs(W-Z).max():.1e}   λ_min(素数側)={ew.min():+.3e}  λ_max={ew.max():.2e}   λ_min(偽零点入り)={ef_.min():+.3e}")

fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=150)
ax = axes[0]
for L, col in zip([2.0, 3.0, 5.0], ["#9ec3ee", "#4f8fd6", "#1f5fbf"]):
    ax.plot(np.sort(spectra[L])[::-1], "o-", ms=3, lw=1, color=col, label=f"窓の半径 L={L}")
ax.axhline(0, color="k", lw=0.8)
ax.set_yscale("symlog", linthresh=1e-3)
ax.set_xlabel("固有値の番号（大きい順）"); ax.set_ylabel("固有値（対数目盛）")
ax.set_title("素数だけから作った「影の行列」W の固有値：すべて正", fontsize=10.5, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25)

ax = axes[1]
ax.plot(Ls, mins_W, "o-", color="#1f5fbf", lw=2, label="最小固有値（素数側 W）")
ax.plot(Ls, mins_Z, "s--", color="#2a9d4f", lw=1, ms=4, label="最小固有値（零点側 Z、検算）")
ax.plot(Ls, mins_F, "^-", color="#d9541e", lw=2, label="軸外の偽零点（実部 1/2±0.3, 虚部 20）を1組混ぜた場合")
ax.axhline(0, color="k", lw=0.8)
ax.set_yscale("symlog", linthresh=1e-5)
ax.set_xlabel("窓の半径 L  （素数 n ≤ e^{2L} が参加）"); ax.set_ylabel("最小固有値")
ax.set_title("窓を広げると正値性はどこまで保たれるか", fontsize=10.5, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/weil_positivity.png")
