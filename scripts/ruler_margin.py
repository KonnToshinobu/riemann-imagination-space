"""
定規の余裕：ロウソクは穴にどれだけの余裕で勝っているか
  模様 f(u) = Σ_k a_k φ_k(u),  φ_k = ガウス帯（幅 σ）, 中心 u_k ∈ [−L, L]
  ヴェイユの定規   Q_W(f) = Σ_ρ |f̂(γ)|²·(重み)  ＝ ロウソク項 Q_C(f) − 穴（素数）項 Q_P(f)   （明示公式）
  普通の定規       ‖f‖² = ∫ |f|² du
  余裕(L) = min_f Q_W(f)/‖f‖²  （一般化固有値の最小。基底の取り方に依らない）
  W, C, P, G はすべて d = u_j − u_k だけの関数（テプリッツ）なので、w(d) を一度計算すれば足りる。
"""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.special import digamma
from scipy.linalg import eigh, toeplitz
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

sigma, spacing = 0.08, 0.08
Ls = [1, 2, 3, 4, 5, 6, 7, 8]
NMAX = int(np.exp(2*max(Ls) + 15*sigma)) + 1
sieve = np.ones(NMAX+1, bool); sieve[:2] = False
for i in range(2, int(NMAX**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0]; del sieve
Ln, wn = [np.log(primes)], [np.log(primes)/np.sqrt(primes)]
for p in primes[primes <= int(NMAX**0.5)]:
    pk = p*p
    while pk <= NMAX:
        Ln.append([np.log(pk)]); wn.append([np.log(p)/np.sqrt(pk)]); pk *= p
Ln = np.concatenate(Ln); wn = np.concatenate(wn); o = np.argsort(Ln); Ln, wn = Ln[o], wn[o]
print(f"素数冪 {len(Ln):,} 個 (n ≤ {NMAX:,})")

rgrid = np.linspace(-120, 120, 48001); repsi = np.real(digamma(0.25 + 0.5j*rgrid))
env = np.exp(-sigma**2*rgrid**2)
zeros = np.array(json.load(open("data/zeros_1500.json")))

def parts(d):
    pole = 2*np.pi*sigma**2*np.exp(sigma**2/4)*2*np.cosh(d/2)
    g0 = sigma*np.sqrt(np.pi)*np.exp(-d*d/(4*sigma**2))
    gam = np.trapezoid(2*np.pi*sigma**2*env*np.cos(rgrid*d)*repsi, rgrid)/(2*np.pi)
    lo = np.searchsorted(Ln, abs(d) - 14*sigma); hi = np.searchsorted(Ln, abs(d) + 14*sigma)
    prim = np.sum(wn[lo:hi]*sigma*np.sqrt(np.pi)*(np.exp(-(Ln[lo:hi]-d)**2/(4*sigma**2)) + np.exp(-(-Ln[lo:hi]-d)**2/(4*sigma**2))))
    zside = 2*np.pi*sigma**2*np.sum(np.exp(-sigma**2*zeros**2)*2*np.cos(zeros*d))
    return pole - g0*np.log(np.pi) + gam, prim, g0, zside     # ロウソク C, 穴 P, 普通の定規 G, 零点側 Z

mmax = int(round(2*max(Ls)/spacing))
ds = np.arange(0, mmax+1)*spacing
C, P, G, Z = np.array([parts(d) for d in ds]).T
W = C - P


res = []; spectra = {}; fmin = {}
for L in Ls:
    n = int(round(2*L/spacing)) + 1
    TW, TC, TP, TG, TZ = (toeplitz(v[:n]) for v in (W, C, P, G, Z))
    ew, vw = eigh(TW, TG); ec = eigh(TC, TG, eigvals_only=True); ep = eigh(-TP, TG, eigvals_only=True); ez = eigh(TZ, TG, eigvals_only=True)
    thr = 1e-7*ew.max()
    kern = int(np.sum(np.abs(ew) < thr)); kern_z = int(np.sum(np.abs(ez) < thr))
    # 重み e^{-σ²γ²} が 1e-7 を越える零点の個数（cos, sin で 2 本ずつ拘束）
    neff = int(np.sum(np.exp(-sigma**2*zeros**2) > 1e-7))
    res.append((L, n, ew.min(), ew.max(), ec.min(), ep.min(), ez.min(), kern, kern_z, neff))
    spectra[L] = ew; fmin[L] = vw[:, 0]
    print(f"L={L}  基底 {n:3d}  余裕 min Q_W/‖f‖² = {ew.min():+.2e}（零点側 {ez.min():+.2e}）  max {ew.max():.2f}   "
          f"ロウソクのみ min {ec.min():+.2e}  穴のみ min {ep.min():+.2e}   長さ≈0 の模様の数 {kern}（零点側 {kern_z}） / 基底 {n}  ↔ 基底 − 2×有効零点数 = {n-2*neff}")
res = np.array(res)

# ---------- 図 ----------
fig = plt.figure(figsize=(15, 9.5), dpi=150, constrained_layout=True)
gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1])
ax = fig.add_subplot(gs[0, 0])
for L, col in zip([2, 4, 6, 8], ["#9ec3ee", "#6aa3e0", "#3b7fd0", "#1f5fbf"]):
    ax.plot(np.arange(1, len(spectra[L])+1), np.sort(spectra[L])[::-1], "o-", ms=3, lw=1, color=col, label=f"L={L}（基底 {len(spectra[L])}）")
ax.axhline(0, color="k", lw=0.8); ax.set_yscale("symlog", linthresh=1e-8)
ax.set_xlabel("固有値の番号（大きい順）"); ax.set_ylabel("長さ比 Q_W(f)/‖f‖²（symlog）")
ax.set_title("ヴェイユの定規で測った長さ比のスペクトル：正の塊のあと、崖の下は「長さゼロ」（1e-8 以下）", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5); ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[1, 0])
ax.plot(res[:, 0], res[:, 1], "s-", color="#999", lw=1.5, label="基底の数（模様の自由度）")
ax.plot(res[:, 0], res[:, 7], "o-", color="#1f5fbf", lw=2.4, ms=7, label="長さ≈0 に測られる模様の数（素数から作った定規）")
ax.plot(res[:, 0], res[:, 8], "x", color="#2a9d4f", ms=9, mew=2, label="同（零点 1500 個から作った定規）")
ax.plot(res[:, 0], res[:, 1] - 2*res[:, 9], "--", color="#d9541e", lw=1.5, label="予言：基底 − 2×(重みが 1e-7 を越える零点の数)")
ax.set_xlabel("窓の半径 L"); ax.set_ylabel("個数")
ax.set_title("定規に見えない模様は、窓を広げた分だけ増える（見えるのは零点の周波数の成分だけ）", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="upper left"); ax.grid(alpha=0.25)

L0 = Ls[-1]; n = len(fmin[L0]); centers = np.linspace(-L0, L0, n)
ug = np.linspace(-L0-0.5, L0+0.5, 4000)
f = np.sum(fmin[L0][:, None]*np.exp(-(ug[None, :]-centers[:, None])**2/(2*sigma**2)), axis=0)
f /= np.abs(f).max()
ax = fig.add_subplot(gs[0, 1])
ax.plot(ug, f, color="#1f5fbf", lw=1)
ax.set_xlabel("u = log x"); ax.set_ylabel("f(u)"); ax.set_title(f"定規が「長さゼロ」と測る模様の一つ（L={L0}）", fontsize=10.5, loc="left"); ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[1, 1])
gam_grid = np.linspace(0, 60, 6000)
fh = np.abs(np.trapezoid(f[None, :]*np.exp(-1j*gam_grid[:, None]*ug[None, :]), ug, axis=1))**2
ax.plot(gam_grid, fh/fh.max(), color="#1f5fbf", lw=1.2, label="|F(γ)|²  （模様 f の周波数成分）")
for g in zeros[zeros < 60]: ax.axvline(g, color="#d9541e", lw=0.8, alpha=0.8)
ax.plot([], [], color="#d9541e", lw=0.8, label="零点 γ")
ax.set_yscale("log"); ax.set_ylim(1e-10, 2); ax.set_xlabel("周波数 γ"); ax.set_ylabel("|F(γ)|²（対数）")
ax.set_title("その模様は零点の周波数をすべて避けている（谷が零点に一致）", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="lower left"); ax.grid(alpha=0.25)
fig.suptitle("定規の余裕を測ったら、余裕はゼロだった ― ヴェイユの定規は零点の周波数しか見ない（半正定値、決して正定値ではない）", fontsize=12.5)
fig.savefig("figures/ruler_margin.png")
