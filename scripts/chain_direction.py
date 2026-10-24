"""
「ほどける方向」：窓 L を広げながら、素数 p の手を一本抜いたときの正値性の壊れ方 |λ_min| を比べる
"""
import numpy as np, matplotlib, time
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.special import digamma
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

Ls = [2, 3, 4, 5, 6, 7, 8]
sigma = 0.08
NMAX = int(np.exp(2*max(Ls) + 15*sigma)) + 1
t0 = time.time()
sieve = np.ones(NMAX+1, bool); sieve[:2] = False
for i in range(2, int(NMAX**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0]; del sieve
Ln, wn, pn = [], [], []
for p in primes:
    lp = np.log(p); pk = p
    while pk <= NMAX:
        Ln.append(np.log(pk)); wn.append(lp/np.sqrt(pk)); pn.append(p); pk *= p
Ln = np.array(Ln); wn = np.array(wn); pn = np.array(pn)
order = np.argsort(Ln); Ln, wn, pn = Ln[order], wn[order], pn[order]
print(f"素数冪 {len(Ln)} 個 (n ≤ {NMAX:,}), 準備 {time.time()-t0:.1f}s")

rgrid = np.linspace(-120, 120, 48001)
repsi = np.real(digamma(0.25 + 0.5j*rgrid)); env = np.exp(-sigma**2*rgrid**2)

def W_eigmin(L, mask, spacing=0.08):
    centers = np.arange(-L, L + 1e-9, spacing); n = len(centers)
    ds = np.round(centers[:, None] - centers[None, :], 6)
    uniq = np.unique(ds)
    vals = {}
    w = np.where(mask, wn, 0.0)
    for d in uniq:
        pole = 2*np.pi*sigma**2*np.exp(sigma**2/4)*2*np.cosh(d/2)
        g0 = sigma*np.sqrt(np.pi)*np.exp(-d*d/(4*sigma**2))
        gam = np.trapezoid(2*np.pi*sigma**2*env*np.cos(rgrid*d)*repsi, rgrid)/(2*np.pi)
        lo = np.searchsorted(Ln, abs(d) - 14*sigma); hi = np.searchsorted(Ln, abs(d) + 14*sigma)
        prim = np.sum(w[lo:hi]*sigma*np.sqrt(np.pi)*(np.exp(-(Ln[lo:hi]-d)**2/(4*sigma**2)) + np.exp(-(-Ln[lo:hi]-d)**2/(4*sigma**2))))
        vals[d] = pole - g0*np.log(np.pi) + gam - prim
    W = np.vectorize(vals.get)(ds)
    return np.linalg.eigvalsh(W).min()

test_primes = [2, 3, 5, 7, 11, 13]
res = {p: [] for p in test_primes}; base = []
for L in Ls:
    t0 = time.time()
    base.append(W_eigmin(L, np.ones_like(wn, bool)))
    line = f"L={L}  本物 λ_min={base[-1]:+.1e}  |"
    for p in test_primes:
        m = W_eigmin(L, pn != p); res[p].append(m)
        line += f"  抜{p}: {m:+.3f}"
    print(line + f"   ({time.time()-t0:.0f}s)")

fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=150)
cols = ["#1f5fbf", "#2a9d4f", "#d9541e", "#8a5fbf", "#c99a2e", "#555"]
ax = axes[0]
for p, c in zip(test_primes, cols):
    ax.plot(Ls, -np.array(res[p]), "o-", color=c, lw=2, label=f"素数 {p} の手を抜く")
ax.set_xlabel("窓の半径 L   （参加する素数 n ≤ e^{2L}）"); ax.set_ylabel("壊れ方の大きさ  −λ_min")
ax.set_title("手を一本抜いたとき、窓を広げると壊れ方はどう育つか", fontsize=10.5, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25)
ax = axes[1]
pp = np.array(test_primes); sat = -np.array([res[p][-1] for p in test_primes])
pred = 2*np.pi*sigma**2/0.08 * 2*np.log(pp)/(np.sqrt(pp)+1)   # /0.08 = 帯の間隔（行列の目盛）
pgrid = np.arange(2, 60); ax.plot(pgrid, 2*np.pi*sigma**2/0.08*2*np.log(pgrid)/(np.sqrt(pgrid)+1), "-", color="#999", lw=1.5, label="予言：オイラー因子の反共鳴の深さ ∝ 2 log p/(√p+1)")
ax.plot(pp, sat, "o", color="#d9541e", ms=8, label=f"実測（L={Ls[-1]} での壊れ方）")
ax.set_xlabel("抜いた素数 p"); ax.set_ylabel("壊れ方の大きさ  −λ_min")
ax.set_title("壊れ方は窓に依らず、抜いた素数自身のオイラー因子で決まる", fontsize=10.5, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25)
print("実測/予言 =", np.round(sat/pred, 3))
fig.tight_layout()
fig.savefig("figures/chain_direction.png")
