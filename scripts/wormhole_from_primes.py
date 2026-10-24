"""
素数だけからワームホール（ランプ）を組み立てる — ベリーの対角近似 (1985)
  零点から測る形状因子:  K(τ) = |Σ w_n e^{2πiτ u_n}|² / Σ w_n²     （u_n: 平均間隔1に正規化した零点）
  素数から組み立てる:     K_diag(τ) = (1/Lr²) Σ_n (Λ(n)²/n) · δ_Δ(τ − log n / Lr),   Lr = log(T/2π)
     （各素数冪 n = p^k の閉軌道が、時間 τ = log n / Lr に Λ(n)²/n の重みで寄与。同じ軌道同士の干渉＝対角）
  素数定理 Σ_{n≤x} Λ(n)²/n ≈ (log x)²/2 から  K_diag(τ) → τ   （ランプ）
  τ>1 のプラトーは、異なる素数同士の干渉（非対角・ハーディ＝リトルウッド素数対）がないと出ない。
"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

# --- 零点側（form_factor.py と同じ） ---
z = np.array(json.load(open("data/zeros_1500.json"))); z = z[z > 600]
u = z/(2*np.pi)*np.log(z/(2*np.pi*np.e)); N = len(u); w = np.hanning(N)
taus = np.linspace(0.0, 2.0, 2001)
K = np.array([abs(np.sum(w*np.exp(2j*np.pi*t*u)))**2 for t in taus])/np.sum(w**2)
def smooth(y, n):
    k = np.ones(n)/n; return np.convolve(y, k, mode="same")
Ks = smooth(K, 121)
Lr = np.log(np.sqrt(z.min()*z.max())/(2*np.pi))

# --- 素数側 ---
NMAX = int(np.exp(2.0*Lr)) + 10            # τ=2 まで必要な n
sieve = np.ones(NMAX+1, bool); sieve[:2] = False
for i in range(2, int(NMAX**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0]
ln, wt = [], []
for p in primes:
    pk = p; lp = np.log(p)
    while pk <= NMAX:
        ln.append(np.log(pk)); wt.append(lp*lp/pk); pk *= p
ln = np.array(ln); wt = np.array(wt); tau_n = ln/Lr

Delta = 0.06                                 # 零点側の移動平均幅と同じ
def K_diag(t):
    return np.sum(wt*np.exp(-(t-tau_n)**2/(2*Delta**2)))/(np.sqrt(2*np.pi)*Delta)/Lr**2
Kd = np.array([K_diag(t) for t in taus])
Delta_raw = 0.004
Kd_raw = np.array([np.sum(wt*np.exp(-(t-tau_n)**2/(2*Delta_raw**2)))/(np.sqrt(2*np.pi)*Delta_raw)/Lr**2 for t in taus])

m = (taus > 0.15) & (taus < 0.85)
print(f"零点 {N} 個, Lr={Lr:.2f}, 素数冪 {len(ln)} 個 (n≤{NMAX})")
print("ランプ領域 0.15<τ<0.85 :  <K_zeros/τ> =", np.mean(Ks[m]/taus[m]), "  <K_primes/τ> =", np.mean(Kd[m]/taus[m]))
print("ランプ領域での相関係数 corr(K_zeros, K_primes) =", np.corrcoef(Ks[m], Kd[m])[0, 1])
m2 = (taus > 1.2) & (taus < 2.0)
print("プラトー領域 1.2<τ<2 :  <K_zeros> =", np.mean(Ks[m2]), "  <K_primes(対角)> =", np.mean(Kd[m2]))

fig, axes = plt.subplots(2, 1, figsize=(13, 9.5), dpi=150, sharex=True)
ax = axes[0]
ax.plot(taus, K, color="#9ec3ee", lw=0.5, alpha=0.6, label="零点から測った形状因子（生）")
ax.plot(taus, Ks, color="#1f5fbf", lw=2.4, label="零点から測った形状因子（移動平均）")
ax.plot(taus, Kd, color="#d9541e", lw=2.2, ls="--", label="素数だけから組み立てたワームホール  (1/Lr²) Σ Λ(n)²/n  δ(τ − log n/Lr)")
ax.plot(taus, np.minimum(taus, 1), color="#888", lw=1, ls=":", label="ランダム行列の予言")
ax.axvline(1, color="k", lw=0.6)
ax.set_ylim(0, 2.0); ax.set_ylabel("K(τ)")
ax.set_title("素数だけでワームホール（ランプ）を編めるか：τ<1 で零点の測定と一致、τ>1 では素数の対角和は伸び続ける", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="upper left"); ax.grid(alpha=0.25)
ax.text(1.45, 1.75, "プラトー：異なる素数同士の干渉\n（素数対の相関）がないと出ない", fontsize=9, color="#7a3a00", ha="center")

ax = axes[1]
ax.plot(taus, K, color="#9ec3ee", lw=0.6, label="零点から測った形状因子（生）")
ax.plot(taus, Kd_raw*0.3, color="#d9541e", lw=0.9, label="素数の閉軌道の位置と重み Λ(n)²/n （縮尺 0.3）")
for n, lab in [(2,"2"),(3,"3"),(4,"2²"),(5,"5"),(7,"7"),(8,"2³"),(9,"3²"),(11,"11"),(13,"13"),(16,"2⁴"),(17,"17"),(19,"19"),(23,"23"),(25,"5²"),(29,"29"),(31,"31")]:
    ax.text(np.log(n)/Lr, 1.85, lab, ha="center", fontsize=7, color="#444")
ax.set_xlim(0, 2.0); ax.set_ylim(0, 2.0)
ax.set_xlabel("τ   （平均間隔=1 の単位の時間）"); ax.set_ylabel("K(τ)")
ax.set_title("小さい τ では個々の素数の軌道が、零点の形状因子の生データのピークにそのまま対応する", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="upper right"); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/wormhole_from_primes.png")
