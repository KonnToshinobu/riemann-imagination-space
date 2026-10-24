import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=150)

# (1) 編み目の増加率 → 周の長さ
ax = axes[0]
rows = np.arange(0, 13)
for q, c, lab in [(1.0, "#777", "増やさない（平面：周 ∝ 列数）"), (1.1, "#2a9d4f", "毎列 1.1 倍"), (1.2, "#d9541e", "毎列 1.2 倍"), (1.3, "#1f5fbf", "毎列 1.3 倍")]:
    circ = 6*(rows+1) if q == 1.0 else 6*q**rows
    ax.plot(rows, circ, "o-", color=c, lw=2, ms=4, label=lab)
ax.set_yscale("log")
ax.set_xlabel("列（中心からの距離 r）"); ax.set_ylabel("その列の目の数（＝周の長さ、対数目盛）")
ax.set_title("① 目を毎列 q 倍に増やすと周は e^{κr} で伸びる\n   → 曲率 −κ²（κ = log q / 列の高さ）", fontsize=10, loc="left")
ax.legend(fontsize=8.5); ax.grid(alpha=0.25)

# (2) 双曲面での測地線の重み vs 素数の重み
ax = axes[1]
ell = np.linspace(0.3, 8, 400)
ax.plot(ell, ell/(2*np.sinh(ell/2)), color="#1f5fbf", lw=2.2, label="セルバーグ：輪の重み  ℓ / (2 sinh(ℓ/2))")
ax.plot(ell, ell*np.exp(-ell/2), color="#1f5fbf", lw=1.0, ls=":", label="その漸近形  ℓ e^{−ℓ/2}")
NMAX = 3000
sieve = np.ones(NMAX+1, bool); sieve[:2] = False
for i in range(2, int(NMAX**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
ps = np.nonzero(sieve)[0]
pk, wk, kk = [], [], []
for p in ps:
    n = p; k = 1
    while n <= NMAX:
        pk.append(np.log(n)); wk.append(np.log(p)/np.sqrt(n)); kk.append(k); n *= p; k += 1
pk, wk, kk = map(np.array, (pk, wk, kk))
ax.scatter(pk[kk == 1], wk[kk == 1], s=14, color="#d9541e", zorder=3, label="素数：Λ(p)/√p = log p · p^{−1/2}  （ℓ = log p）")
ax.scatter(pk[kk > 1], wk[kk > 1], s=14, color="#c99a2e", marker="s", zorder=3, label="素数冪：Λ(p^k)/√(p^k) = log p · p^{−k/2}  （k 周）")
for p in [2, 3, 5, 7]:
    ax.annotate(str(p), (np.log(p), np.log(p)/np.sqrt(p)), textcoords="offset points", xytext=(4, 4), fontsize=8)
ax.set_xlabel("輪の長さ ℓ  （素数なら log p^k）"); ax.set_ylabel("輪一本の重み")
ax.set_title("② 双曲面の閉じた輪の重みと、素数の重みは同じ形\n   （どちらも e^{−ℓ/2}：布の量 e^{ℓ} の平方根で割る）", fontsize=10, loc="left")
ax.legend(fontsize=8); ax.grid(alpha=0.25)

# (3) 符号の違い：跡公式の輪の項
ax = axes[2]
T = np.linspace(0.3, 4.2, 2000); a = 0.0003
def bump(c): return np.exp(-(T-c)**2/(4*a))/(2*np.sqrt(np.pi*a))
sel = sum((np.log(p)/np.sqrt(n))*bump(np.log(n)) for n, p in [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11),(13,13),(16,2),(17,17),(19,19),(23,23),(25,5),(27,3),(29,29),(31,31),(32,2),(37,37),(41,41),(43,43),(47,47),(49,7),(53,53),(59,59),(61,61)])
ax.plot(T, sel, color="#1f5fbf", lw=1.6, label="双曲面（セルバーグ）：輪の項は ＋  → 零点は発光線、RH は定理")
ax.plot(T, -sel, color="#d9541e", lw=1.6, label="素数（ヴェイユ）：輪の項は −  → 零点は吸収線、RH は未解決")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("輪の長さ ℓ = log n"); ax.set_ylabel("跡公式の輪の項")
ax.set_title("③ 二つの跡公式の違いは符号一つ\n   「符号が逆のセルバーグ公式を実現する空間」が未発見", fontsize=10, loc="left")
ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.25)
ax.set_ylim(-16, 16)
fig.tight_layout()
fig.savefig("figures/hyperbolic_crochet.png")
