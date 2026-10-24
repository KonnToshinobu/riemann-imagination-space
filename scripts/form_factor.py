"""
スペクトル形状因子（spectral form factor）
  K(τ) = |Σ_n w_n e^{2πi τ u_n}|² / Σ w_n²      u_n: 平均間隔1に正規化した零点, w_n: ハン窓
ランダム行列(GUE)の予言（連結部分）: K(τ) = τ (τ<1, ランプ),  1 (τ≥1, プラトー)
ベリー(1985): ランプは素数の閉軌道の寄与 Σ Λ(n)²/n を「対角近似」で足し合わせたものに等しい
"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

z = np.array(json.load(open("data/zeros_1500.json")))
lo = 600.0
z = z[z > lo]
u = z/(2*np.pi)*np.log(z/(2*np.pi*np.e))
N = len(u)
w = np.hanning(N)                      # 端の効果を抑える窓
taus = np.linspace(0.0, 2.0, 4001)
K = np.array([abs(np.sum(w*np.exp(2j*np.pi*t*u)))**2 for t in taus]) / np.sum(w**2)

# τ方向の移動平均（1サンプルの形状因子は指数分布で揺らぐので平滑化が必要）
def smooth(y, n):
    k = np.ones(n)/n
    return np.convolve(y, k, mode="same")
Ks = smooth(K, 121)   # 幅 0.06

rho_mid = np.log(np.sqrt(z.min()*z.max())/(2*np.pi))      # 平均密度 (1/2π)log(γ/2π) の代表値 ×2π
print(f"零点 {N} 個 (γ∈[{z.min():.0f},{z.max():.0f}]),  log(γ/2π) ≈ {rho_mid:.2f}")

fig, ax = plt.subplots(figsize=(13, 5.4), dpi=150)
ax.plot(taus, K, color="#9ec3ee", lw=0.6, alpha=0.7, label="形状因子 K(τ)（零点だけから計算、生データ）")
ax.plot(taus, Ks, color="#1f5fbf", lw=2.2, label="K(τ) の移動平均（幅 0.06）")
ax.plot(taus, np.minimum(taus, 1.0), color="#d9541e", lw=2.2, ls="--", label="ランダム行列(GUE)の予言：ランプ τ → プラトー 1")
for n, lab in [(2,"2"),(3,"3"),(4,"2²"),(5,"5"),(7,"7"),(8,"2³"),(9,"3²"),(11,"11"),(13,"13")]:
    t = np.log(n)/rho_mid
    ax.axvline(t, color="#666", lw=0.7, alpha=0.7, ymax=0.62)
    ax.text(t, 1.26, lab, ha="center", va="bottom", fontsize=8, color="#333")
ax.text(0.02, 1.38, "縦線＝素数の閉軌道の位置  τ = log n / log(γ/2π)   （生データのピークと一致）", fontsize=8.5, color="#333")
ax.set_xlim(0, 2.0); ax.set_ylim(0, 2.0)
ax.set_xlabel("τ  （時間。平均間隔=1 の単位でのフーリエ共役変数）")
ax.set_ylabel("K(τ)")
ax.set_title(f"零点 {N} 個の形状因子：ランプ（ワームホール）とプラトーは見えるか", fontsize=11, loc="left")
ax.legend(loc="upper left", fontsize=9, framealpha=0.95)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/form_factor.png")

# 数値的な要約：ランプ領域 (0.2<τ<0.8) での K(τ)/τ の平均、プラトー領域 (1.2<τ<2) の平均
m1 = (taus > 0.2) & (taus < 0.8); m2 = (taus > 1.2) & (taus < 2.0)
print("ランプ領域  <K(τ)/τ> =", np.mean(Ks[m1]/taus[m1]), " (GUEなら 1)")
print("プラトー領域 <K(τ)>  =", np.mean(Ks[m2]), " (GUEなら 1)")
