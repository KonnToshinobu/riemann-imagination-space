"""
零点の「もつれ」：最近接間隔分布と対相関関数。
  零点を平均間隔 1 に正規化 (unfolding): γ̃ = (γ/2π) log(γ/2πe)
  比較対象:
    独立（ポアソン）      : 間隔分布 e^{-s},           対相関 = 1 （相関なし）
    GUE（量子カオス/RMT） : 間隔分布 (32/π²) s² e^{-4s²/π}, 対相関 1 − (sin πr/πr)²
"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

z = np.array(json.load(open("data/zeros_1500.json")))
u = z/(2*np.pi)*np.log(z/(2*np.pi*np.e))        # unfolding: 平均間隔 1
u = u[z > 100]                                    # 小さい零点は漸近公式が粗いので外す
s = np.diff(u)
print("零点数", len(u), " 平均間隔", s.mean())

# 対相関
R = 3.0
diffs = []
for i in range(len(u)):
    d = u[i+1:] - u[i]
    diffs.append(d[d < R])
diffs = np.concatenate(diffs)

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), dpi=150)
x = np.linspace(0, 3.5, 400)
ax = axes[0]
ax.hist(s, bins=35, range=(0, 3.5), density=True, color="#9ec3ee", edgecolor="white", label="リーマン零点（1500個）")
ax.plot(x, np.exp(-x), color="#888", ls="--", lw=1.6, label="独立なら（ポアソン）")
ax.plot(x, 32/np.pi**2*x**2*np.exp(-4*x**2/np.pi), color="#d9541e", lw=2.2, label="量子カオス（GUE）の予言")
ax.set_xlabel("隣り合う零点の間隔 s（平均=1 に正規化）"); ax.set_ylabel("確率密度")
ax.set_title("近くの零点は互いを避ける：間隔 0 の確率がゼロ", fontsize=10.5, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25)

ax = axes[1]
ax.hist(diffs, bins=45, range=(0, R), density=True, weights=np.full(len(diffs), R/len(u)*len(diffs)/len(diffs)), color="#9ec3ee", edgecolor="white", label="リーマン零点の対相関")
# density=True の正規化を「平均1」に合わせ直す
cnt, edges = np.histogram(diffs, bins=45, range=(0, R))
ax.cla()
ax.bar(edges[:-1], cnt/(len(u)*(edges[1]-edges[0])), width=edges[1]-edges[0], align="edge", color="#9ec3ee", edgecolor="white", label="リーマン零点の対相関")
r = np.linspace(0.001, R, 500)
ax.plot(r, np.ones_like(r), color="#888", ls="--", lw=1.6, label="独立なら（相関なし）")
ax.plot(r, 1-(np.sin(np.pi*r)/(np.pi*r))**2, color="#d9541e", lw=2.2, label="モンゴメリー＝GUE の予言 1−(sin πr/πr)²")
ax.set_xlabel("二つの零点の距離 r（平均間隔=1）"); ax.set_ylabel("距離 r に零点の対を見つける確率")
ax.set_title("近くでは強く相関し、遠くでは無相関（=1）に戻る", fontsize=10.5, loc="left")
ax.set_ylim(0, 1.4)
ax.legend(fontsize=9, loc="lower right"); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/zero_correlation.png")
