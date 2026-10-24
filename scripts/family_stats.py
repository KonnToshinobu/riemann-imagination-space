"""色つき螺旋（mod 3,4,5,7 の実指標）の L 関数の零点の統計 vs ゼータの零点 vs GUE"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

def unfold(z, q):
    return z/(2*np.pi)*np.log(q*z/(2*np.pi*np.e))

data = {}
for q in [3, 4, 5, 7]:
    z = np.array(json.load(open(f"data/Lzeros_q{q}.json")))
    data[q] = z[z < 440]          # 440 より上は数値精度が落ちて零点を取りこぼすので使わない
zeta = np.array(json.load(open("data/zeros_1500.json")))

def stats(u, R=3.0):
    s = np.diff(u)
    diffs = np.concatenate([(u[i+1:]-u[i])[(u[i+1:]-u[i]) < R] for i in range(len(u))])
    return s, diffs, len(u)

fig, axes = plt.subplots(2, 3, figsize=(16, 9), dpi=150)
x = np.linspace(0, 3.5, 400); r = np.linspace(0.001, 3, 500)
cols = {3: "#1f5fbf", 4: "#2a9d4f", 5: "#d9541e", 7: "#8a5fbf"}

# 上段：各 q の間隔分布
pooled_s, pooled_d, pooled_n = [], [], 0
for k, q in enumerate([3, 4, 5, 7]):
    z = data[q]; u = unfold(z, q); u = u[z > 30]
    s, d, n = stats(u)
    pooled_s.append(s); pooled_d.append(d); pooled_n += n
    ax = axes[0, k % 3] if k < 3 else None
    if ax is not None:
        ax.hist(s, bins=28, range=(0, 3.5), density=True, color=cols[q], alpha=0.5, edgecolor="white", label=f"mod {q} の色つき螺旋：零点 {n} 個")
        ax.plot(x, 32/np.pi**2*x**2*np.exp(-4*x**2/np.pi), color="k", lw=1.8, label="GUE")
        ax.plot(x, np.exp(-x), color="#888", ls="--", lw=1.2, label="独立なら")
        ax.set_xlabel("隣接零点の間隔（平均 1）"); ax.set_title(f"mod {q}（平均間隔 {s.mean():.3f}）", fontsize=10, loc="left")
        ax.legend(fontsize=8); ax.grid(alpha=0.25)

# 下段左：mod 7 の間隔
z = data[7]; u = unfold(z, 7); u = u[z > 30]; s7, d7, n7 = stats(u)
ax = axes[1, 0]
ax.hist(s7, bins=28, range=(0, 3.5), density=True, color=cols[7], alpha=0.5, edgecolor="white", label=f"mod 7 の色つき螺旋：零点 {n7} 個")
ax.plot(x, 32/np.pi**2*x**2*np.exp(-4*x**2/np.pi), color="k", lw=1.8, label="GUE")
ax.plot(x, np.exp(-x), color="#888", ls="--", lw=1.2)
ax.set_xlabel("隣接零点の間隔（平均 1）"); ax.set_title(f"mod 7（平均間隔 {s7.mean():.3f}）", fontsize=10, loc="left")
ax.legend(fontsize=8); ax.grid(alpha=0.25)

# 下段中：族をまとめた間隔 vs ゼータ
ax = axes[1, 1]
S = np.concatenate(pooled_s)
uz = zeta/(2*np.pi)*np.log(zeta/(2*np.pi*np.e)); uz = uz[zeta > 100]; sz = np.diff(uz)
ax.hist(S, bins=35, range=(0, 3.5), density=True, color="#c99a2e", alpha=0.55, edgecolor="white", label=f"色つき螺旋 4 本の族（零点 {pooled_n} 個）")
ax.hist(sz, bins=35, range=(0, 3.5), density=True, histtype="step", color="#1f5fbf", lw=2, label=f"ゼータ（色なし、零点 {len(uz)} 個）")
ax.plot(x, 32/np.pi**2*x**2*np.exp(-4*x**2/np.pi), color="k", lw=1.8, label="GUE")
ax.set_xlabel("隣接零点の間隔（平均 1）"); ax.set_title("族とゼータは同じ統計に従う", fontsize=10, loc="left")
ax.legend(fontsize=8); ax.grid(alpha=0.25)

# 下段右：対相関（族）
ax = axes[1, 2]
D = np.concatenate(pooled_d)
cnt, edges = np.histogram(D, bins=36, range=(0, 3))
ax.bar(edges[:-1], cnt/(pooled_n*(edges[1]-edges[0])), width=edges[1]-edges[0], align="edge", color="#c99a2e", alpha=0.6, edgecolor="white", label="色つき螺旋の族の対相関")
ax.plot(r, 1-(np.sin(np.pi*r)/(np.pi*r))**2, color="k", lw=1.8, label="モンゴメリー＝GUE  1−(sin πr/πr)²")
ax.axhline(1, color="#888", ls="--", lw=1)
ax.set_ylim(0, 1.4); ax.set_xlabel("零点の対の距離 r"); ax.set_title("族の対相関：近くで避け合い、遠くで無相関", fontsize=10, loc="left")
ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.25)
fig.suptitle("色つき螺旋（ディリクレ指標 mod 3,4,5,7）の L 関数の零点は、ゼータの零点と同じ「もつれ」の統計に従う", fontsize=12)
fig.tight_layout()
fig.savefig("figures/family_stats.png")
for q in [3, 4, 5, 7]:
    z = data[q]; print(f"mod {q}: 零点 {len(z)} 個, 最初の零点 {z[0]:.3f}, 最後 {z[-1]:.1f}")
