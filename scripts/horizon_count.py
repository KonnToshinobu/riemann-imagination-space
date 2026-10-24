"""
地平線が零点を作る：ベリー＝キーティングの xp 模型
  古典軌道 x p = E は双曲線（地平線なしでは無限に伸び、スペクトルは連続）
  位相空間に切り落とし |x| ≥ l_x, |p| ≥ l_p（地平線）を入れると、軌道 xp=E が囲む面積は有限:
      A(E) = E log(E/(l_x l_p)) − E     （l_x l_p = 2π ととる）
  準位の個数 = 面積 / 2π  （半古典量子化）  →  N_BK(E) = (E/2π)(log(E/2π) − 1) + 7/8
  本物の零点の個数（リーマン＝フォン・マンゴルト）: N(T) = (T/2π) log(T/2πe) + 7/8 + S(T)
  同じ式。違いは S(T)（零点の「ゆらぎ」）だけ。
"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

zeros = np.array(json.load(open("data/zeros_1500.json")))
T = np.linspace(2, zeros[-1], 20000)
N_true = np.searchsorted(zeros, T)                      # 階段：T 以下の零点の個数
N_BK = T/(2*np.pi)*(np.log(T/(2*np.pi)) - 1) + 7/8      # 地平線つき xp の半古典計数
S = N_true - N_BK

fig = plt.figure(figsize=(14, 9), dpi=150)
gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1])

# 位相空間の絵
ax = fig.add_subplot(gs[0, 0])
lx = lp = np.sqrt(2*np.pi)
x = np.linspace(lx, 30, 400)
for E, c in [(20, "#9ec3ee"), (60, "#4f8fd6"), (120, "#1f5fbf")]:
    xx = x[x*lp <= E];
    ax.plot(xx, E/xx, color=c, lw=2, label=f"軌道 xp = E  (E={E})")
    ax.fill_between(xx, lp, E/xx, color=c, alpha=0.12)
ax.axvline(lx, color="#d9541e", lw=2.2); ax.axhline(lp, color="#d9541e", lw=2.2)
ax.text(lx+0.4, 27, "地平線 x = l_x", color="#d9541e", fontsize=9)
ax.text(18, lp+0.6, "地平線 p = l_p", color="#d9541e", fontsize=9)
ax.set_xlim(0, 30); ax.set_ylim(0, 30); ax.set_xlabel("x"); ax.set_ylabel("p")
ax.set_title("① xp 模型に地平線（切り落とし）を入れる\n   囲まれた面積 ÷ 2π ＝ 準位の個数", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="upper right"); ax.grid(alpha=0.25)

# 計数の比較
ax = fig.add_subplot(gs[0, 1])
ax.step(T, N_true, where="post", color="#1f5fbf", lw=1.2, label="本物：T 以下のゼータ零点の個数（階段）")
ax.plot(T, N_BK, color="#d9541e", lw=1.8, ls="--", label="地平線つき xp の半古典計数  (T/2π)(log(T/2π) − 1) + 7/8")
ax.set_xlim(0, 200); ax.set_ylim(0, 80)
ax.set_xlabel("T（零点の虚部）"); ax.set_ylabel("個数")
ax.set_title("② 地平線が作る準位の個数と、本物の零点の個数", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="upper left"); ax.grid(alpha=0.25)
axin = ax.inset_axes([0.55, 0.08, 0.42, 0.42])
axin.step(T, N_true, where="post", color="#1f5fbf", lw=1.0); axin.plot(T, N_BK, color="#d9541e", lw=1.4, ls="--")
axin.set_xlim(1800, 1981); axin.set_ylim(1300, 1500); axin.tick_params(labelsize=7); axin.set_title("T≈1900 付近", fontsize=8)

# ゆらぎ S(T)
ax = fig.add_subplot(gs[1, :])
ax.plot(T, S, color="#1f5fbf", lw=0.8)
ax.axhline(0, color="k", lw=0.6)
ax.set_xlim(0, zeros[-1]); ax.set_ylim(-1.5, 1.5)
ax.set_xlabel("T"); ax.set_ylabel("S(T) = 本物 − 半古典")
ax.set_title(f"③ 残りは零点の「ゆらぎ」S(T) だけ：平均 {S.mean():+.3f}、標準偏差 {S.std():.3f}、|S| の最大 {np.abs(S).max():.2f}  （地平線が個数の本体を決め、ゆらぎが素数の対相関を担う）", fontsize=10, loc="left")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/horizon_count.png")
print("S(T): mean", S.mean(), "std", S.std(), "max|S|", np.abs(S).max())
