"""
垂直に交わる二つのトーラス ＝ 光円錐の二本の光の線
  E×E 上で  H²=0, V²=0, H·V=1   →  H, V はヌル方向（光の線）
  時間軸 H+V : (H+V)² = +2        空間軸 H−V : (H−V)² = −2
  Δ' = Δ−H−V, Γ' = Γ−pH−V は時間軸に直交し、長さは −2, −2p（空間的）
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

p, a = 101, -6                                   # 例：p=101 の a_p = −6（y²=x³+x+1 で #E=108）
G = np.array([[0, 1, 1, 1], [1, 0, 1, p], [1, 1, 0, p+1-a], [1, p, p+1-a, 0]], float)   # H, V, Δ, Γ
names = ["H", "V", "Δ", "Γ"]
def ip(x, y): return x @ G @ y
H, V, D, Gm = np.eye(4)
T = H + V; S = H - V; Dp = D - H - V; Gp = Gm - p*H - V
print("(H+V)² =", ip(T, T), " (H−V)² =", ip(S, S), " H² =", ip(H, H), " V² =", ip(V, V))
print("Δ'·(H+V) =", ip(Dp, T), " Γ'·(H+V) =", ip(Gp, T), " Δ'² =", ip(Dp, Dp), " Γ'² =", ip(Gp, Gp), " Δ'·Γ' =", ip(Dp, Gp))
print("固有値の符号:", np.sign(np.linalg.eigvalsh(G)))

fig = plt.figure(figsize=(15, 6.5), dpi=150)

# 左：垂直に交わる二つのトーラス
ax = fig.add_subplot(1, 2, 1, projection="3d")
u = np.linspace(0, 2*np.pi, 80); v = np.linspace(0, 2*np.pi, 40); U, Vv = np.meshgrid(u, v)
R, r = 1.6, 0.45
X = (R + r*np.cos(Vv))*np.cos(U); Y = (R + r*np.cos(Vv))*np.sin(U); Z = r*np.sin(Vv)
ax.plot_surface(X, Y, Z, color="#1f5fbf", alpha=0.35, linewidth=0)
ax.plot_surface(X + R, Z, Y, color="#d9541e", alpha=0.35, linewidth=0)          # 90° 回して、互いに貫通させる
ax.set_box_aspect((1.6, 1, 1)); ax.set_axis_off(); ax.view_init(elev=22, azim=-55)
ax.set_title("垂直に交わる二つのトーラス H（青）と V（橙）\n自分自身とは交わらず（H²=V²=0）、互いに交わる（H·V=1）", fontsize=10)

# 右：H, V が張るミンコフスキー平面と光円錐
ax = fig.add_subplot(1, 2, 2)
# 座標：時間 t = (x_H + x_V)/√2 方向、空間 s = (x_H − x_V)/√2 方向、計量 = 2 x_H x_V = t² − s²
lim = 2.6
ax.fill_between([-lim, lim], [lim, -lim], [lim, lim], color="#ffe9a8", alpha=0.5)
ax.fill_between([-lim, lim], [-lim, lim], [-lim, -lim], color="#ffe9a8", alpha=0.5)
ax.plot([-lim, lim], [-lim, lim], color="#1f5fbf", lw=2.2, label="光の線 H（長さ 0）")
ax.plot([-lim, lim], [lim, -lim], color="#d9541e", lw=2.2, label="光の線 V（長さ 0）")
def arrow(xh, xv, c, lab, dy=0.12):
    t = (xh + xv)/np.sqrt(2); s = (xh - xv)/np.sqrt(2)
    ax.annotate("", xy=(s, t), xytext=(0, 0), arrowprops=dict(arrowstyle="->", lw=2.2, color=c))
    ax.text(s, t + dy, lab, color=c, fontsize=9.5, ha="center")
arrow(1, 1, "#2a9d4f", "時間軸 H+V   長さ² = +2")
arrow(1, -1, "#8a5fbf", "空間軸 H−V   長さ² = −2", dy=-0.3)
arrow(1.5, 0, "#1f5fbf", "H")
arrow(0, 1.5, "#d9541e", "V")
ax.text(0, 2.15, "時間的（正）：中心＝ワームホール", ha="center", fontsize=9, color="#7a5a00")
ax.text(1.9, -0.05, "空間的（負）\nΔ', Γ' はこちら", ha="center", fontsize=9, color="#555")
ax.text(-1.9, -0.05, "空間的（負）", ha="center", fontsize=9, color="#555")
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_xlabel("空間方向 (H−V)/√2"); ax.set_ylabel("時間方向 (H+V)/√2")
ax.set_title(f"H, V が張る平面は光円錐（計量 t² − s²）\n時間軸に直交する Δ' と Γ' の長さ² は −2, −2p（空間的）→ |a_p| ≤ 2√p   （例 p={p}, a_p={a}）", fontsize=10)
ax.legend(fontsize=8.5, loc="lower right"); ax.grid(alpha=0.2)
fig.tight_layout()
fig.savefig("figures/light_cone.png")
