import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

# 楕円断面のトーラス: 断面の横半径 a, 縦半径 b, 中心からの距離 R
R, a, b = 3.0, 1.2, 0.7
u = np.linspace(0, 2*np.pi, 400)          # 断面の角度（0=外側, π=内側）
v = np.linspace(0, 2*np.pi, 200)          # 大円の角度
U, V = np.meshgrid(u, v, indexing="ij")
X = (R + a*np.cos(U))*np.cos(V); Y = (R + a*np.cos(U))*np.sin(V); Z = b*np.sin(U)

# 断面曲線 (a cos u, b sin u) の曲率 κ、法線方向の cos で Gauss 曲率
# K = κ · cos(φ) / (R + a cos u),  φ = 断面の法線と水平面の角
den = (a*a*np.sin(u)**2 + b*b*np.cos(u)**2)**1.5
kappa = a*b/den
cosphi = b*np.cos(u)/np.sqrt(a*a*np.sin(u)**2 + b*b*np.cos(u)**2)
K = kappa*cosphi/(R + a*np.cos(u))
dA_du = np.sqrt(a*a*np.sin(u)**2 + b*b*np.cos(u)**2)*(R + a*np.cos(u))*2*np.pi   # dA = |r_u| (R+a cos u) du dv, v 積分済み
total = np.trapezoid(K*dA_du, u)
pos = np.trapezoid(np.where(K > 0, K, 0)*dA_du, u); neg = np.trapezoid(np.where(K < 0, K, 0)*dA_du, u)
print(f"外側(正)の曲率総量 = {pos:+.6f},  内側(負) = {neg:+.6f},  合計 = {total:+.2e}  (ガウス＝ボネ: 0)")

fig = plt.figure(figsize=(14, 5.2), dpi=150)
ax = fig.add_subplot(1, 2, 1, projection="3d")
Kgrid = np.broadcast_to(K[:, None], U.shape)
norm = plt.Normalize(-np.abs(K).max(), np.abs(K).max())
ax.plot_surface(X, Y, Z, facecolors=plt.cm.RdBu_r(norm(Kgrid)), rstride=4, cstride=4, linewidth=0, antialiased=False, shade=False)
ax.set_box_aspect((1, 1, 0.35)); ax.set_axis_off(); ax.view_init(elev=32, azim=35)
ax.set_title("楕円断面のトーラス：赤＝曲率が正（お椀）、青＝負（鞍）", fontsize=10.5)

ax = fig.add_subplot(1, 2, 2)
ax.plot(np.degrees(u), K, color="#1f5fbf", lw=2.2)
ax.fill_between(np.degrees(u), K, 0, where=K > 0, color="#d9541e", alpha=0.35, label=f"正の部分（外側）：総量 +{pos:.3f}")
ax.fill_between(np.degrees(u), K, 0, where=K < 0, color="#1f5fbf", alpha=0.35, label=f"負の部分（内側）：総量 {neg:.3f}")
ax.axhline(0, color="k", lw=0.8)
for ang, lab in [(0, "外側"), (90, "一番上"), (180, "内側（穴）"), (270, "一番下"), (360, "外側")]:
    ax.axvline(ang, color="#999", lw=0.6); ax.text(ang, ax.get_ylim()[1]*0.95 if False else K.max()*1.05, lab, ha="center", fontsize=8.5)
ax.set_xlabel("断面上の角度（度）"); ax.set_ylabel("ガウス曲率 K")
ax.set_title(f"正と負は厳密に釣り合う（合計 {total:+.1e}）：ガウス＝ボネ ∫K dA = 2πχ = 0", fontsize=10.5, loc="left")
ax.legend(fontsize=9, loc="lower right"); ax.grid(alpha=0.25)
ax.set_xlim(0, 360)
fig.tight_layout()
fig.savefig("figures/torus_curvature.png")
