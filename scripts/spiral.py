"""
素数の螺旋：ベクトル v_n(γ) = Λ(n) n^{-1/2} w(log n) e^{-iγ log n} を頭から尾へつなぐ。
  w は log n に関するガウス窓（打ち切り）。
明示公式の双対形： Σ_ρ [e^{-b(γρ-γ)²}+e^{-b(γρ+γ)²}] = (ガンマ項) − (2/√(πb)) Σ Λ(n) n^{-1/2} e^{-(log n)²/4b} cos(γ log n)
→ 素数の螺旋和の実部は、γ が零点に近いときだけ大きく負になる。
"""
import numpy as np, matplotlib, json
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.special import digamma
import explicit_formula as ef
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

b = 0.5
ns = np.arange(2, ef.NMAX + 1)
lam = ef.lam[2:]
mask = lam > 0
ns, lam = ns[mask], lam[mask]
L = np.log(ns)
w = np.exp(-L**2/(4*b))
amp = lam/np.sqrt(ns)*w

bp = 6.0                                   # 螺旋の道用の広い窓（素数 ~2000 個が参加）
amp_p = lam/np.sqrt(ns)*np.exp(-L**2/(4*bp))
def spiral_path(g):
    v = amp_p*np.exp(-1j*g*L)
    return np.concatenate([[0], np.cumsum(v)])

def prime_side(g):
    return -(2/np.sqrt(np.pi*b))*np.sum(amp*np.cos(g*L))

def arch(g):
    # h(r)=e^{-b(r-g)^2}+e^{-b(r+g)^2}: h(i/2)+h(-i/2) - g0 log π + (1/2π)∫h Re ψ(1/4+ir/2)
    hi = 2*np.real(np.exp(-b*(0.5j-g)**2)+np.exp(-b*(0.5j+g)**2))/2
    hi = np.real(np.exp(-b*(0.5j-g)**2)+np.exp(-b*(0.5j+g)**2))*2/2
    r = np.linspace(-200, 200, 80001)
    h = np.exp(-b*(r-g)**2)+np.exp(-b*(r+g)**2)
    g0 = np.trapezoid(h, r)/(2*np.pi)
    return hi - g0*np.log(np.pi) + np.trapezoid(h*np.real(digamma(0.25+0.5j*r)), r)/(2*np.pi)

def zero_side(g):
    z = ef.zeros
    return float(np.sum(np.exp(-b*(z-g)**2)+np.exp(-b*(z+g)**2)))*2/2*1.0 + float(np.sum(np.exp(-b*(-z-g)**2)+np.exp(-b*(-z+g)**2)))*0

gs = np.linspace(5, 60, 1101)
Pz = np.array([prime_side(g)+arch(g) for g in gs])            # 素数だけ（＋ガンマ項）
Zz = np.array([2*np.sum(np.exp(-b*(ef.zeros-g)**2)+np.exp(-b*(ef.zeros+g)**2))/2*1 for g in gs])  # 零点だけ
# 共役零点込み: Σ_{γρ>0} [h(γρ)+h(-γρ)] = Σ [e^{-b(γρ-g)^2}+e^{-b(γρ+g)^2}] ×2
Zz = np.array([np.sum(np.exp(-b*(ef.zeros-g)**2)+np.exp(-b*(ef.zeros+g)**2))*2 for g in gs])
print("max |素数側−零点側| =", np.max(np.abs(Pz-Zz)))

fig = plt.figure(figsize=(13, 9.5), dpi=150)
gs_ = fig.add_gridspec(2, 3, height_ratios=[1, 1.05])
ax_top = fig.add_subplot(gs_[0, :])
ax_top.plot(gs, Pz, color="#d9541e", lw=1.8, label="素数の螺旋和（素数だけから計算）")
ax_top.plot(gs, Zz, color="#1f5fbf", lw=1.0, ls="--", label="零点の密度（零点だけから計算）")
for z in ef.zeros[ef.zeros < 60]:
    ax_top.axvline(z, color="#1f5fbf", lw=0.4, alpha=0.5)
ax_top.set_xlabel("光源の周波数 γ")
ax_top.set_ylabel("共鳴の強さ")
ax_top.set_title("素数の螺旋ベクトルを足し合わせると、零点の周波数でだけ共鳴する（縦線＝本物の零点）", fontsize=11, loc="left")
ax_top.legend(fontsize=9); ax_top.grid(alpha=0.25)

picks = [(float(ef.zeros[0]), f"γ = {ef.zeros[0]:.4f}（第1零点）"),
         (18.0, "γ = 18.0（零点ではない）"),
         (float(ef.zeros[2]), f"γ = {ef.zeros[2]:.4f}（第3零点）")]
for i, (g, title) in enumerate(picks):
    ax = fig.add_subplot(gs_[1, i])
    z = spiral_path(g)
    ax.plot(z.real, z.imag, color="#444", lw=0.35)
    sc = ax.scatter(z.real[1:], z.imag[1:], c=np.log(ns), cmap="viridis", s=3, zorder=3)
    ax.plot([0, z.real[-1]], [0, z.imag[-1]], color="#d9541e", lw=2.2, label=f"合計ベクトル |Σ| = {abs(z[-1]):.2f}")
    ax.scatter([0], [0], color="k", s=20, zorder=4)
    ax.set_aspect("equal"); ax.grid(alpha=0.25)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc="lower left")
    lim = 4.5
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
fig.text(0.5, 0.005, "下段：各素数 n（色＝log n、暗いほど小さい素数）のベクトル Λ(n) n^{-1/2} e^{-iγ log n} を頭から尾へつないだ道。零点の周波数では道が一方向へ伸び、それ以外では丸まって戻る。", ha="center", fontsize=8.5)
fig.tight_layout(rect=[0, 0.02, 1, 1])
fig.savefig("figures/prime_spiral.png")
