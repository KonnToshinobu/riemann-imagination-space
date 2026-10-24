"""
「回転する光源」の数学的実体：拡大縮小の流れ  (U_t f)(x) = e^{σ t} f(e^t x)
布（関数 f）を速度 t で回したとき、影の総量 ||U_t f||^2 = ∫|U_t f|^2 dx が保存されるのは
指数 σ がいくつのときか。
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

x = np.linspace(1e-6, 60, 400000)
f = np.exp(-(np.log(x) - 1.0)**2 / 0.3) / np.sqrt(x)          # 布：トーラス上の何かの模様（適当な滑らかな関数）
norm0 = np.trapezoid(np.abs(f)**2, x)

ts = np.linspace(-2.0, 2.0, 41)
fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), dpi=150)

for sigma, col in [(0.3, "#d9541e"), (0.5, "#1f5fbf"), (0.7, "#2a9d4f")]:
    bright = []
    for t in ts:
        Uf = np.exp(sigma*t) * np.interp(np.exp(t)*x, x, f, left=0, right=0)
        bright.append(np.trapezoid(np.abs(Uf)**2, x) / norm0)
    axes[0].plot(ts, bright, color=col, lw=2, label=f"σ = {sigma}")
    axes[0].plot(ts, np.exp((2*sigma-1)*ts), color=col, lw=0.8, ls=":")
axes[0].axhline(1, color="k", lw=0.6)
axes[0].set_yscale("log")
axes[0].set_xlabel("回転量 t")
axes[0].set_ylabel("影の総量（回転前を1とする）")
axes[0].set_title("光源を回したとき影の総量は保たれるか：  (U_t f)(x) = e^{σt} f(e^t x)", fontsize=10, loc="left")
axes[0].legend()
axes[0].grid(alpha=0.25)

# メリン変換：この回転の「固有振動」は s = σ + iγ でパラメータ付けされる。ユニタリな回転 σ=1/2 のとき
# ∫|f|^2 dx = (1/2π)∫|F(1/2+iγ)|^2 dγ  (プランシュレル) が成り立つ。
gam = np.linspace(-15, 15, 601)
for sigma, col in [(0.3, "#d9541e"), (0.5, "#1f5fbf"), (0.7, "#2a9d4f")]:
    F = np.array([np.trapezoid(f * x**(sigma + 1j*g - 1), x) for g in gam])
    tot = np.trapezoid(np.abs(F)**2, gam) / (2*np.pi) / norm0
    axes[1].plot(gam, np.abs(F)**2, color=col, lw=2, label=f"σ = {sigma}：固有振動の総量 / 影の総量 = {tot:.3f}")
axes[1].set_xlabel("固有振動の周波数 γ   （零点の虚部に対応）")
axes[1].set_ylabel("|メリン変換 F(σ+iγ)|²")
axes[1].set_title("回転の固有振動に分解したとき、総量が一致するのは σ=1/2 だけ", fontsize=10, loc="left")
axes[1].legend(fontsize=8.5)
axes[1].grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/unitary_half.png")
