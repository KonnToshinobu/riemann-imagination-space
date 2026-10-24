import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import explicit_formula as ef
from math import log, pi

jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

a = ef.a
Ts = np.linspace(0.3, 7.0, 8000)
D = np.array([ef.lhs(T) - ef.rhs_archimedean(T) for T in Ts])   # 本物の零点だけの波（＝素数側）

g0, delta = 40.0, 0.3    # 偽零点の位置: 1/2 ± delta + i*g0 とその共役 (4つ組)

def fake(T, delta):
    # ρ = 1/2 + δ + iγ0 ⇒ h の引数 r = γ0 − iδ. 4つ組の寄与 = 4 Re h(γ0 − iδ)
    r = g0 - 1j*delta
    return 4*np.real(np.exp(-a*r*r)*np.cos(T*r))

F_on  = np.array([fake(T, 0.0)   for T in Ts])
F_off = np.array([fake(T, delta) for T in Ts])

fig, axes = plt.subplots(3, 1, figsize=(13, 11), dpi=150, sharex=True)
titles = [
    "① 本物の零点だけ：波は一度も 0 を上回らない（影は常に「負」）",
    f"② 軸の上に偽の零点を1組追加（実部=1/2, 虚部={g0:.0f}）：ゆらぎは入るが振幅は一定",
    f"③ 軸から外れた偽の零点を1組（4つ組）追加（実部=1/2±{delta}, 虚部={g0:.0f}）：ゆらぎが指数的に成長し、帳簿が破綻",
]
for ax, Y, F, t in zip(axes, [D, D+F_on, D+F_off], [None, F_on, F_off], titles):
    ax.axhline(0, color="k", lw=0.8)
    ax.fill_between(Ts, Y, 0, where=Y > 0, color="#d9541e", alpha=0.35, label="0 を上回った部分（素数側には存在し得ない）")
    ax.plot(Ts, Y, color="#1f5fbf", lw=1.4, label="零点側の波")
    if F is not None:
        ax.plot(Ts, F, color="#888", lw=0.8, ls=":", label="追加した偽零点の寄与だけ")
    ax.set_title(t, fontsize=10.5, loc="left")
    ax.grid(alpha=0.2)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.95)
for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]:
    for ax in axes: ax.axvline(log(p), color="#aaa", lw=0.3)
axes[0].set_ylim(-20, 6)
axes[1].set_ylim(-20, 6)
axes[-1].set_xlabel("T  （= 閉軌道の長さ log n）")
fig.suptitle("偽の零点を混ぜると何が起きるか", fontsize=13)
fig.tight_layout()
fig.savefig("figures/fake_zero.png")
print("本物: max(D) =", D.max())
print("軸上の偽零点: max =", (D+F_on).max(), "  偽零点の振幅 =", np.abs(F_on).max())
print("軸外の偽零点: max =", (D+F_off).max(), "  T=7 での偽零点の振幅 =", np.abs(F_off[-200:]).max())
