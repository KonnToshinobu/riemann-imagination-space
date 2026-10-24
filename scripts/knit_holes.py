"""
穴を編む（コンヌの吸収スペクトル）
  模様 f(x)  →  布 F(x) = Σ_{n≥1} f(n x)   （すべての整数倍に重ねる＝素数の手をすべて使って編む）
  布を周波数に分解（メリン変換、臨界線 s=1/2+ir 上）:
      M F(1/2+ir) = ∫ F(x) x^{-1/2+ir} dx
  事実（コンヌ 1999）:  M F(s) = ζ(s) · M f(s)   なので、どんな f でも F の分解は ζ の零点で消える。
  零点＝「どんな模様で編んでも布が存在できない周波数」＝穴。
  ※ 極 s=1 を避けるため ∫ f dx = 0 の模様を使う（コンヌの空間 S(R)_0）
"""
import numpy as np, json, matplotlib, mpmath as mp
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

zeros = np.array(json.load(open("data/zeros_300.json")))
w = 0.06
def make_pattern(c1, c2):
    # f(x) = G(log x; c1) − e^{c1−c2} G(log x; c2),  ∫ f dx = 0
    def f_u(u):
        return np.exp(-(u-c1)**2/(2*w*w)) - np.exp(c1-c2)*np.exp(-(u-c2)**2/(2*w*w))
    return f_u

patterns = [(0.0, 0.7), (0.5, -0.3), (1.0, 1.6), (-0.4, 0.3)]
U = np.linspace(-12.0, 5.0, 1701)          # u = log x
du = U[1]-U[0]
logn = np.log(np.arange(1, 2500001))

def knit(f_u):
    F = np.zeros_like(U)
    for i, u in enumerate(U):
        F[i] = np.sum(f_u(u + logn))       # Σ_n f(n x),  log(nx) = u + log n
    return F

r = np.linspace(0.5, 60, 2381)
def mellin_half(G):                        # ∫ G(e^u) e^{u/2} e^{iru} du
    E = G*np.exp(U/2)
    return np.array([np.trapezoid(E*np.exp(1j*rr*U), U) for rr in r])

results = []
for (c1, c2) in patterns:
    f_u = make_pattern(c1, c2)
    f = f_u(U); F = knit(f_u)
    Mf = mellin_half(f); MF = mellin_half(F)
    results.append((f, F, Mf, MF))
    print(f"模様 (c1={c1}, c2={c2}):  布の分解 |MF| の最小位置 r = {r[np.argmin(np.abs(MF)[r>10])+np.sum(r<=10)]:.3f}")

# 検算：ゼータを使って  MF = ζ(1/2+ir)·Mf
mp.mp.dps = 15
zeta_line = np.array([complex(mp.zeta(0.5+1j*rr)) for rr in r[::10]])
f, F, Mf, MF = results[0]
err = np.abs(MF[::10] - zeta_line*Mf[::10]).max()/np.abs(MF).max()
print("検算 |MF − ζ·Mf| / max|MF| =", err)

fig, axes = plt.subplots(3, 1, figsize=(13, 13), dpi=150)
cols = ["#1f5fbf", "#d9541e", "#2a9d4f", "#8a5fbf"]

ax = axes[0]
f, F, _, _ = results[0]
ax.plot(U, f, color="#999", lw=1.2, ls="--", label="模様 f（一枚）")
ax.plot(U, F, color=cols[0], lw=1.8, label="布 F = Σ_n f(n·x)：模様を x, 2x, 3x, … すべての整数倍に重ねたもの")
for n in [1, 2, 3, 4, 5, 6, 8]:
    ax.axvline(-np.log(n), color="#bbb", lw=0.5)
    ax.text(-np.log(n), F.max()*1.05, f"×{n}", ha="center", fontsize=7.5, color="#555")
ax.set_xlim(-6, 3); ax.set_xlabel("u = log x"); ax.set_ylabel("振幅")
ax.set_title("① 穴を編む：模様を整数倍の位置に重ねる（左へ行くほど重なりが密になる）", fontsize=10.5, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25)

ax = axes[1]
for k, ((f, F, Mf, MF), c, (c1, c2)) in enumerate(zip(results, cols, patterns)):
    ax.plot(r, np.abs(MF)**2/np.max(np.abs(MF)**2), color=c, lw=2.2-0.3*k, zorder=5-k, label=f"模様 (c1={c1}, c2={c2}) で編んだ布の分解 |MF|²")
for z in zeros[zeros < 60]:
    ax.axvline(z, color="k", lw=0.6, ls=":")
ax.set_yscale("log"); ax.set_xlim(0, 60); ax.set_ylim(1e-5, 2)
ax.set_xlabel("周波数 r   （点線＝ゼータの零点）"); ax.set_ylabel("布の光の強さ（正規化、対数目盛）")
ax.set_title("② 違う模様で編んでも、光が消える穴の位置は同じ：穴＝零点", fontsize=10.5, loc="left")
ax.legend(fontsize=8, ncol=2, loc="lower left"); ax.grid(alpha=0.25)

ax = axes[2]
f, F, Mf, MF = results[0]
ratio = np.abs(MF)/np.abs(Mf)
ax.plot(r, ratio, color="#1f5fbf", lw=2, label="布の分解 ÷ 模様の分解   |MF| / |Mf|   （模様に依らない）")
ax.plot(r[::10], np.abs(zeta_line), "o", color="#d9541e", ms=4, label="|ζ(1/2+ir)|（検算：ゼータを直接計算）")
for z in zeros[zeros < 60]:
    ax.axvline(z, color="k", lw=0.6, ls=":")
ax.set_xlim(0, 60); ax.set_ylim(0, 4)
ax.set_xlabel("周波数 r"); ax.set_ylabel("|MF| / |Mf|")
ax.set_title("③ 布と模様の比をとると、編み方の形そのものが残る：それがゼータ関数の絶対値", fontsize=10.5, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/knit_holes.png")
