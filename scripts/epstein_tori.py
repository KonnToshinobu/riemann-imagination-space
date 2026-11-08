"""
ホッジのトーラスの影：格子点を数えたゼータ（エプシュタイン）の零点
  Q(m,n) = a m² + b mn + c n²  (判別式 D = b²−4ac < 0),  τ = (b + i√|D|)/(2a),  y = Im τ
  Z_Q(s) = Σ'_{(m,n)} Q(m,n)^{−s} = (2/√|D|)^s F(τ,s),
  F(τ,s) = Σ' y^s/|m+nτ|^{2s}
         = 2ζ(2s) y^s + 2√π Γ(s−½)/Γ(s) ζ(2s−1) y^{1−s} + (8π^s√y/Γ(s)) Σ_{n≥1} n^{s−½} σ_{1−2s}(n) K_{s−½}(2πny) cos(2πnx)
  完備化 Λ_Q(s) = (√|D|/2π)^s Γ(s) Z_Q(s) = π^{−s}Γ(s)F(τ,s) は s ↔ 1−s で対称（鏡）。
  D=−4  (x²+y²)   : Z = 4 ζ(s) L(s,χ_{−4})   ← 素数の布（オイラー積）あり → 零点は線上
  D=−20 (x²+5y²)  : 類数 2、オイラー積なし → ダヴェンポート＝ハイルブロン（1936）：線から外れた零点が無限個
"""
import json, numpy as np, mpmath as mp, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]
mp.mp.dps = 40

def divisor_sigma(n, z):
    return sum(mp.mpf(d)**z for d in range(1, n+1) if n % d == 0)

def Lambda(s, a, b, c, nmax=10):
    s = mp.mpc(s)
    D = b*b - 4*a*c; y = mp.sqrt(-D)/(2*a); x = mp.mpf(b)/(2*a)
    F = 2*mp.zeta(2*s)*y**s + 2*mp.sqrt(mp.pi)*mp.gamma(s - mp.mpf(1)/2)/mp.gamma(s)*mp.zeta(2*s - 1)*y**(1 - s)
    tail = 0
    for n in range(1, nmax+1):
        term = n**(s - mp.mpf(1)/2)*divisor_sigma(n, 1 - 2*s)*mp.besselk(s - mp.mpf(1)/2, 2*mp.pi*n*y)*mp.cos(2*mp.pi*n*x)
        tail += term
    F += 8*mp.pi**s*mp.sqrt(y)/mp.gamma(s)*tail
    return mp.pi**(-s)*mp.gamma(s)*F

forms = {"D=−4  (x²+y²)：素数の布あり": (1, 0, 1), "D=−20  (x²+5y²)：素数の布なし": (1, 0, 5)}
T = 60.0
sig = np.linspace(0.02, 0.98, 49); tt = np.linspace(0.5, T, 240)
import os
results = {}
for name, (a, b, c) in forms.items():
    cache = f"data/epstein_grid_D{b*b-4*a*c}.npy"
    if os.path.exists(cache): grid = np.load(cache)
    else:
        grid = np.array([[float(abs(Lambda(mp.mpc(sg, t), a, b, c))) for sg in sig] for t in tt]); np.save(cache, grid)
    # (1) 線上の零点：Λ(1/2+it) は実数 → 符号変化を二分法で精密化
    tf = np.arange(0.3, T, 0.05)
    vals = np.array([float(mp.re(Lambda(mp.mpc(0.5, t), a, b, c))) for t in tf])
    on = []
    for k in np.nonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]:
        z = mp.findroot(lambda t: mp.re(Lambda(mp.mpc(0.5, t), a, b, c)), (tf[k], tf[k+1]), solver="bisect", tol=1e-14)
        on.append(mp.mpc(0.5, z))
    # (2) 線から外れた零点：2 次元格子の局所最小（線から離れた列）から複素ニュートン法
    off = []
    for i_ in range(1, len(tt)-1):
        for j_ in range(1, len(sig)-1):
            if abs(sig[j_] - 0.5) < 0.06: continue
            if grid[i_, j_] <= grid[i_-1:i_+2, j_-1:j_+2].min():
                try: z = mp.findroot(lambda s: Lambda(s, a, b, c), mp.mpc(sig[j_], tt[i_]), tol=1e-16, maxsteps=60)
                except Exception: continue
                if 0 < z.real < 1 and 0.3 < z.imag < T and abs(z.real - 0.5) > 1e-6 and abs(Lambda(z, a, b, c)) < 1e-10:
                    if all(abs(z - w) > 1e-6 for w in off): off.append(z)
    # (3) 偏角の原理で総数を検算：長方形 [−0.3,1.3]×[0.3,T]
    def count_zeros():
        pts = []
        pts += [mp.mpc(x, 0.3) for x in np.linspace(-0.3, 1.3, 400)]
        pts += [mp.mpc(1.3, y) for y in np.linspace(0.3, T, 2400)]
        pts += [mp.mpc(x, T) for x in np.linspace(1.3, -0.3, 400)]
        pts += [mp.mpc(-0.3, y) for y in np.linspace(T, 0.3, 2400)]
        args = np.unwrap([float(mp.arg(Lambda(p, a, b, c))) for p in pts])
        return (args[-1] - args[0])/(2*np.pi)
    ncount = count_zeros()
    zs = sorted(on + off, key=lambda z: float(z.imag))
    results[name] = (grid, zs)
    print(f"{name}: 線上 {len(on)} 個 + 線外 {len(off)} 個 = {len(on)+len(off)}；偏角の原理による総数 = {ncount:.3f}  (0.3<t<{T:.0f})")
    for z in off: print(f"   線から外れた零点  s = {float(z.real):.6f} + {float(z.imag):.6f} i")

zeta_zeros = np.array(json.load(open("data/zeros_300.json")))
zs4 = np.array([float(z.imag) for z in results[list(forms)[0]][1]])
print("D=−4 の零点のうち ζ の零点と一致するもの:", int(np.sum([np.min(np.abs(zeta_zeros - g)) < 1e-4 for g in zs4])), "/", len(zs4), "（残りは L(s,χ_{−4}) の零点）")

# ---------- 図 ----------
fig = plt.figure(figsize=(15, 10), dpi=150, constrained_layout=True)
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.1])
for col, (name, (a, b, c)) in enumerate(forms.items()):
    grid, zs = results[name]
    # 上：格子
    ax = fig.add_subplot(gs[0, col])
    D = b*b - 4*a*c; tau = complex(b/(2*a), np.sqrt(-D)/(2*a))
    pts = np.array([[m + n*tau.real, n*tau.imag] for m in range(-6, 7) for n in range(-4, 5)])
    ax.scatter(pts[:, 0], pts[:, 1], s=14, color="#1f5fbf")
    ax.plot([0, 1, 1+tau.real, tau.real, 0], [0, 0, tau.imag, tau.imag, 0], color="#d9541e", lw=2)
    ax.set_aspect("equal"); ax.set_xlim(-4.5, 4.5); ax.set_ylim(-4.6, 4.6)
    ax.set_title(f"{name}\n周期の格子 ℤ + τℤ（橙がトーラス一枚分）", fontsize=10, loc="left"); ax.grid(alpha=0.2)
    # 下：|Λ_Q(s)| の地図と零点
    ax = fig.add_subplot(gs[1, col])
    ax.pcolormesh(sig, tt, np.log10(grid + 1e-30), cmap="Greys_r", shading="auto", rasterized=True)
    ax.axvline(0.5, color="#d9541e", lw=1.5, alpha=0.8, label="鏡 Re s = 1/2")
    on = [z for z in zs if abs(float(z.real) - 0.5) <= 1e-8]; off = [z for z in zs if abs(float(z.real) - 0.5) > 1e-8]
    ax.scatter([float(z.real) for z in on], [float(z.imag) for z in on], s=28, facecolor="none", edgecolor="#7ee0ff", lw=1.3, label=f"線上の零点 {len(on)} 個")
    if off: ax.scatter([float(z.real) for z in off], [float(z.imag) for z in off], s=60, marker="x", color="#ff5050", lw=2, label=f"線から外れた零点 {len(off)} 個")
    ax.set_xlim(0, 1); ax.set_ylim(0, T); ax.set_xlabel("Re s"); ax.set_ylabel("Im s")
    ax.set_title("格子点を数えたゼータ |Λ_Q(s)|（暗い点＝零点）", fontsize=10, loc="left"); ax.legend(fontsize=8, loc="upper right")

ax = fig.add_subplot(gs[:, 2])
for name, col_, mk in zip(forms, ["#1f5fbf", "#d9541e"], ["o", "x"]):
    zs = results[name][1]
    ax.scatter([float(z.real) for z in zs], [float(z.imag) for z in zs], s=40 if mk == "o" else 70, marker=mk, color=col_, label=name, lw=2 if mk == "x" else 1, facecolor="none" if mk == "o" else col_)
ax.axvline(0.5, color="k", lw=1)
ax.set_xlim(0, 1); ax.set_ylim(0, T); ax.set_xlabel("零点の実部"); ax.set_ylabel("零点の高さ")
ax.set_title("二つのトーラスの影の零点：素数の布があれば 1/2 に乗り、なければ外れる", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="upper right"); ax.grid(alpha=0.25)
fig.suptitle("ホッジのトーラスの影 ― 格子点を数えただけのゼータは、素数で織られていない限り、リーマン予想を破る（ダヴェンポート＝ハイルブロン 1936）", fontsize=12)
fig.savefig("figures/epstein_tori.png")
