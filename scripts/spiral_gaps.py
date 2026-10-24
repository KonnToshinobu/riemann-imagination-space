"""
螺旋を巻いて残る隙間（エラトステネスの篩）
 ① 2,3,5,7 の螺旋を周長 210 のカップに巻く。覆われない隙間に 7 より大きい素数がすべて座る。
 ② 螺旋を増やすほど隙間は狭まる：隙間の割合 Π(1−1/p) ≈ e^{−γ}/log y（メルテンス）→ 大きな素数は稀
 ③ 隙間模様の自己相関（r と r+h が共に隙間である割合）が、昨日数えた素数対のリズム S(h) に一致する
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

fig = plt.figure(figsize=(15, 10), dpi=150)
gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1])

# ---------- ① 周長 210 のカップに 2,3,5,7 の螺旋 ----------
ax = fig.add_subplot(gs[0, 0])
M = 210; sp = [2, 3, 5, 7]; cols = ["#1f5fbf", "#2a9d4f", "#d9541e", "#8a5fbf"]
theta = 2*np.pi*np.arange(M)/M
covered = np.zeros(M, int)
for k, (p, c) in enumerate(zip(sp, cols)):
    r = 1.0 - 0.09*k
    idx = np.arange(0, M, p)
    ax.scatter(r*np.cos(theta[idx]), r*np.sin(theta[idx]), s=9, color=c, label=f"{p} の螺旋（{p} の倍数を覆う）", zorder=3)
    covered[idx] = 1
gaps = np.nonzero(covered == 0)[0]
r = 0.60
ax.scatter(r*np.cos(theta[gaps]), r*np.sin(theta[gaps]), s=12, facecolors="none", edgecolors="#444", label=f"隙間（{len(gaps)} 個 / 210）", zorder=3)
isp = np.ones(M, bool); isp[:2] = False
for i in range(2, 15):
    if isp[i]: isp[i*i::i] = False
big = [n for n in range(11, M) if isp[n]]
ax.scatter(r*np.cos(theta[big]), r*np.sin(theta[big]), s=14, color="#c99a2e", label=f"7 より大きい素数（{len(big)} 個）：すべて隙間に座る", zorder=4)
ax.set_aspect("equal"); ax.set_axis_off()
ax.set_title("① 周長 210 のカップに 2,3,5,7 の螺旋を巻く：残った隙間にしか大きな素数は座れない", fontsize=10, loc="left")
ax.legend(fontsize=8, loc="center", frameon=False)
print("① 隙間", len(gaps), "個, 7より大きい素数", len(big), "個, 素数で隙間でないもの:", [n for n in big if covered[n]])

# ---------- ② 隙間の割合 ----------
ax = fig.add_subplot(gs[0, 1])
primes = [p for p in range(2, 200) if all(p % q for q in range(2, int(p**0.5)+1))]
frac = np.cumprod([1 - 1/p for p in primes])
ax.plot(primes, frac, "o-", color="#1f5fbf", ms=4, lw=1.5, label="螺旋を p まで巻いた後の隙間の割合  Π(1−1/p)")
ax.plot(primes, np.exp(-0.5772)/np.log(primes), "--", color="#d9541e", lw=1.8, label="メルテンス：e^{−γ}/log p")
ax.set_xscale("log"); ax.set_xlabel("巻いた螺旋の最大の素数 p"); ax.set_ylabel("隙間の割合")
ax.set_title("② 螺旋を巻くほど隙間は狭まる：大きな素数が稀になる理由（密度 ∝ 1/log）", fontsize=10, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25)

# ---------- ③ 隙間模様の自己相関 vs 素数対のリズム ----------
Msv = 2*3*5*7*11*13          # 30030
gap = np.ones(Msv, bool)
for p in [2, 3, 5, 7, 11, 13]: gap[::p] = False
H = 60
auto = np.array([np.mean(gap & np.roll(gap, -h)) for h in range(1, H+1)])/np.mean(gap)**2

NA = 10_000_000
isp = np.ones(NA+1, bool); isp[:2] = False
for i in range(2, int(NA**0.5)+1):
    if isp[i]: isp[i*i::i] = False
counts = np.array([np.count_nonzero(isp[:-h] & isp[h:]) for h in range(1, H+1)])

ax = fig.add_subplot(gs[1, :])
hs = np.arange(1, H+1)
ax.bar(hs, counts/counts[1], color="#9ec3ee", edgecolor="white", label="数えた素数対 (p, p+h) の個数（p ≤ 1000万、h=2 を 1）")
ax.plot(hs, auto/auto[1], "o", color="#d9541e", ms=6, label="2,3,5,7,11,13 の螺旋が残す隙間模様の自己相関（r と r+h が共に隙間）")
ax.set_xlabel("距離 h"); ax.set_ylabel("現れやすさ（h=2 を 1）")
ax.set_title("③ 素数対のリズムは、螺旋の隙間模様の自己相関そのもの：もつれの正体は「隙間の形」", fontsize=10, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25, axis="y"); ax.set_xlim(0, H+1)
even = hs % 2 == 0
print("③ 偶数 h での 観測/隙間模型:", np.round((counts/counts[1])[even]/(auto/auto[1])[even], 3))
fig.tight_layout()
fig.savefig("figures/spiral_gaps.png")
