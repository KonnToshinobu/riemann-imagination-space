"""
妄想1：素数のワームホールはトーラスの中心を通り、素数とそれ以外が均衡している
   → 素数と「素数でない数」の均衡＝メビウス関数 μ(n) の和 M(x) = Σ μ(n)
      （μ = +1: 偶数個の異なる素数の積, −1: 奇数個, 0: 平方因子あり）
      素数定理 ⇔ M(x)/x → 0（中心を通る）
      リーマン予想 ⇔ |M(x)| ≤ x^{1/2+ε}（中心からのずれが √x の幅に収まる）  ← 1/2 はここにも現れる
      同じく ψ(x) = Σ_{n≤x} Λ(n) について、RH ⇔ |ψ(x) − x| ≤ √x log²x
妄想2：次の素数は、前の素数に素数を足して現れる
   (a) 一つの素数を足す： p_{k+1} − p_k が素数か
   (b) 二つの素数を足す： p_{k+1} − p_k = q + q'（ゴールドバッハ）
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

N = 10_000_000
# 線形篩で μ(n), Λ(n) を同時に
mu = np.ones(N+1, np.int8); lp = np.zeros(N+1, np.int32)
isp = np.ones(N+1, bool); isp[:2] = False
for i in range(2, int(N**0.5)+1):
    if isp[i]:
        isp[i*i::i] = False
primes = np.nonzero(isp)[0]
# μ: 平方因子で 0、素因数の個数の偶奇で ±1
mu = np.ones(N+1, np.int8)
for p in primes:
    mu[p::p] *= -1
    pp = p*p
    if pp <= N: mu[pp::pp] = 0
mu[0] = 0
M = np.cumsum(mu.astype(np.int64))
x = np.arange(N+1)
# ψ(x) − x
lam = np.zeros(N+1)
for p in primes:
    pk = p
    while pk <= N: lam[pk] = np.log(p); pk *= p
psi = np.cumsum(lam) - x

fig = plt.figure(figsize=(15, 11), dpi=150)
gs = fig.add_gridspec(3, 2)

ax = fig.add_subplot(gs[0, :])
xs = x[1000:]
ax.plot(xs, M[1000:]/np.sqrt(xs), color="#1f5fbf", lw=0.5, label="M(x)/√x  ：素数の偶奇の均衡（メビウス）の、中心からのずれ")
ax.plot(xs, psi[1000:]/np.sqrt(xs), color="#d9541e", lw=0.5, alpha=0.8, label="(ψ(x) − x)/√x  ：素数の個数の、期待値からのずれ")
ax.axhline(0, color="k", lw=0.8); ax.axhline(1, color="#888", lw=0.8, ls="--"); ax.axhline(-1, color="#888", lw=0.8, ls="--")
ax.set_xscale("log"); ax.set_ylim(-2, 2)
ax.set_xlabel("x"); ax.set_ylabel("ずれ ÷ √x")
ax.set_title(f"妄想1：ワームホールは中心を通るか  — 1000万まで、ずれは √x の幅に収まる（max |M|/√x = {np.max(np.abs(M[1000:])/np.sqrt(xs)):.3f}, max |ψ−x|/√x = {np.max(np.abs(psi[1000:])/np.sqrt(xs)):.3f}）  RH ⇔ この幅が x^{{1/2}} で止まる", fontsize=9.5, loc="left")
ax.legend(fontsize=8.5, loc="lower left"); ax.grid(alpha=0.25)

# ---------- 妄想2 ----------
gaps = np.diff(primes)                   # p_{k+1} − p_k
gap_is_prime = isp[gaps]                 # (a)
frac_a = np.mean(gap_is_prime)
# (b) ゴールドバッハ：各偶数 h の表し方の数 r(h) = #{(q,q'): q ≤ q', q+q' = h}
H = int(gaps.max())
small = primes[primes <= H]
r = np.zeros(H+1, np.int64)
for q in small:
    s = q + small[small >= q]
    s = s[s <= H]
    np.add.at(r, s, 1)
r_gap = r[gaps]
print(f"素数 {len(primes)} 個 (≤ {N:,}), 隙間の最大 {H}")
bad_a = np.nonzero(~gap_is_prime & (gaps > 1))[0]
print(f"(a) 隙間自体が素数（＝一つの素数を足す）: {frac_a*100:.1f}%   最初の反例 {primes[bad_a[0]]} → {primes[bad_a[0]+1]} (隙間 {gaps[bad_a[0]]})")
ok_b = gap_is_prime | (r_gap > 0)
print(f"(b) 隙間が素数、または二つの素数の和（ゴールドバッハ）: 反例 {np.sum(~ok_b)} 個  → {primes[np.nonzero(~ok_b)[0]]} の次（隙間 {gaps[~ok_b]}）")

ax = fig.add_subplot(gs[1, 0])
vals, cnt = np.unique(gaps, return_counts=True)
colors = ["#d9541e" if isp[v] else "#9ec3ee" for v in vals]
ax.bar(vals, cnt/len(gaps), color=colors, width=1.6)
ax.set_xlim(0, 120); ax.set_xlabel("隙間 p_{k+1} − p_k"); ax.set_ylabel("割合")
ax.set_title(f"妄想2(a)：隙間が素数なのは {frac_a*100:.0f}%（橙）。残り {100-frac_a*100:.0f}% は合成数（青）→ 一つの素数では届かない（最初は 7→11）", fontsize=9.5, loc="left")
ax.grid(alpha=0.25, axis="y")

ax = fig.add_subplot(gs[1, 1])
ax.scatter(gaps[::37], r_gap[::37], s=3, color="#1f5fbf", alpha=0.5, label="隙間 h を二つの素数の和に書く方法の数 r(h)（ゴールドバッハ）")
hh = np.arange(4, H+1, 2)
C2 = 0.6601618
def Sfun(h):
    v = 2*C2
    for p in small:
        if p > 2 and h % p == 0: v *= (p-1)/(p-2)
    return v
S_h = np.array([Sfun(int(h)) for h in hh])
ax.plot(hh, S_h*hh/(np.log(hh)**2)/2, color="#d9541e", lw=1.2, label="ハーディ＝リトルウッドの予言 S(h)·h/(2 log²h)  ＝ 昨日のリズム S(h) と同じ")
ax.set_xlim(0, 160); ax.set_ylim(0, 25)
ax.set_xlabel("隙間 h"); ax.set_ylabel("r(h)")
ax.set_title("妄想2(b)：素数一つ、または素数二つの和なら必ず届く（2→3 以外に反例なし）。方法の数は昨日のリズム S(h) で決まる", fontsize=9.5, loc="left")
ax.legend(fontsize=8); ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[2, :])
L = np.cumsum(np.where(mu == 0, 0, 1)*0)  # placeholder
# リウヴィル λ(n) = (−1)^{素因数の個数（重複込み）}：ポリアの予想 L(x) ≤ 0 は 906,150,257 で破れる
Omega = np.zeros(N+1, np.int8)
for p in primes:
    pk = p
    while pk <= N:
        Omega[pk::pk] += 1; pk *= p
lam_liou = np.where(Omega % 2 == 0, 1, -1).astype(np.int64); lam_liou[0] = 0
Lx = np.cumsum(lam_liou)
ax.plot(xs, Lx[1000:]/np.sqrt(xs), color="#2a9d4f", lw=0.5, label="L(x)/√x ：リウヴィル（素因数の個数の偶奇）の均衡")
ax.axhline(0, color="k", lw=0.8)
ax.set_xscale("log"); ax.set_ylim(-2, 0.5)
ax.set_xlabel("x"); ax.set_ylabel("L(x)/√x")
ax.set_title("注意：均衡は「常に一方に傾く」のではない。ポリアは L(x) ≤ 0 と予想したが 906,150,257 で破れた（ハーゼルグローブ 1958）。RH が言うのは「傾きの幅が √x」だけ", fontsize=9.5, loc="left")
ax.legend(fontsize=8.5, loc="lower left"); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/center_and_dream.png")
