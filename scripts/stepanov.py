"""
ステパノフの方法（1969）：曲面を作らずに、点を数える ― 兄弟の世界のリーマン予想の「もう一つの証明」
  曲線 y² = f(x) （f = x³+x+1）, 体 F_p。 S = {x ∈ F_p : f(x) が 0 でない平方}, N = |S|。
  本当の値：N = (p − #根 + Σ_x χ(f(x)))/2、ハッセ＝ヴェイユ |Σχ| ≤ 2√p より N ≤ p/2 + √p。
  目標：N ≤ p/2 + C√p を、補助多項式の次数だけから示す（曲面 E×E も交点数も使わない）。
  補助多項式  R(x) = Σ_{j<J} x^{pj} [ r_j(x) + s_j(x) f(x)^{(p−1)/2} ],   deg r_j, s_j < m
  仕掛け：
   (i)  S 上では f^{(p−1)/2} = 1（オイラーの判定法）、x^{pj} = x^j（フロベニウス：p 歩で自分に戻る）
   (ii) 標数 p では x^{pj} のハッセ微分は 0（k < p）→ D^k R = Σ x^{pj} D^k[r_j + s_j f^e]
   (iii) D^b(f^e) = f^{e−b} c_b（c_b は次数 ≤ 2b の多項式）なので、S 上で
           f^k D^k R = U_k(x) := Σ_j x^j [ f^k D^k r_j + Σ_{a+b=k} D^a s_j · f^a c_b ]   （次数 ≤ 2k+J+m−2）
        「U_k ≡ 0 (k < M)」は S を知らずに書ける線型条件：M(J+m−1) + M(M−1) 本。未知数は 2Jm 個。
        未知数 > 条件 なら非零解があり、R は S のすべての点で重複度 ≥ M で消える。
   (iv) R ≢ 0 なら  N·M ≤ deg R ≤ p(J−1) + m − 1 + 3(p−1)/2。
        J ≈ M/2 + 1, M ≈ c√p と取ると  N ≤ p/2 + C√p。
  ここでは F_p 上で実際に線型方程式を解いて R を作り、(a) R ≢ 0、(b) S の全点で重複度 ≥ M、(c) N·M ≤ deg R を検算する。
"""
import numpy as np, math, matplotlib, time
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

def pmul(a, b, p): return np.convolve(a, b) % p
def ppow(a, e, p):
    r = np.array([1], dtype=np.int64)
    for _ in range(e): r = pmul(r, a, p)
    return r
def pdivmod(a, b, p):
    a = a.copy() % p; db = len(b) - 1; inv = pow(int(b[-1]), p-2, p)
    q = np.zeros(max(len(a) - db, 1), dtype=np.int64)
    for i in range(len(a) - 1, db - 1, -1):
        c = (a[i]*inv) % p
        if c: q[i - db] = c; a[i-db:i+1] = (a[i-db:i+1] - c*b) % p
    return q, a[:db] % p
def hasse(a, k, p):
    c = np.array([math.comb(n, k) % p for n in range(len(a))], dtype=np.int64)
    out = (a*c) % p
    return out[k:] if len(out) > k else np.zeros(1, dtype=np.int64)
def trim(a):
    nz = np.nonzero(a)[0]; return a[:nz[-1]+1] if len(nz) else np.zeros(1, dtype=np.int64)
def monomial(d, c): return np.r_[np.zeros(d, dtype=np.int64), [c]]

def kernel_mod_p(A, p):
    """F_p 上のガウス消去で、核の非零ベクトルを一つ返す"""
    A = A.copy() % p; rows, cols = A.shape; piv = []; r = 0
    for c in range(cols):
        if r >= rows: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: A[[r, i]] = A[[i, r]]
        A[r, c:] = (A[r, c:]*pow(int(A[r, c]), p-2, p)) % p
        others = np.nonzero(A[:, c])[0]; others = others[others != r]
        if len(others): A[others, c:] = (A[others, c:] - np.outer(A[others, c], A[r, c:])) % p
        piv.append(c); r += 1
    pivset = set(piv); free = [c for c in range(cols) if c not in pivset]
    if not free: return None, len(piv)
    fc = free[0]; v = np.zeros(cols, dtype=np.int64); v[fc] = 1
    for i, pc in enumerate(piv): v[pc] = (-A[i, fc]) % p
    return v, len(piv)

def choose_params(p):
    """上界 deg R / M を最小にする (M, J, m)。条件：未知数 > 条件、2m+1 < p（R ≢ 0 の補題の仮定）"""
    best = None; e = (p-1)//2
    for M in range(3, int(1.5*np.sqrt(p)) + 2):
        for J in range(M//2 + 1, M + 3):
            den = 2*J - M
            if den <= 0: continue
            m = (M*(J-1) + M*(M-1))//den + 1                    # 2Jm > M(J+m−1) + M(M−1)
            if 2*m + 1 >= p: continue
            bound = (p*(J-1) + m - 1 + 3*e)/M
            if best is None or bound < best[0]: best = (bound, M, J, m)
    return best

def stepanov(p, f=(1, 1, 0, 1)):
    f = np.array(f, dtype=np.int64); e = (p-1)//2
    param_bound, M, J, m = choose_params(p)
    xs = np.arange(p); fx = np.polyval(f[::-1], xs) % p
    leg = np.array([pow(int(v), e, p) for v in fx]); S = xs[leg == 1]; N = len(S)
    fe = ppow(f, e, p); fpow = [ppow(f, k, p) for k in range(M+1)]
    cb = []
    for b in range(M):
        q, rem = pdivmod(hasse(fe, b, p), ppow(f, e-b, p), p); assert not rem.any(); cb.append(trim(q))
    nU = [2*k + J + m - 1 for k in range(M)]; off = np.cumsum([0] + nU)
    ncol = 2*J*m; A = np.zeros((off[-1], ncol), dtype=np.int64)
    def put(k, col, poly):
        assert len(poly) <= nU[k] or not poly[nU[k]:].any()
        L = min(len(poly), nU[k]); A[off[k]:off[k]+L, col] = (A[off[k]:off[k]+L, col] + poly[:L]) % p
    col = 0
    for j in range(J):
        for n in range(m):                                       # r_j の x^n の係数
            for k in range(min(n, M-1) + 1):
                put(k, col, pmul(fpow[k], monomial(j+n-k, math.comb(n, k) % p), p))
            col += 1
        for n in range(m):                                       # s_j の x^n の係数
            for k in range(M):
                for a in range(min(n, k) + 1):
                    put(k, col, pmul(pmul(fpow[a], cb[k-a], p), monomial(j+n-a, math.comb(n, a) % p), p))
            col += 1
    t0 = time.time(); v, rank = kernel_mod_p(A, p); dt = time.time() - t0
    R = np.zeros(p*(J-1) + m + 3*e + 1, dtype=np.int64)
    for j in range(J):
        rj = v[(2*j)*m:(2*j+1)*m]; sj = v[(2*j+1)*m:(2*j+2)*m]
        g = pmul(sj, fe, p); g[:m] = (g[:m] + rj) % p
        R[p*j:p*j+len(g)] = (R[p*j:p*j+len(g)] + g) % p
    R = trim(R); degR = len(R) - 1
    # 検算：S の各点 x0 でのテイラー係数 D^k R(x0), k = 0..M（ホーナー法、S 全体を一度に）
    orders = np.full(N, M + 1)
    for k in range(M, -1, -1):
        Dk = hasse(R, k, p); val = np.zeros(N, dtype=np.int64)
        for c in Dk[::-1]: val = (val*S + c) % p
        orders[val != 0] = k
    return dict(p=p, M=M, J=J, m=m, N=N, nonzero=bool(R.any()), degR=degR, minorder=int(orders.min()), rows=A.shape[0], cols=ncol,
                rank=rank, time=dt, bound=degR/M, param_bound=param_bound, orders=orders)

results = []
for p in [101, 199, 307, 401]:
    r = stepanov(p); results.append(r)
    print(f"p={p}: M={r['M']}, J={r['J']}, m={r['m']}  未知数 {r['cols']} > 条件 {r['rows']}（階数 {r['rank']}） →  R ≢ 0: {r['nonzero']},  deg R = {r['degR']},  "
          f"S の全点（{r['N']} 個）での最小重複度 = {r['minorder']}（≥ M = {r['M']}）")
    print(f"      本当の N = {r['N']}   実際の R の上界 deg R / M = {r['bound']:.1f}   保証（パラメータ）= {r['param_bound']:.1f}   p/2 = {p/2:.1f}   ハッセ＝ヴェイユ p/2+√p = {p/2 + np.sqrt(p):.1f}   自明 p = {p}"
          f"   （N·M = {r['N']*r['M']} ≤ deg R = {r['degR']}: {r['N']*r['M'] <= r['degR']}）  [{r['time']:.0f}s]")

# ---------- 図 ----------
fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), dpi=150, constrained_layout=True)
ax = axes[0]
ps = [r["p"] for r in results]; x = np.arange(len(ps))
ax.bar(x - 0.3, [r["N"] for r in results], 0.28, color="#1f5fbf", label="本当の N（f(x) が平方になる x の個数）")
ax.bar(x, [r["param_bound"] for r in results], 0.28, color="#f3c1ae", label="ステパノフの保証（パラメータだけから）")
ax.bar(x, [r["bound"] for r in results], 0.28, color="#d9541e", label="実際に作った R の deg R / M")
ax.bar(x + 0.3, ps, 0.28, color="#ccc", label="自明な上界 p")
for i_, r in enumerate(results):
    ax.plot([i_-0.45, i_+0.45], [r["p"]/2 + np.sqrt(r["p"])]*2, color="#2a9d4f", lw=1.8)
ax.plot([], [], color="#2a9d4f", lw=1.8, label="p/2 + √p（ハッセ＝ヴェイユ）")
ax.set_xticks(x); ax.set_xticklabels([f"p = {p}" for p in ps]); ax.set_ylabel("個数")
ax.set_title("曲面を作らずに点を数える", fontsize=10.5, loc="left"); ax.legend(fontsize=8); ax.grid(alpha=0.25, axis="y")
ax = axes[1]
ax.plot(ps, [(r["param_bound"] - r["p"]/2)/np.sqrt(r["p"]) for r in results], "o--", color="#f3a58a", lw=2, label="(ステパノフの保証 − p/2)/√p")
ax.plot(ps, [(r["bound"] - r["p"]/2)/np.sqrt(r["p"]) for r in results], "o-", color="#d9541e", lw=2, label="(実際に作った R の上界 − p/2)/√p")
ax.plot(ps, [(r["N"] - r["p"]/2)/np.sqrt(r["p"]) for r in results], "s-", color="#1f5fbf", lw=2, label="(本当の N − p/2)/√p")
ax.axhline(1, color="#2a9d4f", lw=1.5, label="ハッセ＝ヴェイユの定数 1")
ax.plot(ps, [np.sqrt(p)/2 for p in ps], "--", color="#999", label="自明な上界なら √p/2（増え続ける）")
ax.set_xlabel("p"); ax.set_ylabel("√p で割った超過分"); ax.set_ylim(-1.5, 11)
ax.set_title("超過分は √p の定数倍に留まる", fontsize=10.5, loc="left"); ax.legend(fontsize=8); ax.grid(alpha=0.25)
ax = axes[2]
r = results[-1]
ax.hist(r["orders"], bins=np.arange(-0.5, r["M"] + 2.5), color="#8a5fbf", edgecolor="white")
ax.axvline(r["M"] - 0.5, color="k", ls="--", lw=1); ax.text(r["M"] - 0.7, ax.get_ylim()[1]*0.9, f"M = {r['M']}", ha="right")
ax.set_xlabel(f"R が消える重複度（{r['M']+1} は「{r['M']+1} 以上」）"); ax.set_ylabel("S の点の個数")
ax.set_title(f"p = {r['p']}：R は S の {r['N']} 点すべてで重複度 ≥ {r['M']}", fontsize=10.5, loc="left"); ax.grid(alpha=0.25, axis="y")
fig.suptitle("ステパノフの方法 ― すべての点で高い重複度で消える補助多項式 R を、フロベニウス（x^p = x）で安く作り、次数で点の個数を抑える", fontsize=11.5)
fig.savefig("figures/stepanov.png")
