"""
ホッジ予想が止まっている正確な場所：ヴェイユ類（Weil 1977）
  A = C⁴/Λ（四次元複素トーラス ＝ 八次元実トーラス）、K = Q(i) が作用（S = i 倍, S² = −1）
  複素構造 J（J² = −1, JS = SJ）。ヴェイユ型 ＝ S が H^{1,0} 上に固有値 (i,i,−i,−i) を持つ（符号 (2,2)）。
  偏極 θ = Im H（H = diag(1,1,−1,−1) の K-エルミート形式）、J ∈ U(H) として J = P diag(±i,±i,∓i,∓i) P^{-1}, P ∈ U(2,2)。
  ヴェイユ類 W = ∧⁴_K V*：K-行列式 det_K(v1..v4) = D_R + i D_I の実部・虚部（整数係数の 4-形式、2 本）。
  確かめること
   (1) W は符号 (2,2) のときだけ (2,2) 型（ホッジ類＝どちらから見ても同じ影）。(3,1) 型では (3,1)+(1,3) に落ちる。
   (2) 一般の P：有理 (1,1) 類 NS_Q は θ だけ（LLL で探索）→ 布の積 D⁴ = Q·θ²（1 次元）。
       有理ホッジ類 B⁴_Q は θ², D_R, D_I の 3 次元（LLL）。→ W の 2 次元は「影はあるが布が見つからない」。
   (3) 特別な P（Q(i) 係数 → A ～ E_i⁴）：NS_Q は 16 次元、W ⊂ Sym² NS（布の積で書ける＝代数的）。
   (4) 光の定規（ホッジ＝リーマン）は W 上で正定値：光は W を布と同じ符号で見ている。
  現状：判別式 1 のヴェイユ型は Markman (2023), Floccari–Fu (2025) が代数的と証明。他の判別式は未解決。
"""
import itertools, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.linalg import expm, null_space
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]
rng = np.random.default_rng(3)

n = 8
# 実座標 (x1..x4, y1..y4), z_k = x_k + i y_k。S = i 倍
S = np.block([[np.zeros((4, 4)), -np.eye(4)], [np.eye(4), np.zeros((4, 4))]])
Hk = np.diag([1., 1., -1., -1.])                                   # K-エルミート形式、符号 (2,2)
def real_of(M):                                                     # C⁴ 上の C-線型写像 M を R⁸ の行列に
    return np.block([[M.real, -M.imag], [M.imag, M.real]])
def E_form(Hm):                                                     # θ(v,w) = Im H(v,w) の実 8×8 反対称行列
    A = Hm.real; B = Hm.imag
    return np.block([[B, A], [-A, B]]) if False else np.block([[Hm.imag, Hm.real], [-Hm.real, Hm.imag]])
theta2 = E_form(Hk.astype(complex))                                  # 2-形式（反対称行列）
def make_J(P, sig=(2, 2)):
    D = np.diag([1j]*sig[0] + [-1j]*sig[1])
    J = real_of(P @ D @ np.linalg.inv(P))
    # 偏極の正値性 θ(v, Jv) > 0 になるように向きを選ぶ
    T = theta2 @ J
    if np.linalg.eigvalsh((T + T.T)/2).min() < 0: J = -J
    return J
def random_U22():
    A = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4)); A = (A + A.conj().T)/2
    return expm(1j*np.linalg.inv(Hk) @ A)
def gaussian_U22():                                                 # Q(i) 係数の U(2,2) 元（ケイリー変換）
    while True:
        A = rng.integers(-2, 3, size=(4, 4)) + 1j*rng.integers(-2, 3, size=(4, 4)); A = (A + A.conj().T)/2
        X = 1j*np.linalg.inv(Hk) @ A                                 # X† H + H X = 0
        try: P = (np.eye(4) - X) @ np.linalg.inv(np.eye(4) + X)      # ∈ U(H)
        except np.linalg.LinAlgError: continue
        if np.linalg.cond(P) < 50: return P

# 外積代数の道具
def subsets(k): return list(itertools.combinations(range(n), k))
S2, S4 = subsets(2), subsets(4)
i2 = {s: a for a, s in enumerate(S2)}; i4 = {s: a for a, s in enumerate(S4)}
def ext_power(A, k):
    """A: V*→V*（n×n）の ∧^k V* への持ち上げ（基底 e_I, I = k-部分集合）"""
    Sk = subsets(k); M = np.zeros((len(Sk), len(Sk)), dtype=A.dtype)
    for b, Jb in enumerate(Sk):
        for a, Ia in enumerate(Sk):
            M[a, b] = np.linalg.det(A[np.ix_(Ia, Jb)])
    return M
def form2_to_vec(T):                                                # 反対称行列 → ∧² の座標
    return np.array([T[i, j] for (i, j) in S2])
def wedge22(u, v):                                                  # ∧²×∧² → ∧⁴
    out = np.zeros(len(S4))
    for (a, Ia) in enumerate(S2):
        for (b, Jb) in enumerate(S2):
            if set(Ia) & set(Jb): continue
            K = tuple(sorted(Ia + Jb)); perm = list(Ia + Jb)
            sgn = np.linalg.det(np.eye(4)[np.argsort(perm)])
            out[i4[K]] += sgn*u[a]*v[b]
    return out
Q44 = np.zeros((70, 70))                                            # ∧⁴×∧⁴ → ∧⁸ の対（定規）
for a, Ia in enumerate(S4):
    for b, Jb in enumerate(S4):
        if set(Ia) & set(Jb): continue
        perm = list(Ia + Jb); Q44[a, b] = np.linalg.det(np.eye(8)[perm])

# ヴェイユ類：det_K の実部・虚部
Zc = np.hstack([np.eye(4), 1j*np.eye(4)])                           # e_a の K-座標 z(e_a) ∈ C⁴
detK = np.array([np.linalg.det(Zc[:, list(I)]) for I in S4])
W_R, W_I = detK.real.round().astype(int), detK.imag.round().astype(int)
print("ヴェイユ類 W_R, W_I：整数係数の 4-形式、非零成分", np.count_nonzero(W_R), np.count_nonzero(W_I))

def hodge_data(J):
    Jt = J.T                                                         # V* 上の作用
    C2 = ext_power(Jt, 2); C4 = ext_power(Jt, 4)                     # ヴェイユ作用素 i^{p−q}
    # H^{1,0} = {φ : J^T φ = iφ}
    w, v = np.linalg.eig(Jt); H10 = v[:, np.abs(w - 1j) < 1e-9]
    Omega = np.array([np.linalg.det(H10[list(I), :]) for I in S4])   # (4,0) のプリュッカー座標
    # (2,2)_R = {u : C4 u = u, Q44(u, Re Ω) = Q44(u, Im Ω) = 0}
    cond22 = np.vstack([C4 - np.eye(70), (Q44 @ Omega.real)[None, :], (Q44 @ Omega.imag)[None, :]])
    H22 = null_space(cond22)
    H11 = null_space(C2 - np.eye(28))
    return C2, C4, Omega, cond22, H22, H11

def lll_kernel(cond, scale=10**8):
    """整数ベクトル v（成分小）で cond v ≈ 0 を満たすものを LLL で探す。見つかった短いベクトルの張る階数と基底を返す"""
    m, d = cond.shape
    B = np.hstack([np.eye(d, dtype=object), np.rint(cond.T*scale).astype(object)])
    Mtx = DomainMatrix([[ZZ(int(x)) for x in row] for row in B.tolist()], (d, d+m), ZZ)
    red = np.array(Mtx.lll().to_Matrix().tolist(), dtype=float)
    vs = [r[:d] for r in red if 0 < np.abs(r[:d]).max() <= 100 and np.abs(cond @ r[:d]).max() < 1e-6]
    return np.array(vs).T if vs else np.zeros((d, 0))

def analyze(J, label):
    C2, C4, Omega, cond22, H22, H11 = hodge_data(J)
    th = form2_to_vec(theta2)
    out = {}
    out["dim22"] = H22.shape[1]
    # (1) W の型
    for name, w in [("W_R", W_R.astype(float)), ("W_I", W_I.astype(float))]:
        proj = H22 @ (H22.T @ w)                                     # (2,2) への直交射影（基底は正規直交）
        out[name + "_in22"] = np.linalg.norm(proj)/np.linalg.norm(w)
    # (2) 有理 (1,1) 類 NS_Q と布の積 D⁴
    NS = lll_kernel(C2 - np.eye(28))
    out["rankNS"] = np.linalg.matrix_rank(NS) if NS.size else 0
    out["theta_in_NS"] = np.abs((C2 - np.eye(28)) @ th).max() < 1e-9
    prods = [wedge22(NS[:, a], NS[:, b]) for a in range(NS.shape[1]) for b in range(a, NS.shape[1])] if NS.size else []
    D4 = np.array(prods).T if prods else wedge22(th, th)[:, None]
    out["rankD4"] = np.linalg.matrix_rank(D4)
    # 有理ホッジ類 B⁴_Q
    B4 = lll_kernel(cond22)
    out["rankB4"] = np.linalg.matrix_rank(B4) if B4.size else 0
    # W ⊂ D⁴ ?
    res = [np.linalg.norm(D4 @ np.linalg.lstsq(D4, w, rcond=None)[0] - w)/np.linalg.norm(w) for w in [W_R.astype(float), W_I.astype(float)]]
    out["W_in_D4_resid"] = max(res)
    # (4) 光の定規 Q44 を W 上で
    Wm = np.stack([W_R, W_I], 1).astype(float)
    G = Wm.T @ Q44 @ Wm
    out["HR_eigs"] = np.linalg.eigvalsh(G)
    out["theta_wedge_W"] = max(np.abs(wedge22(th, th) @ Q44 @ W_R.astype(float)), np.abs(wedge22(th, th) @ Q44 @ W_I.astype(float)))
    print(f"[{label}]  dim (2,2)_R = {out['dim22']},  W の (2,2) 成分の割合 = {out['W_R_in22']:.3f}, {out['W_I_in22']:.3f}")
    print(f"      有理 (1,1) 類 NS_Q の階数 = {out['rankNS']}（θ ∈ NS: {out['theta_in_NS']}）, 布の積 D⁴ の階数 = {out['rankD4']}, 有理ホッジ類 B⁴_Q の階数 = {out['rankB4']}")
    print(f"      W は布の積で書けるか：残差 = {out['W_in_D4_resid']:.2e}  → {'書ける（代数的）' if out['W_in_D4_resid'] < 1e-6 else '書けない'}")
    print(f"      光の定規 ∫W∧W の固有値 = {np.round(out['HR_eigs'], 3)}（正定値なら光は W を布と同じ符号で見ている）, θ²·W = {out['theta_wedge_W']:.1e}")
    return out

print("=== (1) 符号の違い：ヴェイユ型 (2,2) と非ヴェイユ型 (3,1) ===")
P = random_U22()
resA = analyze(make_J(P, (2, 2)), "一般のヴェイユ型 (2,2)")
resB = analyze(make_J(P, (3, 1)), "同じ P で符号 (3,1)")
print("=== (3) 特別な点：P が Q(i) 係数（A ～ E_i⁴）===")
resC = analyze(make_J(np.eye(4, dtype=complex), (2, 2)), "特別なヴェイユ型（A = E_i⁴）")

# ---------- 図 ----------
fig = plt.figure(figsize=(15, 9), dpi=150, constrained_layout=True)
gs = fig.add_gridspec(2, 3)
# H⁴ のホッジ分解と有理ホッジ類
ax = fig.add_subplot(gs[0, :2])
types = [(4, 0), (3, 1), (2, 2), (1, 3), (0, 4)]; dims = [1, 16, 36, 16, 1]
ax.bar(range(5), dims, color=["#ccc", "#ccc", "#9ec3ee", "#ccc", "#ccc"], width=0.6)
ax.bar([2], [1], color="#1f5fbf", width=0.6, label="布の積 θ²（代数的、1 次元）")
ax.bar([2], [2], bottom=[1], color="#d9541e", width=0.6, label="ヴェイユ類 W（ホッジ類だが布が見つからない、2 次元）")
ax.text(2, 4.2, "有理ホッジ類 B⁴_Q：3 次元\n（36 次元の影のうち、整数の目に見えるのは 3 本だけ）", ha="center", fontsize=9)
ax.set_xticks(range(5)); ax.set_xticklabels([f"H^{{{p},{q}}}\n{d} 次元" for (p, q), d in zip(types, dims)])
ax.set_ylabel("次元"); ax.set_ylim(0, 40)
ax.set_title("一般のヴェイユ型アーベル四次元多様体の H⁴（70 次元）：影の分解と、整数の点", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="upper right"); ax.grid(alpha=0.25, axis="y")

ax = fig.add_subplot(gs[0, 2])
labs = ["ヴェイユ型\n(2,2)", "非ヴェイユ型\n(3,1)"]
vals = [resA["W_R_in22"], resB["W_R_in22"]]
ax.bar(labs, vals, color=["#d9541e", "#999"])
ax.set_ylim(0, 1.15); ax.set_ylabel("W の (2,2) 成分の割合")
for i_, v in enumerate(vals): ax.text(i_, v + 0.03, f"{v:.3f}", ha="center")
ax.set_title("(1) W がホッジ類（どちらから見ても同じ影）\nになるのは符号 (2,2) のときだけ", fontsize=10, loc="left"); ax.grid(alpha=0.25, axis="y")

ax = fig.add_subplot(gs[1, :2])
cats = ["有理 (1,1) 類 NS_Q", "布の積 D⁴", "有理ホッジ類 B⁴_Q"]
gen = [resA["rankNS"], resA["rankD4"], resA["rankB4"]]; spc = [resC["rankNS"], resC["rankD4"], resC["rankB4"]]
x = np.arange(3)
ax.bar(x - 0.2, gen, 0.38, color="#d9541e", label="一般のヴェイユ型（End = Q(i)）")
ax.bar(x + 0.2, spc, 0.38, color="#1f5fbf", label="特別な点（E_i⁴ と同種）")
for i_, (g_, s_) in enumerate(zip(gen, spc)):
    ax.text(i_ - 0.2, g_ + 0.4, str(g_), ha="center", fontsize=9); ax.text(i_ + 0.2, s_ + 0.4, str(s_), ha="center", fontsize=9)
ax.set_xticks(x); ax.set_xticklabels(cats); ax.set_ylabel("階数（LLL で見つけた整数ベクトルの張る次元）")
ax.set_title(f"(2)(3) 一般の点では B⁴_Q = θ² ⊕ W で W は布の積に入らない（残差 {resA['W_in_D4_resid']:.2f}）。"
             f"特別な点では布が増え、W は布の積で書ける（残差 {resC['W_in_D4_resid']:.0e}）", fontsize=10, loc="left")
ax.legend(fontsize=8.5); ax.grid(alpha=0.25, axis="y")

ax = fig.add_subplot(gs[1, 2])
ev = resA["HR_eigs"]
ax.bar(["固有値 1", "固有値 2"], ev, color="#c99a2e")
ax.axhline(0, color="k", lw=0.8); ax.set_ylabel("∫ W∧W の固有値")
ax.set_title("(4) 光の定規は W 上で正定値：\n光は W を「布と同じ符号」で見ている（θ²·W = 0、原始的）", fontsize=10, loc="left"); ax.grid(alpha=0.25, axis="y")
fig.suptitle("ホッジ予想が止まっている場所：ヴェイユ類 ― 影はあり、光の定規も布と同じ符号で見ているのに、布が見つからない（判別式 1 は 2023 年に解決、他は未解決）", fontsize=11.5)
fig.savefig("figures/weil_classes.png")
