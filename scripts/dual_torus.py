"""
向こう側から見ると巻き方が逆になる：双対トーラスへの射影（フーリエ＝向井変換）
  A = E_i × E_i = C²/Λ,  Λ の基底 (1,0),(i,0),(0,1),(0,i),  H^k(A) = ∧^k V*（V* = R^4, 基底 e1..e4）, 全体 16 次元
  上へ巻く L = ω∧·（ω = e12 + e34：光＝ケーラー類）,  下へ折り返す Λ = L の随伴,  H = (2 − k) on H^k（[Λ, L] = H, [H, Λ] = 2Λ）
  双対トーラス Â：H^k(Â) = ∧^k V。射影 Φ: ∧^k V* → ∧^{4−k} V,  e_J ↦ sgn(J, J^c) f_{J^c}（ポアンカレ束 e^{c₁(P)} で積分）
  主偏極 φ_ω: V → V* で Â ≅ A と同一視し Φ' = ∧(φ_ω)∘Φ : H^k(A) → H^{4−k}(A)
  確かめること
   (1) Φ' L Φ'^{-1} = ±Λ   …… 向こう側では「巻く」が「折り返す」になる（レフシェッツ型標準予想の、トーラスでの証明の筋）
   (2) Φ'(布) = 布         …… NS(A)（整数 (1,1) 類、4 本）が NS(A) に写る
   (3) ヴェイユ作用素 C = ∧²J*（(p,q) を i^{p−q} 倍：向こう側から見る）は布を固定し、(2,0)+(0,2) の影を反転する
"""
import itertools, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

n = 4
subsets = [frozenset(c) for k in range(n+1) for c in itertools.combinations(range(n), k)]   # 16 個、次数順
idx = {S: i for i, S in enumerate(subsets)}
deg = np.array([len(S) for S in subsets])
def sgn_merge(A, B):
    """e_A ∧ e_B = sgn · e_{A∪B}（A∩B=∅）"""
    a = sorted(A); b = sorted(B); s = 1
    for x in a:
        s *= (-1)**sum(1 for y in b if y < x)
    return s
def wedge_matrix(two_form):
    """two_form: dict {(i,j): coeff} (i<j)  → 左から ∧ する 16×16 行列"""
    M = np.zeros((16, 16))
    for S in subsets:
        for (i, j), c in two_form.items():
            if i in S or j in S: continue
            T = S | {i, j}
            M[idx[T], idx[S]] += c*sgn_merge({i, j}, S)
    return M
omega = {(0, 1): 1.0, (2, 3): 1.0}
L = wedge_matrix(omega); Lam = L.T; Hw = np.diag(2.0 - deg)
print("(0) sl2 の交換関係 [Λ, L] = H :", np.allclose(Lam@L - L@Lam, Hw), "   [H, Λ] = 2Λ :", np.allclose(Hw@Lam - Lam@Hw, 2*Lam))
# ハード・レフシェッツ：L^k : H^{2−k} → H^{2+k} が同型
for k in [1, 2]:
    src = deg == 2-k; dst = deg == 2+k
    Mk = np.linalg.matrix_power(L, k)[np.ix_(dst, src)]
    print(f"    ハード・レフシェッツ L^{k}: H^{2-k} → H^{2+k} の階数 {np.linalg.matrix_rank(Mk)} / {src.sum()}")

# 射影 Φ : ∧^k V* → ∧^{4−k} V,  e_J ↦ sgn(J, J^c) f_{J^c}
Phi = np.zeros((16, 16))
full = frozenset(range(n))
for S in subsets:
    Sc = full - S
    Phi[idx[Sc], idx[S]] = sgn_merge(S, Sc)
# 主偏極 φ_ω : V → V*,  f_i ↦ ω(f_i, ·)
Om = np.zeros((n, n))
for (i, j), c in omega.items(): Om[i, j] = c; Om[j, i] = -c
def ext_power_map(A):
    """線型写像 A: V→V* (n×n) の外積代数への持ち上げ ∧A : ∧^k V → ∧^k V*（基底の添字で同一視）"""
    M = np.zeros((16, 16))
    for S in subsets:
        k = len(S); s = sorted(S)
        for T in subsets:
            if len(T) != k: continue
            t = sorted(T)
            M[idx[T], idx[S]] = np.linalg.det(A[np.ix_(t, s)]) if k else 1.0
    return M
PhiP = ext_power_map(Om) @ Phi                       # Φ' : H^k(A) → H^{4−k}(A)
conj = PhiP @ L @ np.linalg.inv(PhiP)
c1 = np.sum(conj*Lam)/np.sum(Lam*Lam)
print(f"(1) 向こう側へ射影すると巻き方が逆になる： Φ' L Φ'^{{-1}} = ({c1:+.0f})·Λ  …… 一致 {np.allclose(conj, c1*Lam)}")
print(f"    Φ' は次数 k を 4−k に写す：", all(np.allclose(deg[np.nonzero(PhiP[:, j])[0]], 4-deg[j]) for j in range(16)))

# 複素構造と H^{1,1}, NS
Pi = np.array([[1, 1j, 0, 0], [0, 0, 1, 1j]])
pairs = [tuple(sorted(S)) for S in subsets if len(S) == 2]; i2 = [idx[frozenset(p)] for p in pairs]
om = np.array([Pi[0, i]*Pi[1, j] - Pi[0, j]*Pi[1, i] for (i, j) in pairs])       # H^{2,0} のプリュッカー座標
Q = np.zeros((6, 6))
for a, (i, j) in enumerate(pairs):
    for b, (k, l) in enumerate(pairs):
        if len({i, j, k, l}) == 4: Q[a, b] = np.linalg.det(np.eye(4)[[i, j, k, l]])
Qo = Q if np.real(om@Q@np.conj(om)) > 0 else -Q
Mcond = np.vstack([Qo@om.real, Qo@om.imag])
g = np.indices([7]*6).reshape(6, -1).T - 3; g = g[np.any(g != 0, axis=1)]
NS = g[np.abs(g@Mcond.T).max(axis=1) < 1e-9].T.astype(float)
rho = np.linalg.matrix_rank(NS)
# 整数基底（LLL の代わりに、見つかったベクトルから階段状に 4 本選ぶ）
basis = []
for v in NS.T[np.argsort(np.abs(NS.T).sum(axis=1))]:
    if np.linalg.matrix_rank(np.array(basis + [v])) > len(basis): basis.append(v)
    if len(basis) == rho: break
NSb = np.array(basis).T
names = {(0, 1): "e12", (0, 2): "e13", (0, 3): "e14", (1, 2): "e23", (1, 3): "e24", (2, 3): "e34"}
def show(v): return " ".join(f"{c:+.0f}{names[p]}" for c, p in zip(v, pairs) if abs(c) > 1e-9)
P2 = PhiP[np.ix_(i2, i2)]                            # Φ' の H² → H² 部分
imgs = P2 @ NSb
print(f"(2) 布 ρ = {rho} 本。Φ'(布) が布の格子に入るか：", all(np.abs(np.linalg.lstsq(NSb, imgs[:, j], rcond=None)[0] - np.round(np.linalg.lstsq(NSb, imgs[:, j], rcond=None)[0])).max() < 1e-9 for j in range(rho)),
      "  階数 [NS, Φ'NS] =", np.linalg.matrix_rank(np.hstack([NSb, imgs])))
for j in range(rho): print(f"      {show(NSb[:, j]):18s} →  {show(imgs[:, j])}")
# H^{1,1} と H^{2,0}+H^{0,2} の Φ' による像
_, _, Vt = np.linalg.svd(Mcond); H11 = Vt[2:].T; H20 = Vt[:2].T
print("    Φ'(H^{1,1}) ⊂ H^{1,1}:", np.abs(Mcond @ (P2@H11)).max() < 1e-9, "   Φ'(H^{2,0}+H^{0,2}) ⊂ H^{2,0}+H^{0,2}:", np.abs((H11.T @ (P2@H20))).max() < 1e-9,
      "   光 ω ↦", show(P2 @ np.array([1, 0, 0, 0, 0, 1.])))

# (3) ヴェイユ作用素 C = ∧² J*
Jr = np.zeros((4, 4))                               # 複素構造 J: V → V （z ↦ iz）
Rp = np.vstack([Pi.real, Pi.imag])                  # x ↦ (Re z, Im z)
Jc = np.array([[0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0]], float)   # (Re,Im) ↦ (−Im, Re)
Jr = np.linalg.inv(Rp) @ Jc @ Rp
C = ext_power_map(Jr.T)[np.ix_(i2, i2)]             # V* 上は転置、∧² へ
evC = np.linalg.eigvalsh((C + C.T)/2)
print("(3) ヴェイユ作用素 C（向こう側から見る）の H² 上の固有値:", np.round(evC, 6), "  C² = 1:", np.allclose(C@C, np.eye(6)))
print("    布は固定される C v = v:", np.allclose(C@NSb, NSb), "   (2,0)+(0,2) は反転 C v = −v:", np.allclose(C@H20, -H20))

# ---------- 図 ----------
fig = plt.figure(figsize=(15, 9), dpi=150, constrained_layout=True)
gs = fig.add_gridspec(2, 3)
def mat(ax, M, title):
    im = ax.imshow(M, cmap="RdBu_r", vmin=-1.5, vmax=1.5)
    for b in np.cumsum([1, 4, 6, 4])-0.5: ax.axhline(b, color="k", lw=0.5); ax.axvline(b, color="k", lw=0.5)
    ax.set_xticks([0, 2.5, 7.5, 12.5, 15]); ax.set_xticklabels(["H0", "H1", "H2", "H3", "H4"], fontsize=8)
    ax.set_yticks([0, 2.5, 7.5, 12.5, 15]); ax.set_yticklabels(["H0", "H1", "H2", "H3", "H4"], fontsize=8)
    ax.set_title(title, fontsize=10, loc="left")
mat(fig.add_subplot(gs[0, 0]), L, "上へ巻く L = ω∧·（光を一枚重ねる）")
mat(fig.add_subplot(gs[0, 1]), Lam, "下へ折り返す Λ（L の随伴）  [Λ, L] = H")
mat(fig.add_subplot(gs[0, 2]), conj, f"向こう側で見た L：Φ' L Φ'^(-1) = {c1:+.0f}·Λ\n（双対トーラスへ射影すると、巻く↔折り返す）")

ax = fig.add_subplot(gs[1, 0])
dims = [1, 4, 6, 4, 1]
for k, d in enumerate(dims):
    ax.scatter([k]*d, np.arange(d) - (d-1)/2, s=60, color="#1f5fbf", zorder=3)
    if k < 4: ax.annotate("", xy=(k+0.85, 0.6), xytext=(k+0.15, 0.6), arrowprops=dict(arrowstyle="->", color="#d9541e", lw=1.8))
    if k > 0: ax.annotate("", xy=(k-0.85, -0.6), xytext=(k-0.15, -0.6), arrowprops=dict(arrowstyle="->", color="#2a9d4f", lw=1.8))
ax.text(2, 3.3, "L：上へ巻く →", color="#d9541e", ha="center", fontsize=10); ax.text(2, -3.5, "← Λ：下へ折り返す", color="#2a9d4f", ha="center", fontsize=10)
ax.annotate("", xy=(3.5, 1.9), xytext=(0.5, 1.9), arrowprops=dict(arrowstyle="<->", color="#8a5fbf", lw=1.5, connectionstyle="arc3,rad=-0.3"))
ax.text(2, 2.6, "Φ'：向こう側へ射影（k ↔ 4−k）", color="#8a5fbf", ha="center", fontsize=9)
ax.set_xticks(range(5)); ax.set_xticklabels([f"H{k}\n{d} 次元" for k, d in enumerate(dims)]); ax.set_yticks([])
ax.set_xlim(-0.5, 4.5); ax.set_ylim(-4, 4); ax.set_title("包帯の段：ハード・レフシェッツ L^k: H^{2−k} ≅ H^{2+k}（定理）", fontsize=10, loc="left")
ax.text(2, 0, "ω", color="#c99a2e", fontsize=12, ha="center", va="center", fontweight="bold")

ax = fig.add_subplot(gs[1, 1]); ax.axis("off")
rows = [["布（NS, 整数の (1,1) 類）", "向こう側 Φ'(布)"]] + [[show(NSb[:, j]), show(imgs[:, j])] for j in range(rho)]
rows += [["光 ω = e12+e34", show(P2 @ np.array([1, 0, 0, 0, 0, 1.]))]]
t = ax.table(cellText=rows[1:], colLabels=rows[0], loc="center", cellLoc="center"); t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.6)
ax.set_title(f"(2) 布は布に写る：E_i×E_i の布 {rho} 本とその像（すべて整数係数）", fontsize=10, loc="left")

ax = fig.add_subplot(gs[1, 2])
labels = ["布 1", "布 2", "布 3", "布 4", "影 (2,0)+(0,2) ①", "②"]
vals = [1, 1, 1, 1, -1, -1]
ax.bar(range(6), vals, color=["#1f5fbf"]*4 + ["#d9541e"]*2)
ax.set_xticks(range(6)); ax.set_xticklabels(labels, fontsize=8, rotation=20); ax.axhline(0, color="k", lw=0.8)
ax.set_ylabel("C の固有値"); ax.set_ylim(-1.5, 1.5)
ax.set_title("(3) 向こう側から見る（ヴェイユ作用素 C = ∧²J*）：\n布は同じに見え（+1）、(2,0)+(0,2) の影は巻き方が逆（−1）", fontsize=10, loc="left")
fig.suptitle("向こう側から見ると巻き方が逆になる ― 双対トーラスへの射影は L を Λ に、布を布に写す（トーラスで「布の包帯」が締まる理由）", fontsize=12)
fig.savefig("figures/dual_torus.png")
