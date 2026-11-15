"""
複素構造の球面（ツイスター球）を回す：どの側から見ても布である影は何か
  T⁴ = R⁴/Z⁴ ≅ H/Z⁴（四元数）。左乗法 I, J, K が三つの複素構造、単位ベクトル (a,b,c) ごとに J_abc = aI+bJ+cK。
  北極 I = A（E_i×E_i）、南極 −I = Ā（複素共役）、赤道 = J, K, …
  光（ケーラー形式）ω_abc = a ω_I + b ω_J + c ω_K、 ω_J, ω_K = Re Ω, Im Ω（符号まで）
  各点で：布 = 整数の (1,1) 類 = {v ∈ Z⁶ : ∧²J_abc v = v}、本数 ρ(a,b,c)
  予言：ρ は光の方向 (a,b,c) が有理なとき 4、一般（無理）の点で 3。
        どの点でも布であり続ける整数類は、反自己双対な 3 本（H−V, Δ−H−V, Γ'：表裏の差で作った空間的な布）だけ。
        時間的な布（光 ω の方向）は、その光が選ばれている点でだけ布。
  → マークマンの道具「超正則層（どの複素構造でも正則）」が作る類は、この「どの側から見ても布」の類。
"""
import itertools, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

# 四元数 q = x1 + y1 i + x2 j + y2 k の左乗法行列（基底 1, i, j, k）
def lmul(q):
    a, b, c, d = q
    return np.array([[a, -b, -c, -d], [b, a, -d, c], [c, d, a, -b], [d, -c, b, a]], float)
I, J, K = lmul((0, 1, 0, 0)), lmul((0, 0, 1, 0)), lmul((0, 0, 0, 1))
assert np.allclose(I @ J, K) and np.allclose(I @ I, -np.eye(4))
pairs = list(itertools.combinations(range(4), 2))
names = ["dx1∧dy1", "dx1∧dx2", "dx1∧dy2", "dy1∧dx2", "dy1∧dy2", "dx2∧dy2"]
def ext2(A): return np.array([[np.linalg.det(A[np.ix_(list(Ia), list(Kb))]) for Kb in pairs] for Ia in pairs])
def kahler(Jm):                                    # ω(v,w) = g(Jv, w) の 2-形式座標
    return np.array([Jm[j, i] - Jm[i, j] for (i, j) in pairs])/2 * 2  # ω_ij = g(J e_i, e_j) = J_ji ; 反対称化
def kahler_vec(Jm): return np.array([Jm[j, i] for (i, j) in pairs])
wI, wJ, wK = kahler_vec(I), kahler_vec(J), kahler_vec(K)
print("三つの光：ω_I =", wI, " ω_J =", wJ, " ω_K =", wK)

grid = [np.array(c) for c in itertools.product(range(-2, 3), repeat=6) if any(c)]
G = np.array(grid, float)
def cloths_at(Jm):
    C = ext2(Jm.T)
    return np.abs(G @ C.T - G).max(axis=1) < 1e-9
def rho(mask): return np.linalg.matrix_rank(G[mask]) if mask.any() else 0

# 球面を走査
lat = np.linspace(-np.pi/2, np.pi/2, 61); lon = np.linspace(-np.pi, np.pi, 121)
R = np.zeros((len(lat), len(lon)), int)
always = np.ones(len(G), bool)
for i_, la in enumerate(lat):
    for j_, lo in enumerate(lon):
        a, b, c = np.sin(la), np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo)
        m = cloths_at(a*I + b*J + c*K)
        R[i_, j_] = rho(m)
        if (i_*7 + j_*3) % 5 == 0: always &= m           # 一部の点で「常に布」を絞る（全点だと重いので間引き）
# 特別な点
for nm, Jm in [("北極 I（= A）", I), ("南極 −I（= Ā）", -I), ("赤道 J", J), ("赤道 K", K), ("有理な方向 (I+2J+3K)/√14", (I + 2*J + 3*K)/np.sqrt(14)), ("一般の点 (I+√2J+√3K)/√6", (I + np.sqrt(2)*J + np.sqrt(3)*K)/np.sqrt(6))]:
    m = cloths_at(Jm); print(f"{nm:24s} 布の本数 ρ = {rho(m)}")
always_cl = G[always]
rk = np.linalg.matrix_rank(always_cl)
print(f"\n球面上のどの点でも布である整数類（|係数|≤2）：{always.sum()} 本、階数 {rk}")
# 基底を表示
basis = []
for v in always_cl[np.argsort(np.abs(always_cl).sum(axis=1))]:
    if np.linalg.matrix_rank(np.array(basis + [v])) > len(basis): basis.append(v)
    if len(basis) == rk: break
def show(v): return " ".join(f"{c:+.0f}{n}" for c, n in zip(v, names) if abs(c) > 1e-9)
for v in basis: print("    ", show(v), "   光 ω_I との内積 =", f"{float(v @ np.array([1,0,0,0,0,1.]))/1:.0f}" if False else f"{int(round(v[0]+v[5]))}")
# 反自己双対との一致
Q = np.zeros((6, 6))
for a_, (i, j) in enumerate(pairs):
    for b_, (k, l) in enumerate(pairs):
        if len({i, j, k, l}) == 4: Q[a_, b_] = np.linalg.det(np.eye(4)[[i, j, k, l]])
print("    すべて反自己双対（⋆v = −v）か：", all(np.allclose(Q @ v, -v) for v in always_cl))
print("    H−V, Δ−H−V, Γ' の張る空間と一致するか：", np.linalg.matrix_rank(np.vstack([always_cl, [[1,0,0,0,0,-1],[0,0,1,-1,0,0],[0,1,0,0,1,0]]])) == rk)
# 時間的な布 H = dx1∧dy1 が布である点の集合
Hc = np.array([1,0,0,0,0,0.]); Hmask = np.zeros_like(R, bool)
for i_, la in enumerate(lat):
    for j_, lo in enumerate(lon):
        a, b, c = np.sin(la), np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo)
        C = ext2((a*I + b*J + c*K).T); Hmask[i_, j_] = np.allclose(C @ Hc, Hc)
print("横の繊維 H が布である球面上の点：", Hmask.sum(), "/", Hmask.size, "（北極・南極だけ）" if Hmask.sum() <= 4*len(lon) else "")

# ---------- 図 ----------
fig = plt.figure(figsize=(15, 8.5), dpi=150, constrained_layout=True)
gs = fig.add_gridspec(2, 3, width_ratios=[2, 1, 1])
ax = fig.add_subplot(gs[:, 0])
im = ax.pcolormesh(np.degrees(lon), np.degrees(lat), R, cmap="Blues", vmin=2, vmax=4, shading="auto")
ax.scatter([0], [90], s=180, marker="*", color="#d9541e", edgecolor="k", zorder=5, label="北極 I = A（E_i×E_i）：光 ω")
ax.scatter([0], [-90], s=180, marker="*", color="#8a5fbf", edgecolor="k", zorder=5, label="南極 −I = Ā（共役）")
ax.scatter([0, 90, 180, -90, -180], [0]*5, s=90, marker="o", color="#c99a2e", edgecolor="k", zorder=5, label="赤道 ±J, ±K：光 Re Ω, Im Ω（共役のグラフが布になる側）")
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, ticks=[2, 3, 4]); cb.set_label("その複素構造での布の本数 ρ")
ax.set_xlabel("経度（赤道面 J–K の角度）"); ax.set_ylabel("緯度（I 成分）")
ax.set_title("複素構造の球面を回す：どの側から見ても布の本数は 3 以上。特別な点で 4", fontsize=10.5, loc="left")
ax.legend(fontsize=8, loc="lower left"); ax.grid(alpha=0.2)

ax = fig.add_subplot(gs[0, 1:]); ax.axis("off")
rows = [[show(v), "どの側から見ても布", "0"] for v in basis]
rows += [["dx1∧dy1（横の繊維 H）", "北極・南極でだけ布", "1"], ["共役のグラフ {(z,z̄)}", "赤道 K でだけ布", "0"]]
tb = ax.table(cellText=rows, colLabels=["整数類", "球面上で布である場所", "光 ω_I"], loc="center", cellLoc="center"); tb.auto_set_font_size(False); tb.set_fontsize(8.5); tb.scale(1, 1.6)
ax.set_title("どの側から見ても布：表裏の差で作った三本（反自己双対）。時間的な布は一点でだけ布", fontsize=10, loc="left")

ax = fig.add_subplot(gs[1, 1:])
ax.imshow(Hmask, cmap="Oranges", aspect="auto", origin="lower", extent=[-180, 180, -90, 90])
ax.set_xlabel("経度"); ax.set_ylabel("緯度"); ax.set_title("横の繊維 H が布である場所（橙）：北極と南極だけ\n光が回ると、時間的な布は布でなくなる", fontsize=10, loc="left")
fig.suptitle("A と Ā をまたぐ布 ― 複素構造の球面を回すと、光は回り、時間的な布は消え、表裏の差の三本だけがどの側から見ても布のまま（超正則＝マークマンの道具の縮図）", fontsize=11)
fig.savefig("figures/twistor_sphere.png")
