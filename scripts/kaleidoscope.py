"""
万華鏡：ピカール＝レフシェッツ反射で布の像を集める（E_i × E_i, H² = Z⁶, 符号 (3,3)）
  鏡 s_δ(x) = x + Q(x,δ) δ  （δ² = −2 の格子ベクトル δ。整数格子を保ち、Q を保つ）
  (A) 家族の万華鏡 G：偏極 θ = e12+e34 に直交する根 δ ∈ θ^⊥ による反射の群（主偏極アーベル曲面の族のモノドロミー Sp(4,Z) の像）
      → 一本の布 e12（横の繊維）から始めて像を集める。像は「どこかのトーラスの布」。このトーラスの布（NS に入る）は何本か。
  (B) このトーラスの万華鏡 W(NS)：布の格子 NS（階数 4, 符号 (1,3)）の中の根による反射の群
      → 一本の布から、他の布（縦の繊維、対角線、i のグラフ）に届くか。
  (C) 部分族の万華鏡の不変類：E×E'（Sp(2,Z)×Sp(2,Z)）、E×E（対角の Sp(2,Z)）の不変格子は布だけか（二次元ではヴェイユ現象が起きない理由）
"""
import itertools, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

pairs = list(itertools.combinations(range(4), 2)); names = ["e12", "e13", "e14", "e23", "e24", "e34"]
Q = np.zeros((6, 6), int)
for a, (i, j) in enumerate(pairs):
    for b, (k, l) in enumerate(pairs):
        if len({i, j, k, l}) == 4: Q[a, b] = int(round(np.linalg.det(np.eye(4)[[i, j, k, l]])))
q = lambda x, y: int(x @ Q @ y)
theta = np.array([1, 0, 0, 0, 0, 1])
# このトーラス E_i×E_i の布 NS と、超越格子 T
Pi = np.array([[1, 1j, 0, 0], [0, 0, 1, 1j]])
om = np.array([Pi[0, i]*Pi[1, j] - Pi[0, j]*Pi[1, i] for (i, j) in pairs])
Mc = np.vstack([Q @ om.real, Q @ om.imag])
is_cloth = lambda v: np.abs(Mc @ v).max() < 1e-9
NSb = np.array([[-1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1], [0, 0, -1, 1, 0, 0], [0, 1, 0, 0, 1, 0]]).T   # −e12, e34, −e14+e23, e13+e24
assert all(is_cloth(NSb[:, j]) for j in range(4))
GNS = NSb.T @ Q @ NSb
print("布の格子 NS のグラム行列（基底 −e12, e34, −e14+e23, e13+e24）:\n", GNS, "\n  符号:", np.sign(np.linalg.eigvalsh(GNS.astype(float))).astype(int))
def show(v): return " ".join(f"{c:+d}{n}" for c, n in zip(v, names) if c)

# 根（δ² = −2）を係数 |·| ≤ 1 の範囲で列挙
def roots(in_sub=None, bound=1):
    out = []
    for c in itertools.product(range(-bound, bound+1), repeat=6):
        v = np.array(c)
        if not v.any() or q(v, v) != -2: continue
        if in_sub is not None and not in_sub(v): continue
        out.append(v)
    return out
def reflect(x, d): return x + q(x, d)*d

def orbit(start, gens, max_coeff=40, max_depth=8):
    seen = {tuple(start): 0}; frontier = [start]
    for depth in range(1, max_depth+1):
        nxt = []
        for x in frontier:
            for d in gens:
                y = reflect(x, d)
                if np.abs(y).max() > max_coeff: continue
                t = tuple(y)
                if t not in seen: seen[t] = depth; nxt.append(y)
        frontier = nxt
        if not frontier: break
    return seen

# (A) 家族の万華鏡：θ^⊥ の根
rootsA = roots(in_sub=lambda v: q(v, theta) == 0)
print(f"\n(A) 家族の万華鏡：θ に直交する根（|係数|≤1）{len(rootsA)} 枚の鏡")
orbA = orbit(np.array([1, 0, 0, 0, 0, 0]), rootsA, max_coeff=12, max_depth=4)
vecsA = np.array(list(orbA.keys())); depA = np.array(list(orbA.values()))
clothA = np.array([is_cloth(v) for v in vecsA])
print(f"    e12 の像：{len(vecsA)} 本（深さ ≤4, |係数|≤12）。すべて Q(v,v)=0, Q(v,θ)=1:", all(q(v, v) == 0 and q(v, theta) == 1 for v in vecsA))
print(f"    そのうち、このトーラス E_i×E_i の布であるもの：{clothA.sum()} 本")
for v in vecsA[clothA][:12]: print("       ", show(v))
print("    残りは「別のトーラスの布」：像を一つ ", show(vecsA[~clothA][0]), " など")

# (B) このトーラスの万華鏡：NS の中の根
rootsB = [v for v in roots(bound=2) if is_cloth(v)]
print(f"\n(B) このトーラスの万華鏡：NS の中の根（|係数|≤2）{len(rootsB)} 枚の鏡。例：", ", ".join(show(v) for v in rootsB[:4]))
targets = {"横の繊維 H = e12": np.array([1, 0, 0, 0, 0, 0]), "縦の繊維 V = e34": np.array([0, 0, 0, 0, 0, 1]),
           "対角線 Δ = H+V−e14+e23": np.array([1, 0, -1, 1, 0, 1]), "i のグラフ Γ = H+V+e13+e24": np.array([1, 1, 0, 0, 1, 1]),
           "反対角線 H+V+e14−e23": np.array([1, 0, 1, -1, 0, 1])}
orbB = orbit(np.array([1, 0, 0, 0, 0, 0]), rootsB, max_coeff=20, max_depth=6)
print(f"    e12 の像：{len(orbB)} 本、すべて布:", all(is_cloth(np.array(t)) for t in orbB))
for nm, t in targets.items():
    hit = tuple(t) in orbB or tuple(-t) in orbB
    print(f"    {nm:18s} に届くか：{'届く（深さ ' + str(orbB.get(tuple(t), orbB.get(tuple(-t)))) + '）' if hit else '届かない'}")
# 届かない理由：W(NS) の不変量（θ との内積、Q(v,v)、2 で割った剰余類…）
print("    各布の (Q(v,v), Q(v,θ)):", {nm: (q(t, t), q(t, theta)) for nm, t in targets.items()})

# (C) 部分族の万華鏡（本物のモノドロミー行列）の不変格子
from scipy.linalg import null_space
def ext2(A):
    M = np.zeros((6, 6), int)
    for b_, (k, l) in enumerate(pairs):
        for a_, (i_, j_) in enumerate(pairs):
            M[a_, b_] = int(round(np.linalg.det(A[np.ix_([i_, j_], [k, l])])))
    return M
Tm = np.array([[1, 1], [0, 1]]); Sm = np.array([[0, -1], [1, 0]]); I2 = np.eye(2, dtype=int)
blk = lambda A, B: np.block([[A, np.zeros((2, 2), int)], [np.zeros((2, 2), int), B]])
gens1 = [ext2(blk(Tm, I2)), ext2(blk(Sm, I2)), ext2(blk(I2, Tm)), ext2(blk(I2, Sm))]      # E×E' : SL2(Z)×SL2(Z)
gens2 = [ext2(blk(Tm, Tm)), ext2(blk(Sm, Sm))]                                              # E×E  : 対角 SL2(Z)
def invariants(gs): return null_space(np.vstack([g - np.eye(6, dtype=int) for g in gs]).astype(float))
inv1 = invariants(gens1); inv2 = invariants(gens2)
print(f"\n(C) 部分族 E×E' の万華鏡（SL2(Z)×SL2(Z)）の不変格子：階数 {inv1.shape[1]}（e12, e34 ＝ 繊維、どちらも布）")
print(f"    部分族 E×E の万華鏡（対角 SL2(Z)）の不変格子：階数 {inv2.shape[1]}（e12, e34, Δ ＝ すべて布）")
print("    → 二次元では、万華鏡で動かない影はすべて布（レフシェッツ (1,1) 定理）。四次元のヴェイユ類は、動かないのに布が見つからない最初の例。")

# ---------- 図 ----------
fig = plt.figure(figsize=(15, 9), dpi=150, constrained_layout=True)
gs = fig.add_gridspec(2, 3)
# (A) 像を「このトーラスの布からの距離」平面に落とす
ax = fig.add_subplot(gs[:, :2])
xy = (Mc @ vecsA.T).T
jit = np.random.default_rng(0).normal(scale=0.08, size=xy.shape)
sc = ax.scatter(xy[~clothA, 0] + jit[~clothA, 0], xy[~clothA, 1] + jit[~clothA, 1], c=depA[~clothA], cmap="viridis", s=12, alpha=0.7, label="別のトーラスの布（万華鏡の像）")
ax.scatter(xy[clothA, 0], xy[clothA, 1], s=160, marker="*", color="#d9541e", edgecolor="k", zorder=5, label=f"このトーラス E_i×E_i の布（{clothA.sum()} 本、原点）")
cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01); cb.set_label("鏡に映した回数（深さ）")
ax.set_xlabel("Q(v, Re Ω)   （この値が 0 でないと、このトーラスでは布でない）"); ax.set_ylabel("Q(v, Im Ω)")
ax.set_title(f"(A) 家族の万華鏡：横の繊維 e12 の像 {len(vecsA)} 本。どれも「どこかのトーラスの布」だが、このトーラスの布は {clothA.sum()} 本だけ", fontsize=10, loc="left")
ax.legend(fontsize=9, loc="upper right"); ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[0, 2]); ax.axis("off")
rows = [[nm, ("届く" if (tuple(t) in orbB or tuple(-t) in orbB) else "届かない"), f"{q(t, t)}, {q(t, theta)}"] for nm, t in targets.items()]
tb = ax.table(cellText=rows, colLabels=["布（長さ² はすべて 0）", "H から鏡像で", "(v², v·θ)"], loc="center", cellLoc="center"); tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1, 1.7)
ax.set_title(f"(B) このトーラスの万華鏡（NS の中の鏡 {len(rootsB)} 枚）\ne12 の像 {len(orbB)} 本はすべて布。他の布に届くか", fontsize=10, loc="left")

ax = fig.add_subplot(gs[1, 2])
labs = ["E×E'\nSL2×SL2", "E×E\n対角 SL2", "ヴェイユ型 4 次元\n（実験 28）"]
inv = [inv1.shape[1], inv2.shape[1], 3]; alg = [inv1.shape[1], inv2.shape[1], 1]
x = np.arange(3)
ax.bar(x - 0.2, inv, 0.38, color="#9ec3ee", label="万華鏡で動かない影（不変なホッジ類）")
ax.bar(x + 0.2, alg, 0.38, color="#1f5fbf", label="そのうち布（または布の積）")
for i_, (a_, b_) in enumerate(zip(inv, alg)): ax.text(i_ - 0.2, a_ + 0.1, str(a_), ha="center"); ax.text(i_ + 0.2, b_ + 0.1, str(b_), ha="center")
ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=8.5); ax.set_ylim(0, 4.2); ax.legend(fontsize=8)
ax.set_title("(C) 万華鏡で動かない影は布か：二次元では全部、四次元で初めて外れる", fontsize=10, loc="left"); ax.grid(alpha=0.25, axis="y")
fig.suptitle("万華鏡（ピカール＝レフシェッツ反射）に一本の布を入れる ― 像はすべて「どこかのトーラスの布」。このトーラスの布は小さな万華鏡でつながり、四次元で初めてほどける", fontsize=11.5)
fig.savefig("figures/kaleidoscope.png")
