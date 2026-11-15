"""
石鹸膜の実演：整数の影の縁に膜は張れるか（E_i × E_i, 平坦計量）
  座標 u = (x1, y1, x2, y2), 格子 Z⁴, 光 ω = dx1∧dy1 + dx2∧dy2, 複素構造 J: (x,y) ↦ (−y, x)
  影のクラス v ∈ ∧² Z⁴（6 成分 v_ij）。平らな 2-トーラス span(a,b) の面積 = |a∧b|（プリュッカー座標のノルム）。
  膜の最小面積（質量ノルム）：∧²R⁴ = Λ⁺ ⊕ Λ⁻（自己双対／反自己双対）に分けて  mass(v) = √2·max(|v⁺|, |v⁻|)
      （分解可能 ⇔ |v⁺| = |v⁻| ⇔ 一枚の平らなトーラスで表せる。そのとき mass = 面積）
  ヴィルティンガー：mass(v) ≥ |∫_v ω|、等号 ⇔ 複素曲線（布）が代表。
  光円錐：v² = Q(v,v)。θ = ω 方向が時間的。
  予言：(1,1) 類（布の候補）のうち、光円錐の内側 (v² ≥ 0, ∫ω > 0) ⇔ 一枚の膜で等号。
        外側（空間的、原始的）⇔ 等号不可能、二枚の膜の差としてだけ現れる。
        (1,1) でない類 ⇔ 等号は決して成り立たない（光に届かない膜）。
"""
import itertools, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

pairs = list(itertools.combinations(range(4), 2))        # (0,1)=x1y1, (0,2)=x1x2, (0,3)=x1y2, (1,2)=y1x2, (1,3)=y1y2, (2,3)=x2y2
names = ["dx1∧dy1", "dx1∧dx2", "dx1∧dy2", "dy1∧dx2", "dy1∧dy2", "dx2∧dy2"]
def wedge(a, b): return np.array([a[i]*b[j] - a[j]*b[i] for (i, j) in pairs], float)
omega = np.array([1, 0, 0, 0, 0, 1.])
# 交点形式 Q（∧²×∧² → ∧⁴）、ホッジ星 ⋆、自己双対分解
Q = np.zeros((6, 6))
for a_, (i, j) in enumerate(pairs):
    for b_, (k, l) in enumerate(pairs):
        if len({i, j, k, l}) == 4: Q[a_, b_] = np.linalg.det(np.eye(4)[[i, j, k, l]])
star = Q.copy()                                           # 正規直交基底では ⋆ の行列 = Q
Pp, Pm = (np.eye(6) + star)/2, (np.eye(6) - star)/2
def mass(v): return np.sqrt(2)*max(np.linalg.norm(Pp @ v), np.linalg.norm(Pm @ v))
def light(v): return float(omega @ Q @ v) if False else float(v[0] + v[5])      # ∫ ω = v_{x1y1} + v_{x2y2}
J = np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]], float)
C = np.array([[np.linalg.det(J[np.ix_(list(I), list(K))]) for K in pairs] for I in pairs])   # ∧²J（ヴェイユ作用素）
is11 = lambda v: np.allclose(C @ v, v)
q = lambda v: float(v @ Q @ v)
assert np.allclose(omega @ Q @ omega, 2) and np.isclose(mass(omega), 2) and np.isclose(light(omega), 2)

# 名のある類：張る 2 本の格子ベクトルから
span = {
    "横の繊維 H = E×{pt}":        [(1, 0, 0, 0), (0, 1, 0, 0)],
    "縦の繊維 V = {pt}×E":        [(0, 0, 1, 0), (0, 0, 0, 1)],
    "対角線 Δ = {(z,z)}":        [(1, 0, 1, 0), (0, 1, 0, 1)],
    "i のグラフ {(z,iz)}":        [(1, 0, 0, 1), (0, 1, -1, 0)],
    "共役のグラフ {(z,z̄)}":      [(1, 0, 1, 0), (0, 1, 0, -1)],
    "x1x2 平面 {(x, x')}":        [(1, 0, 0, 0), (0, 0, 1, 0)],
}
classes = {k: wedge(np.array(a, float), np.array(b, float)) for k, (a, b) in span.items()}
H, V, D = classes["横の繊維 H = E×{pt}"], classes["縦の繊維 V = {pt}×E"], classes["対角線 Δ = {(z,z)}"]
classes["H − V（布の差）"] = H - V
classes["Δ − H − V（対角線の原始部分）"] = D - H - V
classes["Re(dz1∧dz2)（純粋な影）"] = np.array([0, 1, 0, 0, -1, 0.])
print(f"{'類':34s} {'(1,1)?':7s} {'光 ∫ω':>7s} {'膜の最小面積':>12s} {'v²':>5s}  判定")
rows = []
for k, v in classes.items():
    m, l, s = mass(v), light(v), q(v)
    dec = np.isclose(np.linalg.norm(Pp @ v), np.linalg.norm(Pm @ v))
    if is11(v) and np.isclose(m, abs(l)) and l > 0: verdict = "布：一枚の膜で等号"
    elif is11(v) and l >= -1e-9 and m > abs(l) + 1e-9: verdict = "布の候補だが等号不可：二枚の膜の差"
    elif is11(v): verdict = "布（向きが逆）"
    else: verdict = "布でない：どの膜も光に届かない"
    rows.append((k, is11(v), l, m, s, dec, verdict))
    print(f"{k:34s} {'○' if is11(v) else '×':7s} {l:7.2f} {m:12.3f} {s:5.0f}  {verdict}{'' if dec else '（一枚の平らなトーラスでは表せない）'}")

# すべての小さい整数類で：光円錐と等号の関係
pts = []
for c in itertools.product(range(-2, 3), repeat=6):
    v = np.array(c, float)
    if not v.any(): continue
    pts.append((light(v), mass(v), q(v), is11(v)))
pts = np.array(pts)
h = pts[:, 3] == 1
eq = np.isclose(pts[:, 1], np.abs(pts[:, 0]))
print(f"\n|係数|≤2 の整数類 {len(pts)} 個：(1,1) 型 {h.sum()} 個。")
print("  等号 mass = |∫ω| が成り立つのは (1,1) 型だけか：", bool(np.all(h[eq])))
print("  (1,1) 型のうち等号が成り立つもの ⇔ 光円錐の内側 (v² ≥ 0)：", bool(np.all((pts[h, 2] >= -1e-9) == eq[h])))
print("  (1,1) 型で v² < 0（空間的）のものは全部、等号不可（差でしか作れない）：", bool(np.all(~eq[h & (pts[:, 2] < 0)])))

# ---------- 図 ----------
fig = plt.figure(figsize=(15, 9), dpi=150, constrained_layout=True)
gs = fig.add_gridspec(2, 3)
ax = fig.add_subplot(gs[:, :2])
jit = np.random.default_rng(1).normal(scale=0.03, size=(len(pts), 2))
nh = ~h
ax.scatter(np.abs(pts[nh, 0]) + jit[nh, 0], pts[nh, 1] + jit[nh, 1], s=6, color="#bbb", alpha=0.5, label="(1,1) でない類：光に届かない（線の上に来ない）")
sp = h & (pts[:, 2] < -1e-9); tl = h & (pts[:, 2] >= -1e-9)
ax.scatter(np.abs(pts[sp, 0]) + jit[sp, 0], pts[sp, 1] + jit[sp, 1], s=14, color="#1f5fbf", alpha=0.8, label="布の候補・光円錐の外（v²<0）：等号不可 → 二枚の膜の差")
ax.scatter(np.abs(pts[tl, 0]), pts[tl, 1], s=22, color="#d9541e", alpha=0.9, label="布の候補・光円錐の内（v²≥0）：一枚の膜で等号（＝布）")
xx = np.linspace(0, 6, 10); ax.plot(xx, xx, "k-", lw=1.2, label="石鹸膜の下限  面積 = |∫ω|（ヴィルティンガー）")
offs = {"横の繊維 H = E×{pt}": (8, -12), "対角線 Δ = {(z,z)}": (8, -12), "共役のグラフ {(z,z̄)}": (8, 6), "x1x2 平面 {(x, x')}": (8, 2),
        "H − V（布の差）": (8, 18), "Re(dz1∧dz2)（純粋な影）": (8, 30)}
for k, is11_, l, m, s, dec, verdict in rows:
    if k in offs: ax.annotate(k.split("（")[0].split(" =")[0], (abs(l), m), textcoords="offset points", xytext=offs[k], fontsize=7.5, color="#333",
                              arrowprops=dict(arrowstyle="-", color="#999", lw=0.6))
ax.set_xlim(-0.3, 5.2); ax.set_ylim(0, 7.5)
ax.set_xlabel("光の積分 |∫_v ω|   （布なら、これが面積）"); ax.set_ylabel("その影を縁に持つ膜の最小面積  mass(v)")
ax.set_title("石鹸膜は光の下限に届くか：|係数|≤2 の整数類 15,624 個（E_i×E_i, 平坦計量）", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="upper left"); ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[0, 2])
labs = ["Δ\n(面積 2)", "− H\n(1)", "− V\n(1)", "= Δ−H−V\n光 0, 面積 2"]
vals = [2, -1, -1, 0]
ax.bar(range(4), [2, 1, 1, 2], color=["#d9541e", "#d9541e", "#d9541e", "#1f5fbf"], alpha=[1, 0.5, 0.5, 1] if False else 1)
ax.bar(range(3), [2, -1, -1], color=["#d9541e", "#9ec3ee", "#9ec3ee"])
ax.bar([3], [2], color="#1f5fbf")
ax.axhline(0, color="k", lw=0.8); ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=8.5)
ax.set_ylabel("光の積分 ∫ω（＝膜の面積、符号つき）")
ax.set_title("空間的な布は二枚の差：\n光 2 の膜から光 1+1 の膜を引く → 光 0、影だけ残る", fontsize=10, loc="left"); ax.grid(alpha=0.25, axis="y")

ax = fig.add_subplot(gs[1, 2])
cats = ["時間的 (v²>0)\n例 Δ, H+V", "光的 (v²=0)\n例 H, V", "空間的 (v²<0)\n例 H−V, Δ−H−V"]
gap = [np.mean(pts[h & (pts[:, 2] > 1e-9), 1] - np.abs(pts[h & (pts[:, 2] > 1e-9), 0])),
       np.mean(pts[h & np.isclose(pts[:, 2], 0), 1] - np.abs(pts[h & np.isclose(pts[:, 2], 0), 0])),
       np.mean(pts[sp, 1] - np.abs(pts[sp, 0]))]
ax.bar(cats, gap, color=["#d9541e", "#c99a2e", "#1f5fbf"])
for i_, g_ in enumerate(gap): ax.text(i_, g_ + 0.03, f"{g_:.2f}", ha="center")
ax.set_ylabel("膜の面積 − 光の積分（平均）"); ax.set_title("光円錐で分かれる：内側は一枚の膜、外側は差", fontsize=10, loc="left"); ax.grid(alpha=0.25, axis="y")
fig.suptitle("石鹸膜の実演 ― 光円錐の内側の布は一枚の膜で光の下限に届き、外側（原始的）の布は二枚の膜の差としてだけ現れる。(1,1) でない影には届く膜がない", fontsize=11.5)
fig.savefig("figures/soap_films.png")
