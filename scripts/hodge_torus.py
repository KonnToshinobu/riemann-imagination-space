"""
ホッジのトーラス：光で積分する定規
  X = C²/Λ（二次元複素トーラス ＝ 四次元実トーラス）, 周期行列 Π (2×4)
  H²(X,Z) = ∧²Λ*（階数 6）, 定規 Q(α,β) = ∫_X α∧β  → 符号 (3,3)
  影の分解:  H² = H^{2,0} ⊕ H^{1,1} ⊕ H^{0,2},   ω = dz₁∧dz₂ の成分 ＝ Π の 2×2 小行列式（プリュッカー座標）
  H^{1,1}_R = {v : Q(v, Re ω) = Q(v, Im ω) = 0}（実 4 次元）。Q をここに制限すると符号 (1,3) ＝ 光円錐（ホッジ＝リーマン）
  布（代数的サイクル）= NS(X) = H²(X,Z) ∩ H^{1,1}。本数 ρ（ピカール数）は格子の数論で決まる：
     一般のトーラス 0、主偏極アーベル曲面 1、E×E'（非同種）2、E×E 3、E_i×E_i（虚数乗法）4
  ホッジ指数定理：NS 上でも符号は (1, ρ−1)
"""
import itertools, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]
rng = np.random.default_rng(7)

pairs = list(itertools.combinations(range(4), 2))          # e_{ij}, i<j
Q = np.zeros((6, 6))
for a, (i, j) in enumerate(pairs):
    for b, (k, l) in enumerate(pairs):
        if len({i, j, k, l}) == 4:
            perm = [i, j, k, l]; sgn = np.linalg.det(np.eye(4)[perm])
            Q[a, b] = sgn
assert sorted(np.round(np.linalg.eigvalsh(Q))) == [-1, -1, -1, 1, 1, 1]

def int_grid(vmax):
    g = np.indices([2*vmax+1]*6).reshape(6, -1).T - vmax
    return g[np.any(g != 0, axis=1)]
_grids = {}
def analyze(Pi, vmax=3):
    om = np.array([Pi[0, i]*Pi[1, j] - Pi[0, j]*Pi[1, i] for (i, j) in pairs])      # ω のプリュッカー座標
    Qo = Q if np.real(om @ Q @ np.conj(om)) > 0 else -Q                             # 向き：∫ω∧ω̄ > 0 になるように
    M = np.vstack([Qo @ om.real, Qo @ om.imag])                                     # H^{1,1} の条件 (2×6)
    _, s_, Vt = np.linalg.svd(M); H11 = Vt[2:].T                                    # H^{1,1}_R の基底 6×4
    ev11 = np.linalg.eigvalsh(H11.T @ Qo @ H11)
    if vmax not in _grids: _grids[vmax] = int_grid(vmax)
    g = _grids[vmax]
    ok = np.abs(g @ M.T).max(axis=1) < 1e-8
    NS = g[ok].T.astype(float)
    rho = np.linalg.matrix_rank(NS) if NS.size else 0
    if rho:
        B = np.linalg.svd(NS)[0][:, :rho]; evNS = np.linalg.eigvalsh(B.T @ Qo @ B)
    else: evNS = np.array([])
    return ev11, evNS, rho

def E_x_E(t1, t2): return np.array([[1, t1, 0, 0], [0, 0, 1, t2]], complex)
def general(Z): return np.hstack([np.eye(2), Z])

tori = [
    ("一般のトーラス（周期はでたらめ）\n布なし・代数的ですらない", general(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))),
    ("主偏極アーベル曲面（Z 対称, Im Z > 0）\n一般のもの", general(np.array([[0.3+1.1j, 0.2+0.4j], [0.2+0.4j, -0.5+1.3j]]))),
    ("E_τ × E_τ'（同種でない二本の楕円曲線）", E_x_E(0.37+1.21j, -0.58+0.83j)),
    ("E_τ × E_τ（同じ楕円曲線を二本）", E_x_E(0.37+1.21j, 0.37+1.21j)),
    ("E_i × E_i（虚数乗法、τ = i）", E_x_E(1j, 1j)),
]
res = []
for name, Pi in tori:
    ev11, evNS, rho = analyze(Pi)
    res.append((name, ev11, evNS, rho))
    print(f"{name.splitlines()[0]:38s} 影の定規 H^{{1,1}} の符号 (+{np.sum(ev11>0)}, −{np.sum(ev11<0)})   布の本数 ρ = {rho}   布の定規 NS の符号 (+{np.sum(evNS>1e-9)}, −{np.sum(evNS<-1e-9)})")

# τ' の平面で、E_τ × E_τ' の布の本数の星座：τ' = (aτ+b)/(cτ+d)（E_τ と同種）なら ρ = 3
tau0 = 0.37+1.21j
pts = {}
for a, b, c, d in itertools.product(range(-3, 4), repeat=4):
    if a*d - b*c == 0: continue
    tau = (a*tau0 + b)/(c*tau0 + d)
    if not (-1.5 <= tau.real <= 1.5 and 0.15 <= tau.imag <= 2.6): continue
    key = (round(tau.real, 6), round(tau.imag, 6))
    if key in pts: continue
    pts[key] = (abs(a*d - b*c), analyze(E_x_E(tau0, tau), vmax=4)[2])
gen = [analyze(E_x_E(tau0, complex(rng.uniform(-1.5, 1.5), rng.uniform(0.2, 2.5))), vmax=4)[2] for _ in range(12)]
keys = np.array(list(pts.keys())); degs = np.array([v[0] for v in pts.values()]); rhos = np.array([v[1] for v in pts.values()])
print("E_τ × E_τ'：τ' がでたらめのとき ρ =", [int(x) for x in gen], "   τ' が E_τ と同種のとき ρ の分布 =", {int(k): int(v) for k, v in zip(*np.unique(rhos, return_counts=True))})
cm = [analyze(E_x_E(1j, (a+1j*b)/c), vmax=4)[2] for a, b, c in [(1,1,1),(1,2,1),(1,1,2),(2,3,1),(1,3,2)]]
print("参考 E_i × E_τ'（τ' ∈ Q(i)、両方とも虚数乗法）のとき ρ =", [int(x) for x in cm])

# ---------- 図 ----------
fig = plt.figure(figsize=(15, 9.5), dpi=150, constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1])
ax = fig.add_subplot(gs[0, 0])
y0 = 0
yt, yl = [], []
for name, ev11, evNS, rho in res:
    xs = np.arange(len(ev11)); ax.bar(xs - 0.2, ev11, 0.38, color=["#d9541e" if e > 0 else "#1f5fbf" for e in ev11], bottom=y0)
    if len(evNS): ax.bar(np.arange(len(evNS)) + 0.2, evNS, 0.38, color=["#d9541e" if e > 0 else "#1f5fbf" for e in evNS], bottom=y0, hatch="//", edgecolor="white")
    ax.axhline(y0, color="#999", lw=0.6)
    ax.text(4.2, y0 + 0.2, f"布 ρ = {rho}本\n影 (+1, −3)  布 (+{int(np.sum(evNS>1e-9))}, −{int(np.sum(evNS<-1e-9))})", fontsize=9, va="bottom")
    yt.append(y0); yl.append(name); y0 += 3.6
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=8.5)
ax.set_xticks(range(4)); ax.set_xticklabels(["固有値 1", "2", "3", "4"]); ax.set_xlim(-0.6, 6.2)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#d9541e", label="正（時間方向＝光／ケーラー類）"), Patch(color="#1f5fbf", label="負（空間方向）"),
                   Patch(facecolor="#777", hatch="//", edgecolor="white", label="斜線＝布（代数的サイクル NS）に制限した定規")], fontsize=8, loc="lower right")
ax.set_ylim(-1.6, y0 + 0.6); ax.set_title("光で積分する定規 Q(α,β)=∫α∧β：影 H^{1,1} ではどのトーラスでも (+1, −3) の光円錐。\n布 NS に制限しても (+1, −(ρ−1))（ホッジ指数定理）。布の本数 ρ だけが格子の数論で変わる", fontsize=10, loc="left")

ax = fig.add_subplot(gs[0, 1])
m3 = rhos == 3; m2 = rhos == 2
ax.scatter(keys[m3, 0], keys[m3, 1], s=120/np.sqrt(degs[m3]) + 8, color="#1f5fbf", label="ρ = 3：τ' が E_τ と同種（τ' = (aτ+b)/(cτ+d)）、大きい点ほど低い次数", zorder=3)
if m2.any(): ax.scatter(keys[m2, 0], keys[m2, 1], s=10, color="#bbb", label="同種だが探索範囲（|係数| ≤ 4）を越えて布が見つからず", zorder=2)
rt = np.array([[rng.uniform(-1.5, 1.5), rng.uniform(0.2, 2.5)] for _ in range(12)])
ax.scatter(rt[:, 0], rt[:, 1], s=30, marker="x", color="#d9541e", label="τ' がでたらめ → ρ = 2", zorder=4)
ax.scatter([tau0.real], [tau0.imag], s=160, marker="*", color="#c99a2e", edgecolor="k", label="τ 自身（E_τ × E_τ, ρ = 3）", zorder=5)
ax.set_xlim(-1.5, 1.5); ax.set_ylim(0, 2.6); ax.set_aspect("equal")
ax.set_xlabel("Re τ'"); ax.set_ylabel("Im τ'")
ax.set_title("E_τ × E_τ' の布の本数：τ' の平面に浮かぶ星座\n布が一本増えるのは、τ' が τ と数論的に結ばれている（同種）ときだけ", fontsize=10, loc="left")
ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1); ax.grid(alpha=0.25)
fig.suptitle("ホッジのトーラス：定規は「数える」のではなく「光で積分する」。定規の符号は幾何で決まり（常に光円錐）、布の本数は数論で決まる", fontsize=12)
fig.savefig("figures/hodge_torus.png")
