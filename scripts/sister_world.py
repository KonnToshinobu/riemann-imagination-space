"""
証明ができた兄弟の世界：有限体上の楕円曲線 E: y² = x³ + x + 1
  各素数 p で点の個数 #E(F_p) = p + 1 − a_p を数える。
  E の「ゼータ」の零点は  α, ᾱ = (a_p ± √(a_p² − 4p)) / 2,   零点の実部 1/2 ⇔ |α| = √p ⇔ |a_p| ≤ 2√p（ハッセ）
  証明の心臓（ホッジ指数定理＝光円錐）：
    曲面 E×E 上の 4 本の曲線 H（横の繊維）, V（縦の繊維）, Δ（対角線）, Γ（フロベニウスのグラフ）の交点数
      H²=V²=0, H·V=1, Δ²=0, Γ²=0, Δ·H=Δ·V=1, Γ·H=1, Γ·V=p, Δ·Γ=#E(F_p)=p+1−a_p
    交点形式の符号は (+,−,−,−)：時間方向が一本（H+V）、残りは空間方向。
    時間方向に直交する Δ'=Δ−H−V, Γ'=Γ−pH−V の Gram 行列 [[−2, −a_p],[−a_p, −2p]] が負半定値
      ⇔ 4p − a_p² ≥ 0 ⇔ |a_p| ≤ 2√p   ← リーマン予想（この世界版）
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

A, B = 1, 1
PMAX = 3000
sv = np.ones(PMAX+1, bool); sv[:2] = False
for i in range(2, int(PMAX**0.5)+1):
    if sv[i]: sv[i*i::i] = False
primes = [int(p) for p in np.nonzero(sv)[0] if p > 3 and (4*A**3 + 27*B**2) % p != 0]

def count_points(p):
    x = np.arange(p, dtype=np.int64)
    f = (x**3 + A*x + B) % p
    # ルジャンドル記号 (f/p) をオイラーの判定法で
    leg = np.array([pow(int(v), (p-1)//2, p) for v in f])
    leg = np.where(leg == p-1, -1, leg)
    return p + 1 + int(np.sum(leg))        # 無限遠点を含む

ap, Np = [], []
for p in primes:
    N = count_points(p); ap.append(p + 1 - N); Np.append(N)
ap = np.array(ap, float); ps = np.array(primes, float)
print(f"素数 {len(ps)} 個 (p ≤ {PMAX}):  |a_p| ≤ 2√p の反例 = {np.sum(np.abs(ap) > 2*np.sqrt(ps))} 個   max |a_p|/(2√p) = {np.max(np.abs(ap)/(2*np.sqrt(ps))):.4f}")

# 零点
disc = ap**2 - 4*ps
alpha = (ap + np.sqrt(disc.astype(complex)))/2
# 光円錐：4×4 Gram 行列の固有値の符号
sig = []
for p, a in zip(ps, ap):
    G = np.array([[0, 1, 1, 1], [1, 0, 1, p], [1, 1, 0, p+1-a], [1, p, p+1-a, 0]], float)   # H, V, Δ, Γ
    ev = np.linalg.eigvalsh(G); sig.append((np.sum(ev > 1e-9), np.sum(ev < -1e-9)))
sig = np.array(sig)
print("交点形式の符号 (正の個数, 負の個数) の種類:", set(map(tuple, sig)))

fig = plt.figure(figsize=(15, 10), dpi=150)
gs = fig.add_gridspec(2, 3)

ax = fig.add_subplot(gs[0, :2])
ax.scatter(ps, ap/(2*np.sqrt(ps)), s=6, color="#1f5fbf", label="a_p / (2√p)   （E: y² = x³ + x + 1 の点を数えた結果）")
ax.axhline(1, color="#d9541e", lw=2); ax.axhline(-1, color="#d9541e", lw=2)
ax.fill_between([0, PMAX], -1, 1, color="#d9541e", alpha=0.07)
ax.set_xlim(0, PMAX); ax.set_ylim(-1.3, 1.3)
ax.set_xlabel("素数 p"); ax.set_ylabel("a_p / (2√p)")
ax.set_title(f"① ハッセの定理（この世界のリーマン予想）：{len(ps)} 個の素数すべてで |a_p| ≤ 2√p。反例 0", fontsize=10.5, loc="left")
ax.legend(fontsize=9, loc="lower right"); ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[0, 2])
th = np.angle(alpha)
ax.hist(np.abs(th), bins=30, range=(0, np.pi), density=True, color="#9ec3ee", edgecolor="white", label="零点の角度 θ（α = √p e^{iθ}）")
t = np.linspace(0, np.pi, 200); ax.plot(t, 2/np.pi*np.sin(t)**2, color="k", lw=1.8, label="佐藤＝テイト分布 (2/π) sin²θ")
ax.set_xlabel("θ"); ax.set_title("零点の角度の分布", fontsize=10.5, loc="left")
ax.legend(fontsize=8); ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[1, 0])
z = alpha/np.sqrt(ps)
ax.scatter(z.real, z.imag, s=5, color="#1f5fbf", label="零点 α/√p（本物の a_p）")
ax.scatter(z.real, -z.imag, s=5, color="#1f5fbf")
circ = np.exp(1j*np.linspace(0, 2*np.pi, 400)); ax.plot(circ.real, circ.imag, color="#d9541e", lw=1.5, label="半径 √p の円（＝実部 1/2）")
# 偽の a_p（ハッセを破る）を入れたら
fake_a = 2.3*np.sqrt(ps[:60]); fd = fake_a**2 - 4*ps[:60]
fz = ((fake_a + np.sqrt(fd))/2)/np.sqrt(ps[:60]); fz2 = ((fake_a - np.sqrt(fd))/2)/np.sqrt(ps[:60])
ax.scatter(np.r_[fz, fz2], np.zeros(120), s=12, marker="x", color="#555", label="もし |a_p| = 2.3√p なら：零点が円から外れる")
ax.set_aspect("equal"); ax.set_xlim(-1.8, 1.8); ax.set_ylim(-1.3, 1.3)
ax.set_title("② 零点はすべて半径 √p の円の上", fontsize=10.5, loc="left")
ax.legend(fontsize=7.5, loc="lower left"); ax.grid(alpha=0.25)

ax = fig.add_subplot(gs[1, 1:])
det = 4*ps - ap**2
ax.scatter(ps, det/(4*ps), s=6, color="#2a9d4f", label="(4p − a_p²)/4p ＝ 光円錐の内側にいる度合い（0 以上が必要）")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlim(0, PMAX); ax.set_ylim(-0.1, 1.05)
ax.set_xlabel("素数 p"); ax.set_ylabel("det / 4p")
ax.set_title(f"③ ホッジ指数定理（光円錐）：交点形式の符号は常に (+1, −3)。時間軸に直交する Gram 行列 [[−2, −a_p],[−a_p, −2p]] の行列式 4p − a_p² ≥ 0", fontsize=9.5, loc="left")
ax.legend(fontsize=9, loc="lower right"); ax.grid(alpha=0.25)
fig.suptitle("証明ができた兄弟の世界：有限体上の楕円曲線では、光円錐（ホッジ指数定理）が零点を半径 √p の円に縛る", fontsize=12)
fig.tight_layout()
fig.savefig("figures/sister_world.png")
