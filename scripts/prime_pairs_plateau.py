"""
A. 素数対 (p, p+h) を数え、ハーディ＝リトルウッドのリズム S(h) が出るか
B. 形状因子を「素数側」から厳密に組み立て（明示公式＝恒等式）、対角（同じ素数同士）と非対角（異なる素数同士）に分ける
C. 非対角項を「Λ(n)Λ(n+h) → S(h)（ハーディ＝リトルウッドの平均）」で置き換えた模型がプラトーを出すか（ボゴモルニー＝キーティング 1996）
"""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

# ---------- A. 素数対を数える ----------
NA = 10_000_000
isp = np.ones(NA+1, bool); isp[:2] = False
for i in range(2, int(NA**0.5)+1):
    if isp[i]: isp[i*i::i] = False
H = 60
counts = np.array([np.count_nonzero(isp[:-h] & isp[h:]) for h in range(1, H+1)])
# ハーディ＝リトルウッド: S(h) = 0 (h 奇数),  2 C2 Π_{p|h, p>2} (p-1)/(p-2) (h 偶数)
small_p = [p for p in range(3, 200) if isp[p]]
C2 = np.prod([1 - 1/(p-1)**2 for p in small_p] + [1 - 1/(p-1)**2 for p in range(211, 20000) if isp[p]])
def S(h):
    if h % 2: return 0.0
    v = 2*C2
    for p in small_p:
        if h % p == 0: v *= (p-1)/(p-2)
    return v
Sh = np.array([S(h) for h in range(1, H+1)])
ratio_obs = counts/counts[1]                   # h=2 を 1 とする
ratio_HL = Sh/Sh[1]
print("A. 素数対 (p ≤ 1e7):  h=2:", counts[1], " h=4:", counts[3], " h=6:", counts[5], " h=30:", counts[29], " h=3:", counts[2])
print("   観測/予言 (偶数 h):", np.round((ratio_obs/np.where(ratio_HL > 0, ratio_HL, np.nan))[1::2], 3))

# ---------- B. 形状因子を素数側から ----------
z = np.array(json.load(open("data/zeros_1500.json")))
a, b = 600.0, z[-1]
zz = z[(z >= a) & (z <= b)]
w_of = lambda r: 0.5*(1 - np.cos(2*np.pi*(r-a)/(b-a)))
wz = w_of(zz)
Lr = np.log(np.sqrt(a*b)/(2*np.pi))
taus = np.linspace(0.3, 2.0, 851)
ts = taus*Lr
K_zeros = np.array([abs(np.sum(wz*np.exp(1j*t*zz)))**2 for t in ts])/np.sum(wz**2)

# ŵ(v) = (1/2π) ∫ w(r) e^{-ivr} dr  を数値で作って補間
rg = np.linspace(a, b, 40001); wr = w_of(rg)
vg = np.linspace(-0.12, 0.12, 4801)
what = np.array([np.trapezoid(wr*np.exp(-1j*v*rg), rg) for v in vg])/(2*np.pi)
def W(v):
    v = np.asarray(v); out = np.zeros(v.shape, complex); m = np.abs(v) < 0.12
    out[m] = np.interp(v[m], vg, what.real) + 1j*np.interp(v[m], vg, what.imag)
    return out

NB = int(np.exp(2.0*Lr + 0.3)) + 2
lam = np.zeros(NB+1); sv = np.ones(NB+1, bool); sv[:2] = False
for i in range(2, int(NB**0.5)+1):
    if sv[i]: sv[i*i::i] = False
for p in np.nonzero(sv)[0]:
    pk = p
    while pk <= NB: lam[pk] = np.log(p); pk *= p
n_all = np.arange(1, NB+1); logn = np.log(n_all); amp = lam[1:]/np.sqrt(n_all)

K_pr, K_diag, K_off = [], [], []
for t in ts:
    m = np.abs(logn - t) < 0.12
    c = amp[m]*W(logn[m] - t)                       # 各 n の寄与（複素）
    tot = abs(np.sum(c))**2; dia = np.sum(np.abs(c)**2)
    K_pr.append(tot); K_diag.append(dia); K_off.append(tot - dia)
norm = np.sum(wz**2)
K_pr = np.array(K_pr)/norm; K_diag = np.array(K_diag)/norm; K_off = np.array(K_off)/norm
print("B. 恒等式の検算 max|K_zeros − K_primes| / max K =", np.max(np.abs(K_zeros - K_pr))/K_zeros.max())

# ---------- C. 非対角をハーディ＝リトルウッドで置き換える ----------
HC = 400
Sh_C = np.array([S(h) for h in range(1, HC+1)])
K_HL = []
for t in ts:
    m = np.abs(logn - t) < 0.12
    nn = n_all[m]; Wn = W(logn[m] - t)
    acc = 0.0
    for h in range(2, HC+1, 2):                     # 奇数 h は S(h)=0
        if h >= len(nn): break
        # Λ(n)Λ(n+h) → S(h);  1/√(n(n+h)) と窓の積を n について足す
        acc += Sh_C[h-1]*np.sum(2*np.real(Wn[:-h]*np.conj(Wn[h:]))/np.sqrt(nn[:-h]*nn[h:]))
    K_HL.append(acc)
K_HL = np.array(K_HL)/norm

def smooth(y, n=61):
    k = np.ones(n)/n; return np.convolve(y, k, mode="same")
sm = slice(40, -40)
m1 = (taus > 1.2) & (taus < 1.9)
print("C. プラトー領域 1.2<τ<1.9:  <K_zeros> =", np.mean(smooth(K_zeros)[m1]), " <対角> =", np.mean(smooth(K_diag)[m1]),
      " <非対角(厳密)> =", np.mean(smooth(K_off)[m1]), " <非対角(HL模型)> =", np.mean(smooth(K_HL)[m1]))
m0 = (taus > 0.4) & (taus < 0.9)
print("   ランプ領域 0.4<τ<0.9:  <非対角(厳密)> =", np.mean(smooth(K_off)[m0]), " <非対角(HL模型)> =", np.mean(smooth(K_HL)[m0]))

# ---------- 図 ----------
fig, axes = plt.subplots(2, 1, figsize=(13, 10.5), dpi=150)
ax = axes[0]
hs = np.arange(1, H+1)
ax.bar(hs, ratio_obs, color="#9ec3ee", edgecolor="white", label="数えた素数対の個数（h=2 を 1 とする）  p ≤ 1000万")
ax.plot(hs, ratio_HL, "o", color="#d9541e", ms=5, label="ハーディ＝リトルウッドのリズム S(h)/S(2)（1923 年の予想）")
for h in [6, 12, 30, 42, 60]:
    ax.text(h, ratio_obs[h-1] + 0.08, str(h), ha="center", fontsize=8)
ax.set_xlabel("素数対の距離 h"); ax.set_ylabel("現れやすさ")
ax.set_title("A. 素数のリズム：距離 h の素数対 (p, p+h) はどれだけ一緒に現れるか  —  奇数は 0、6 の倍数で多い", fontsize=10.5, loc="left")
ax.legend(fontsize=9); ax.grid(alpha=0.25, axis="y"); ax.set_xlim(0, H+1)

ax = axes[1]
ax.plot(taus[sm], smooth(K_zeros)[sm], color="#1f5fbf", lw=2.6, label="零点から測った形状因子")
ax.plot(taus[sm], smooth(K_pr)[sm], color="k", lw=1.0, ls=":", label="素数から厳密に組み立てたもの（明示公式、全部の n）")
ax.plot(taus[sm], smooth(K_diag)[sm], color="#d9541e", lw=2.0, ls="--", label="対角：同じ素数同士の干渉（ワームホール＝ランプ）")
ax.plot(taus[sm], smooth(K_off)[sm], color="#2a9d4f", lw=2.0, label="非対角（厳密）：異なる素数同士の干渉")
ax.plot(taus[sm], smooth(K_HL)[sm], color="#8a5fbf", lw=2.0, ls="--", label="非対角（模型）：Λ(n)Λ(n+h) をリズム S(h) に置き換え")
ax.plot(taus, np.minimum(taus, 1), color="#888", lw=0.8, ls=":")
ax.axhline(0, color="k", lw=0.6); ax.axvline(1, color="k", lw=0.6)
ax.set_xlim(0.35, 1.95); ax.set_ylim(-1.2, 2.0)
ax.set_xlabel("τ"); ax.set_ylabel("K(τ)")
ax.set_title("B/C. プラトーは、素数対のリズム S(h) が作る非対角の打ち消しで出来ている", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, loc="upper left", ncol=2); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/prime_pairs_plateau.png")
