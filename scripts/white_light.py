"""
白色光とコヒーレント光（指標の直交性）
  素数の波  S_a(t) = Σ_{n ≡ a (mod q)} Λ(n) n^{-1/2} w(log n) e^{-it log n}    （余り a の素数だけの螺旋）
  色なし（ゼータ、コヒーレント）:  S_0 = Σ_a S_a            → |S_0|² = Σ_a |S_a|² + Σ_{a≠b} S_a S̄_b
  色つき（指標 χ）:               S_χ = Σ_a χ(a) S_a
  直交性:  (1/φ(q)) Σ_χ |S_χ|² = Σ_a |S_a|²    ← 色を全部足す（白色光）と、違う余り同士の干渉 Σ_{a≠b} が消える
"""
import numpy as np, matplotlib, json
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

q = 5
chars = {  # mod 5 の 4 つの指標（生成元 2: 2^0=1, 2^1=2, 2^2=4, 2^3=3）
    "χ0（色なし）": {1: 1, 2: 1, 3: 1, 4: 1},
    "χ1": {1: 1, 2: 1j, 4: -1, 3: -1j},
    "χ2（実）": {1: 1, 2: -1, 4: 1, 3: -1},
    "χ3": {1: 1, 2: -1j, 4: -1, 3: 1j},
}
NMAX = 200000
lam = np.zeros(NMAX+1); sv = np.ones(NMAX+1, bool); sv[:2] = False
for i in range(2, int(NMAX**0.5)+1):
    if sv[i]: sv[i*i::i] = False
for p in np.nonzero(sv)[0]:
    pk = p
    while pk <= NMAX: lam[pk] = np.log(p); pk *= p
n = np.arange(1, NMAX+1); L = np.log(n)
b = 2.0
w = np.exp(-L**2/(4*b))                      # 窓（明示公式の試験関数と同じ形）
amp = lam[1:]/np.sqrt(n)*w
ts = np.linspace(0.5, 45, 1800)
S = {}
for a in [1, 2, 3, 4]:
    m = (n % q == a)
    S[a] = np.array([np.sum(amp[m]*np.exp(-1j*t*L[m])) for t in ts])
S0 = sum(S.values())                                             # コヒーレント（ゼータ、5 の倍数を除く）
Schi = {name: sum(c[a]*S[a] for a in [1, 2, 3, 4]) for name, c in chars.items()}
incoh = sum(np.abs(S[a])**2 for a in [1, 2, 3, 4])               # 各色の強さの和（白色光）
family = np.mean([np.abs(v)**2 for v in Schi.values()], axis=0)  # 族の平均
cross = np.abs(S0)**2 - incoh                                    # 違う余り同士の干渉
print("直交性の検算  max|族の平均 − Σ_a|S_a|²| =", np.max(np.abs(family - incoh)))

zeros = np.array(json.load(open("data/zeros_300.json")))
fig, axes = plt.subplots(2, 1, figsize=(13, 9), dpi=150, sharex=True)
ax = axes[0]
ax.plot(ts, np.abs(S0)**2, color="#1f5fbf", lw=2.2, label="コヒーレント光（色なし＝ゼータ）：|Σ_a S_a|²  — 違う余りの素数同士も干渉する")
ax.plot(ts, family, color="#c99a2e", lw=2.2, ls="--", label="白色光（4 色の族の平均）：(1/4) Σ_χ |S_χ|²  ＝ Σ_a |S_a|²  — 干渉が消える")
for z in zeros[zeros < 45]:
    ax.axvline(z, color="k", lw=0.6, ls=":")
ax.set_ylabel("光の強さ"); ax.set_title("色を全部足すと干渉縞が消える（指標の直交性）：点線＝ゼータの零点", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5); ax.grid(alpha=0.25)

ax = axes[1]
ax.plot(ts, cross, color="#2a9d4f", lw=1.8, label="干渉項 Σ_{a≠b} S_a S̄_b ＝ コヒーレント − 白色光（違う余りの素数対 p, p+h の寄与）")
ax.axhline(0, color="k", lw=0.6)
for z in zeros[zeros < 45]:
    ax.axvline(z, color="k", lw=0.6, ls=":")
ax.set_xlabel("周波数 t"); ax.set_ylabel("干渉項")
ax.set_title("色を塗ると消えるもの＝この干渉項。零点のもつれ（素数対のリズム）はここに住んでいる", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/white_light.png")
