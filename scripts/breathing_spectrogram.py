"""
素数の呼吸のスペクトログラム
  ψ(x) = x − Σ_ρ x^ρ/ρ − log 2π − ½log(1−x^{−2})   （フォン・マンゴルトの明示公式）
  時間を u = log x で測ると、零点 ρ=1/2+iγ の倍音は  −2 Re[e^{iγu}/ρ]・√x
  → r(u) = (ψ(x) − x + log 2π)/√x  は「基音 x の周りの呼吸」。
     RH ⇔ すべての倍音の包絡線が時間に沿って一定（横縞の明るさが変わらない）
  偽の零点 β+iγ_f (β>1/2) を混ぜると、その縞だけ e^{(β−1/2)u} で明るくなる。
"""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

N = 10_000_000
# Λ(n) と ψ(x)
sieve = np.ones(N + 1, bool); sieve[:2] = False
for i in range(2, int(N**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0]
lam = np.zeros(N + 1); lam[primes] = np.log(primes)
for p in primes[primes <= int(N**0.5)]:
    q = p * p
    while q <= N: lam[q] = np.log(p); q *= p
psi = np.cumsum(lam)
del sieve, lam

# 呼吸 r(u)
u0, u1 = np.log(30.0), np.log(N)
M = 2**17
u = np.linspace(u0, u1, M); du = u[1] - u[0]
x = np.exp(u)
r = (psi[np.floor(x).astype(np.int64)] - x + np.log(2*np.pi)) / np.sqrt(x)

# 偽の零点
beta, gf = 0.7, 54.7                       # 52.97 と 56.45 の間
rho_f = beta + 1j*gf
r_fake = r - 2*np.real(np.exp(rho_f*u) / rho_f) / np.sqrt(x)

# スペクトログラム（ガウス窓 STFT, 零詰め FFT）
sig = 1.3
centers = np.linspace(u0 + 2*sig, u1 - 2*sig, 90)
gmax = 100.0
pad = 16
def stft(y):
    out = []
    for c in centers:
        w = np.exp(-(u - c)**2 / (2*sig**2))
        Y = np.fft.rfft(y * w, n=pad*M) * du
        f = np.fft.rfftfreq(pad*M, d=du) * 2*np.pi
        k = f <= gmax
        out.append(np.abs(Y[k]) / (sig*np.sqrt(2*np.pi)/2))   # 振幅 A の cos が A と表示される正規化
    return f[k], np.array(out).T                              # (freq, time)
f, S = stft(r)
_, Sf = stft(r_fake)
# 倍音 1/|ρ| の自然な減衰を補正（γ/2 倍）: 本物の零点の縞は ≈1 の明るさに揃うはず
wgt = np.maximum(f, 1.0)[:, None] / 2
Sn, Sfn = S*wgt, Sf*wgt

zeros = np.array(json.load(open("data/zeros_300.json")))
zeros = zeros[zeros < gmax]
print(f"ψ(x) を x ≤ {N:,} で計算。呼吸 r(u) の標準偏差 = {r.std():.3f}, 最大 |r| = {np.abs(r).max():.3f}")
# 本物の縞の明るさの時間変化（最初の 10 個の零点で平均）と偽の縞
def stripe(Sm, g):
    i = np.argmin(np.abs(f - g)); return Sm[i]
real_tracks = np.array([stripe(Sn, g) for g in zeros[:10]])
fake_track = stripe(Sfn, gf)
print("本物の零点の縞（10 本平均）: 前半", real_tracks[:, :45].mean().round(3), " 後半", real_tracks[:, 45:].mean().round(3))
print(f"偽の零点 β={beta} の縞:       前半", fake_track[:45].mean().round(3), " 後半", fake_track[45:].mean().round(3),
      f"  （予言：e^{{(β−1/2)Δu}} = {np.exp((beta-0.5)*(centers[-1]-centers[0])):.1f} 倍）")

# ---------- 図 ----------
fig = plt.figure(figsize=(15, 12.5), dpi=150, constrained_layout=True)
gs = fig.add_gridspec(3, 2, height_ratios=[1, 2.1, 2.1], width_ratios=[3.2, 1])
xt = [10**k for k in range(2, 8)]

ax = fig.add_subplot(gs[0, :])
ax.plot(u, r, color="#1f5fbf", lw=0.4)
ax.axhline(0, color="#d9541e", lw=1.5, label="基音 x（ワームホールを通る唯一の呼吸）")
ax.set_xlim(u0, u1); ax.set_xticks(np.log(xt)); ax.set_xticklabels([f"$10^{k}$" for k in range(2, 8)])
ax.set_ylabel("(ψ(x) − x)/√x")
ax.set_title(f"素数の呼吸： ψ(x) が基音 x の周りで揺れる様子（√x で割って時間 log x で見る）、x ≤ {N:,}", fontsize=11, loc="left")
ax.legend(fontsize=9, loc="upper right"); ax.grid(alpha=0.25)

def spec(ax, Sm, title, zmarks=True, fake=None):
    im = ax.pcolormesh(centers, f, Sm, cmap="magma", vmin=0, vmax=1.6, shading="auto", rasterized=True)
    if zmarks:
        for g in zeros: ax.plot([centers[-1] - 0.25, centers[-1] - 0.02], [g, g], color="#7ee0ff", lw=0.9)
        ax.text(centers[-1] - 0.3, 97, "←零点 γ", color="#7ee0ff", fontsize=8, ha="right", va="top")
    if fake is not None:
        ax.annotate(f"偽の零点 β={beta}, γ={fake}", xy=(centers[-1], fake), xytext=(centers[len(centers)//2], fake + 20),
                    color="white", fontsize=9, arrowprops=dict(arrowstyle="->", color="white"))
    ax.set_xlim(centers[0], centers[-1]); ax.set_ylim(0, gmax)
    kk = [k for k in range(2, 8) if centers[0] <= np.log(10**k) <= centers[-1]]
    ax.set_xticks([np.log(10**k) for k in kk]); ax.set_xticklabels([f"$10^{k}$" for k in kk])
    ax.set_ylabel("周波数 γ"); ax.set_title(title, fontsize=10.5, loc="left")
    return im

ax1 = fig.add_subplot(gs[1, 0])
im = spec(ax1, Sn, "本物の呼吸のスペクトログラム：横縞がすべて零点の位置にあり、明るさが時間に沿って変わらない（＝ RH が言っていること）")
ax2 = fig.add_subplot(gs[2, 0])
spec(ax2, Sfn, "偽の零点（線から 0.2 外れた倍音）を一本混ぜた呼吸：その縞だけが時間とともに明るくなり、調和が壊れる", fake=gf)
ax2.set_xlabel("x（時間 u = log x）")
cb = fig.colorbar(im, ax=[ax1, ax2], shrink=0.6, pad=0.01, aspect=40); cb.set_label("倍音の振幅 ×(γ/2)   （本物はどれも ≈1）")

ax3 = fig.add_subplot(gs[1, 1])
Tavg = Sn.mean(axis=1)
ax3.plot(Tavg, f, color="#1f5fbf", lw=1.2)
for g in zeros: ax3.axhline(g, color="#d9541e", lw=0.5, alpha=0.6)
ax3.set_ylim(0, gmax); ax3.set_xlim(0, 1.6)
ax3.set_xlabel("時間平均した振幅"); ax3.set_title("楽譜：ピークは零点（橙線）", fontsize=10.5, loc="left"); ax3.grid(alpha=0.25)

ax4 = fig.add_subplot(gs[2, 1])
for k, g in enumerate(zeros[:10]):
    ax4.plot(centers, real_tracks[k], color="#1f5fbf", lw=1, alpha=0.6, label="本物の零点の縞（10 本）" if k == 0 else None)
ax4.plot(centers, fake_track, color="#d9541e", lw=2.2, label=f"偽の零点の縞 (β={beta})")
ax4.plot(centers, fake_track[0]*np.exp((beta-0.5)*(centers-centers[0])), "k--", lw=1, label="予言 e^{(β−1/2)u}")
ax4.set_xticks([np.log(10**k) for k in range(3, 7)]); ax4.set_xticklabels([f"$10^{k}$" for k in range(3, 7)])
ax4.set_xlabel("x"); ax4.set_ylabel("縞の明るさ"); ax4.set_title("一本の倍音の呼吸：平ら vs 膨らむ", fontsize=10.5, loc="left")
ax4.legend(fontsize=7.5); ax4.grid(alpha=0.25)

fig.suptitle("素数の呼吸を録音して、倍音ごとに分解する ― リーマン予想 ＝ どの倍音も永遠に基音と同じ割合で呼吸する", fontsize=12.5)
fig.savefig("figures/breathing_spectrogram.png")
