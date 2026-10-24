import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import explicit_formula as ef
from math import log

# 零点が「見ている」素数軌道:  D(T) = Σ_ρ h_T(γ) − (アルキメデス項)  ← 零点のみ
# 素数側:                     P(T) = −2 Σ_n Λ(n) n^{-1/2} g_T(log n)       ← 素数のみ
Ts = np.linspace(0.3, 5.0, 6000)
D = np.array([ef.lhs(T) - ef.rhs_archimedean(T) for T in Ts])
P = np.array([ef.prime_sum(T) for T in Ts])

jp = [f.name for f in font_manager.fontManager.ttflist if "Noto" in f.name and ("CJK" in f.name or "JP" in f.name)]
if jp: plt.rcParams["font.family"] = jp[0]

fig, ax = plt.subplots(figsize=(13, 5.2), dpi=150)
ax.plot(Ts, D, color="#1f5fbf", lw=2.2, label="零点だけから計算した波  Σρ h(γ) − (ガンマ項)")
ax.plot(Ts, P, color="#d9541e", lw=1.0, ls="--", label="素数だけから計算した和  −2 Σ Λ(n) n^{-1/2} g(log n)")
for n in range(2, 200):
    pass
for n in list(range(2,40))+[41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181]:
    if ef.lam[n] > 0:
        T = log(n); p = int(round(np.exp(ef.lam[n])))
        k = round(T/ef.lam[n])
        lbl = f"{p}" if k == 1 else f"{p}^{k}"
        ax.annotate(lbl,(T,P[np.argmin(np.abs(Ts-T))]),textcoords="offset points",xytext=(0,-9),ha="center",fontsize=6.5 if n<60 else 5.6,color="#333",rotation=0 if n<60 else 90)
        ax.axvline(T, color="#888", lw=0.4, alpha=0.5)
    if False and n < 60 or (n < 200 and n % 2 == 1 and k == 1 and n % 3):
            ax.text(T, 0.3, lbl, ha="center", va="top", fontsize=6.5, rotation=90, color="#444")
ax.set_xlabel("T  （= 閉軌道の長さ log n）")
ax.set_ylabel("振幅")
ax.set_title(f"ヴェイユの明示公式：零点 {len(ef.zeros)} 個が素数の閉軌道（長さ log p^k、重み Λ(n)/√n）を再構成する", fontsize=11)
lo=min(D.min(),P.min())*1.08; ax.set_ylim(lo, 0.35)
ax.legend(loc="upper right", fontsize=8.5, framealpha=0.95)
ax.grid(alpha=0.2)
fig.tight_layout()
fig.savefig("figures/explicit_formula_orbits.png")
print("max |D−P| =", np.max(np.abs(D-P)))
