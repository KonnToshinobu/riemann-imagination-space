"""色つき螺旋（実ディリクレ指標 mod q）の L 関数の零点を、臨界線上の Z 関数の符号変化で求める"""
import mpmath as mp, numpy as np, json, sys, time
mp.mp.dps = 12
q = int(sys.argv[1]); T = float(sys.argv[2])
chars = {3: [0, 1, -1], 4: [0, 1, 0, -1], 5: [0, 1, -1, -1, 1], 7: [0, 1, 1, -1, 1, -1, -1]}
chi = chars[q]
a = 0 if chi[q-1] == 1 else 1            # χ(−1)=+1 → 偶 (a=0), −1 → 奇 (a=1)
def Z(t):
    s = mp.mpc(0.5, t)
    v = (mp.mpf(q)/mp.pi)**(s/2) * mp.gamma((s+a)/2) * mp.dirichlet(s, chi)
    return v.real                          # 実指標なら ε=1 で Λ(1/2+it) は実
t0 = time.time()
step = 0.08
ts = np.arange(0.5, T, step)
vals = np.array([float(Z(t)) for t in ts])
zeros = []
for i in range(len(ts)-1):
    if vals[i]*vals[i+1] < 0:
        zeros.append(float(mp.findroot(Z, (ts[i], ts[i+1]), solver="bisect", tol=1e-9)))
json.dump(zeros, open(f"data/Lzeros_q{q}.json", "w"))
# 検算：零点の個数の理論値 N(T) ≈ (T/2π) log(qT/(2πe)) + (定数)
NT = T/(2*np.pi)*np.log(q*T/(2*np.pi*np.e))
print(f"q={q} (a={a}): 零点 {len(zeros)} 個 (T≤{T}),  理論 ≈ {NT:.0f},  {time.time()-t0:.0f}s")
