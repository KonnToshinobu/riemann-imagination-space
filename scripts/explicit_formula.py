"""
ヴェイユの明示公式の数値検証
  Σ_ρ h(γ)  =  h(i/2)+h(-i/2) − g(0) log π + (1/2π)∫ h(r) Re ψ(1/4+ir/2) dr
               − 2 Σ_n Λ(n) n^{-1/2} g(log n)
  h(r) = exp(-a r^2) cos(T r),   g(u) = (1/2π)∫ h(r) e^{-iur} dr  (hのフーリエ変換)
左辺 = ゼータ零点だけ、右辺の最後の和 = 素数だけ。
"""
import numpy as np, mpmath as mp, json, os, sys
from math import log, pi, sqrt, exp, cosh

mp.mp.dps = 15
NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 300
cache = f"data/zeros_{NZ}.json"
if os.path.exists(cache):
    zeros = json.load(open(cache))
else:
    zeros = [float(mp.zetazero(n).imag) for n in range(1, NZ + 1)]
    json.dump(zeros, open(cache, "w"))
zeros = np.array(zeros)

# von Mangoldt Λ(n) for n <= NMAX
NMAX = 20000
lam = np.zeros(NMAX + 1)
sieve = np.ones(NMAX + 1, bool); sieve[:2] = False
for i in range(2, int(NMAX**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0]
for p in primes:
    pk = p
    while pk <= NMAX:
        lam[pk] = log(p); pk *= p
ns = np.arange(1, NMAX + 1)

a = 0.0002

def G(u):            # Fourier transform of exp(-a r^2) with 1/2π
    return exp(-u*u/(4*a)) / (2*sqrt(pi*a))

def lhs(T):          # zero side
    # γ>0 の零点のみ保持しているので、共役零点 (γ<0) の分を2倍で補う
    return 2*float(np.sum(np.exp(-a*zeros**2) * np.cos(T*zeros)))

from scipy.special import digamma as _dg
_r = np.linspace(-600, 600, 240001)
_w = np.exp(-a*_r**2) * np.real(_dg(0.25 + 0.5j*_r))
def gamma_term(T):
    return float(np.trapezoid(_w*np.cos(T*_r), _r)) / (2*pi)

def rhs_archimedean(T):      # everything on the right except primes
    g0 = 0.5*(G(-T)+G(T))
    return 2*exp(a/4)*cosh(T/2) - g0*log(pi) + gamma_term(T)

def prime_sum(T):
    gl = 0.5*(np.array([G(u-T)+G(u+T) for u in np.log(ns)]))
    return float(-2*np.sum(lam[1:]/np.sqrt(ns) * gl))

if __name__ == "__main__":
    print(f"零点 {NZ} 個 (γ_max={zeros[-1]:.1f}),  Λ(n) は n≤{NMAX},  a={a}")
    print(f"{'T':>10} {'零点側(左辺)':>14} {'右辺':>14} {'差':>12}   T=log n")
    for n in [2, 3, 4, 5, 7, 11, 13, 17, 29, 31, 97, 101, 127, 128]:
        T = log(n)
        L = lhs(T); R = rhs_archimedean(T) + prime_sum(T)
        print(f"{T:10.5f} {L:14.6f} {R:14.6f} {L-R:12.2e}   n={n}")
    # 素数でない T
    for T in [0.5, 1.0, 2.0, 3.3, 4.0]:
        L = lhs(T); R = rhs_archimedean(T) + prime_sum(T)
        print(f"{T:10.5f} {L:14.6f} {R:14.6f} {L-R:12.2e}   (素数のlogではない)")
