"""外部提案（notes/external_proposal_2026-11-15.md）の二つの主張を数値で確認する"""
import numpy as np
P = 10**6
sieve = np.ones(P+1, bool); sieve[:2] = False
for i in range(2, int(P**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0]
print("もつれ状態 |Ψ(s)> = Σ n^{-s/2}|n>|n> のノルム²（対数）。Im s に依らず Re s のみで決まる")
print("   P      σ=0.3     σ=0.5     σ=0.7     σ=1.0     σ=1.2")
for Pc in [10, 100, 1000, 10**4, 10**5, 10**6]:
    ps = primes[primes <= Pc]
    print(f"{Pc:>7} " + " ".join(f"{np.sum(-np.log1p(-ps**(-sg))):9.2f}" for sg in [0.3, 0.5, 0.7, 1.0, 1.2]))
N = 30
a = np.diag(np.sqrt(np.arange(1, N)), 1); I = np.eye(N)
A = np.kron(a, I); B = np.kron(I, a)
for s in [0.5+14.13j, 0.3+14.13j, 0.9+5j, 2.0+0j]:
    th = np.arctanh(2**(-s/2))
    c = np.cosh(th)*A - np.sinh(th)*B.conj().T
    H = c.conj().T @ c
    print(f"s = {s}:  c†c の最小固有値 = {np.linalg.eigvalsh((H+H.conj().T)/2).min():+.2e}  （A†A ≥ 0 は s に依らない恒等式）")
