from sage.all import *

q  = 12289
n  = 8
k  = 16
ell = 8
dim = k + ell

Zq = Integers(q)
MS = MatrixSpace(Zq, n, dim)
VS = VectorSpace(Zq, dim)

def main():
    A = MS.random_element()
    print("A:")
    print(A)

if __name__ == "__main__":
    main()
