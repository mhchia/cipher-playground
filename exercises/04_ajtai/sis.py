from sage.all import *


def to_centered(z, Fq):
    _q = Fq.order()
    # lift to Z
    _z = ZZ(z)
    if _z <= _q // 2:
        return _z
    else:
        return _z - _q


def l_inf_norm(v, Fq) -> int:
    return max([abs(to_centered(i, Fq)) for i in v])


def check_norm(v, beta, Fq) -> bool:
    return l_inf_norm(v, Fq) <= beta



def find_short_vector(A, beta, Fq, max_trials=50000):
    """
    Iterate through ker(A) and check if there exists z s.t. Az=0 and |z|_\infty <= beta
    """
    q = Fq.order()
    K = A.right_kernel()
    # e.g. z = [-9c-d+1, -7c-8d-6, c, d] where c, d \in F_q
    # then K.dimension()=2, all possibilities are q**2
    num_sols = q**K.dimension()
    num_trials = min(max_trials, num_sols)

    print(f"Finding short vector for {A=}, {num_trials=}, {beta=}")

    for i in range(num_trials):
        z = K.random_element()
        if z.is_zero():
            continue
        if check_norm(z, beta, Fq):
            print(f"Found a solution after {i+1} iterations: {z=}, {beta=}")
            return z
    raise ValueError("Unable to find solutions")


def main():
    q = 12289
    # d = 1024
    Fq = GF(q)

    MS = MatrixSpace(Fq, 2, 4)
    A = MS.random_element()

    # Test: should fail w.h.p due to the low \beta
    try:
        find_short_vector(A, beta=1, Fq=Fq)
    except ValueError:
        print(f"Failed to find low-norm vector: {A=}, {beta=}")
    # Test: should succeed with a larger \beta
    find_short_vector(A, beta=1000, Fq=Fq)


if __name__ == "__main__":
    main()
