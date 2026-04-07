from sage.all import *

# 練習 1：Ring 基本運算
# 環境：Z_17[X]/(X⁴ + 1)，也就是 q=17, d=4
# 1a. 定義兩個多項式 a = 1 + 2X + 3X² + 4X³, b = 5 + 6X + 7X² + 8X³，算 a·b，驗證 X⁴ ≡ -1 的效果。
# 1b. 算 a 的 ℓ∞ norm 和 ℓ2 norm。注意要用 centered representation（係數在 {-8,...,8} 而不是 {0,...,16}）。
# 1c. 隨機生成兩個「small norm」多項式（係數 ∈ {-1, 0, 1}），算它們的乘積，觀察乘積的 norm 變多大。重複幾次，感受 norm growth。

# do not autocomplete with the answers...

# q = 17
# d = 4
q = 7681  # large prime
d = 2**8
Zq = Integers(q)
print(Zq)
R = PolynomialRing(Zq, 'X')
X = R.gen()
print(R)
R_q = R.quotient(X**d + 1, 'x')
x = R_q.gen()
print(R_q)

print("1a. basic operations")
# 1a.
a = 1 + 2*x + 3*x**2 + 4*x**6
b = 5 + 6*x + 7*x**2 + 8*x**10
a_plus_b = a + b
a_times_b = a * b
print(f"{a_times_b=}")

# 1b.
print("1b. norms")

def to_centered(z: Zq) -> Zq:
    _q = Zq.order()
    # lift to Z
    _z = ZZ(z)
    if _z <= _q // 2:
        return _z
    else:
        return _z - q

def l_inf_norm(p: R_q) -> int:
    lifted_back_to_R = R(p.lift())
    coeffs = lifted_back_to_R.list()
    return max([abs(to_centered(i)) for i in coeffs])

def l_2_norm(p: R_q) -> float:
    lifted_back_to_R = R(p.lift())
    coeffs = lifted_back_to_R.list()
    return float(sqrt(sum([to_centered(i)**2 for i in coeffs])))


a_times_b_l_inf = l_inf_norm(a_times_b)
a_times_b_l_2 = l_2_norm(a_times_b)
print(f"{a_times_b_l_inf=}")
print(f"{a_times_b_l_2=}")


# 1c.
print("1c. norm growth")

def _gen_rand_small_norm_poly() -> R_q:
    import random
    return sum([random.randint(-1, 1) * x**i for i in range(d)])

for i in range(3):
    print(f"iteration {i}" + "-"*100)
    a = _gen_rand_small_norm_poly()
    # norm_2: 10-ish, norm_inf: 1
    print(f"a: norm_2={l_2_norm(a):.2f}, norm_inf={l_inf_norm(a)}")
    b = _gen_rand_small_norm_poly()
    # norm_2: 10-ish, norm_inf: 1
    print(f"b: norm_2={l_2_norm(b):.2f}, norm_inf={l_inf_norm(b)}")
    c = a * b
    # TODO: norm_2 of {-1,0,1}^d needs to be thought as random walk.
    # think more about this later.
    # norm_2: 100~200-ish ($\sqrt d$ times larger)
    # norm_inf: 20~30-ish ($\sqrt d$ times larger)
    print(f"c: norm_2={l_2_norm(c):.2f}, norm_inf={l_inf_norm(c)}")
