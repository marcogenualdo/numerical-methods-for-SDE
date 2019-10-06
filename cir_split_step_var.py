from math import sqrt, exp

a = 1
b = 1
sigma = 2
X0 = 1

dt = 0.05

k_end = int(1 // dt) + 1

d = 4 * a / sigma ** 2
lbb = 4 / (sigma ** 2 * dt) * exp(b * dt)

m1 = lambda k: d * sum([lbb ** i for i in range(k)]) + lbb ** k * X0

alpha = lambda k: d **2 + 2 * d + 2 * lbb * (2 + d) * m1(k-1)

m2 = lambda k: alpha(k) * sum([lbb ** (2 * i) for i in range(k)]) + lbb ** k * X0 ** 2

var = m2(k_end) - m1(k_end) ** 2
print(k_end, var)
