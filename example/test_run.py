from src import Polynomial

f = Polynomial([1, 0, -2])   # f(x) = x^2 - 2
g = Polynomial([1, -1])      # g(x) = x - 1

print("f(x) =", f)
print("g(x) =", g)

q, r = f / g
print("f(x) / g(x) = q(x) =", q, ", r(x) =", r)

gcd = f.gcd(g)
print("gcd(f, g) =", gcd)

roots = f.find_real_roots()
print("Approx real roots of f(x):", roots)