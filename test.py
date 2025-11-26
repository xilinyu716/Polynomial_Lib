from polylib import Polynomial, set_fft_threshold

p1 = Polynomial([1, -3, 2])
p2 = Polynomial([1, -1])
print(p1 + p2)
q, r = p1 / p2
print(q, r)
set_fft_threshold(128)