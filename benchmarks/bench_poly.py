import random
import time
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import Polynomial
from src.polynomial import set_fft_threshold

def rand_poly(n):
    return Polynomial([random.uniform(-1, 1) for _ in range(n)])

def time_op(op, reps=50):
    t0 = time.perf_counter()
    for _ in range(reps):
        op()
    return time.perf_counter() - t0

def main():
    sizes = [32, 64, 128, 256, 512, 1024]
    for n in sizes:
        p1 = rand_poly(n)
        p2 = rand_poly(n)
        # Measure with different thresholds
        for thr in [0, 32, 64, 128]:
            set_fft_threshold(thr)
            t_mul = time_op(lambda: p1 * p2)
            print(f"n={n} thr={thr} mul {t_mul:.6f}s")
    p3 = rand_poly(300)
    p4 = rand_poly(300)
    t_gcd = time_op(lambda: p3.gcd(p4))
    print(f"n=300 gcd {t_gcd:.6f}s")
    p5 = rand_poly(200)
    t_sc = time_op(lambda: p5.num_real_roots(-10.0, 10.0))
    print(f"n=200 sturm roots {t_sc:.6f}s")

if __name__ == '__main__':
    main()
