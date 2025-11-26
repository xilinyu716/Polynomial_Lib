import unittest
from polylib import Polynomial, set_fft_threshold

class TestPolyLibPackage(unittest.TestCase):
    def test_add_mul_div_derivative_eval(self):
        p1 = Polynomial([1, -3, 2])
        p2 = Polynomial([1, -1])
        p_add = p1 + p2
        self.assertEqual(p_add.coeffs, [1, -2, 1])
        p_mul = p1 * p2
        self.assertEqual(p_mul.coeffs, [1, -4, 5, -2])
        q, r = p1 / p2
        self.assertEqual(q.coeffs, [1, -2])
        self.assertTrue(all(abs(c) < 1e-9 for c in r.coeffs))
        der = p1.derivative()
        self.assertEqual(der.coeffs, [2, -3])
        val = p1.evaluate(1.0)
        self.assertAlmostEqual(val, 0.0, places=9)

    def test_gcd_roots_bezout(self):
        p1 = Polynomial([1, -3, 2])
        p2 = Polynomial([1, -1])
        g = p1.gcd(p2)
        self.assertEqual(g.coeffs, [1, -1])
        nroots = p1.num_real_roots(-10.0, 10.0)
        self.assertEqual(nroots, 2)
        gcd, s, t = p1.bezout_identity(p2)
        left = s * p1 + t * p2
        self.assertEqual(left.coeffs, gcd.coeffs)
        roots = p1.find_real_roots()
        self.assertTrue(any(abs(r - 1.0) < 1e-6 for r in roots))
        self.assertTrue(any(abs(r - 2.0) < 1e-6 for r in roots))

    def test_fft_threshold_config(self):
        p1 = Polynomial([1] + [0]*200)
        p2 = Polynomial([1] + [0]*200)
        set_fft_threshold(0)
        r1 = (p1 * p2).coeffs
        set_fft_threshold(128)
        r2 = (p1 * p2).coeffs
        self.assertEqual(r1, r2)

if __name__ == '__main__':
    unittest.main()
