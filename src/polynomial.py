"""
Polynomial Library
This library provides operations on polynomials including:
- GCD using Euclidean algorithm
- Sturm sequence
- Root finding
- Bezout's identity
- Polynomial addition and multiplication
"""

import math
import numpy as np
from typing import List, Tuple, Union

class Polynomial:
    def __init__(self, coefficients: List[float]):
        """Initialize a polynomial with given coefficients.
        
        Args:
            coefficients: List of coefficients from highest degree to constant term.
                         For example, [1, 2, 3] represents x² + 2x + 3.
        """
        # Remove leading zeros
        while len(coefficients) > 1 and math.isclose(coefficients[0], 0, abs_tol=1e-10):
            coefficients = coefficients[1:]
        self.coeffs = coefficients
    
    def __repr__(self) -> str:
        terms = []
        for power, coeff in enumerate(reversed(self.coeffs)):
            if math.isclose(coeff, 0, abs_tol=1e-10):
                continue
                
            power = len(self.coeffs) - 1 - power
            
            if power == 0:
                term = f"{coeff:.4f}".rstrip('0').rstrip('.') if '.' in f"{coeff:.4f}" else f"{coeff}"
            else:
                if math.isclose(coeff, 1, abs_tol=1e-10):
                    coeff_str = ""
                elif math.isclose(coeff, -1, abs_tol=1e-10):
                    coeff_str = "-"
                else:
                    coeff_str = f"{coeff:.4f}".rstrip('0').rstrip('.') if '.' in f"{coeff:.4f}" else f"{coeff}"
                
                if power == 1:
                    term = f"{coeff_str}x"
                else:
                    term = f"{coeff_str}x^{power}"
            
            terms.append(term)
        
        if not terms:
            return "0"
        
        # Combine terms with proper signs
        expression = terms[0]
        for term in terms[1:]:
            if term.startswith('-'):
                expression += f" - {term[1:]}"
            else:
                expression += f" + {term}"
        
        return expression
    
    def degree(self) -> int:
        """Return the degree of the polynomial."""
        return len(self.coeffs) - 1
    
    def evaluate(self, x: float) -> float:
        """Evaluate the polynomial at point x using Horner's method."""
        result = 0.0
        for coeff in self.coeffs:
            result = result * x + coeff
        return result
    
    def derivative(self) -> 'Polynomial':
        """Return the derivative of the polynomial."""
        if len(self.coeffs) == 1:
            return Polynomial([0])
        
        new_coeffs = []
        for i, coeff in enumerate(self.coeffs[:-1]):
            new_coeffs.append(coeff * (len(self.coeffs) - 1 - i))
        return Polynomial(new_coeffs)
    
    def __add__(self, other: 'Polynomial') -> 'Polynomial':
        """Add two polynomials."""
        max_degree = max(self.degree(), other.degree())
        new_coeffs = [0] * (max_degree + 1)
        
        for i in range(len(self.coeffs)):
            new_coeffs[max_degree - self.degree() + i] += self.coeffs[i]
        
        for i in range(len(other.coeffs)):
            new_coeffs[max_degree - other.degree() + i] += other.coeffs[i]
        
        return Polynomial(new_coeffs)
    
    def __mul__(self, other: 'Polynomial') -> 'Polynomial':
        """Multiply two polynomials using convolution."""
        degree = self.degree() + other.degree()
        new_coeffs = [0] * (degree + 1)
        
        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(other.coeffs):
                new_coeffs[i + j] += a * b
        
        return Polynomial(new_coeffs)
    
    def __sub__(self, other: 'Polynomial') -> 'Polynomial':
        """Subtract two polynomials."""
        return self + (other * Polynomial([-1]))
    
    def __truediv__(self, other: 'Polynomial') -> Tuple['Polynomial', 'Polynomial']:
        """Polynomial division returning quotient and remainder."""
        if other.degree() == 0 and math.isclose(other.coeffs[0], 0, abs_tol=1e-10):
            raise ZeroDivisionError("Cannot divide by zero polynomial")
            
        if self.degree() < other.degree():
            return Polynomial([0]), self
        
        remainder = Polynomial(self.coeffs.copy())
        divisor = other
        quotient_coeffs = [0] * (self.degree() - divisor.degree() + 1)
        
        while remainder.degree() >= divisor.degree():
            # Compute the next term of the quotient
            if remainder.is_zero():
                break
            leading_coeff = remainder.coeffs[0] / divisor.coeffs[0]
            power = remainder.degree() - divisor.degree()
            
            # Create the term polynomial
            term_coeffs = [0] * (power + 1)
            term_coeffs[0] = leading_coeff
            term = Polynomial(term_coeffs)
            
            # Update quotient and remainder
            quotient_coeffs[power] = leading_coeff
            remainder = remainder - (term * divisor)
        
        quotient = Polynomial(quotient_coeffs)
        
        # Remove any minuscule coefficients due to floating point errors
        cleaned_remainder = []
        for coeff in remainder.coeffs:
            if abs(coeff) < 1e-10:
                cleaned_remainder.append(0.0)
            else:
                cleaned_remainder.append(coeff)
        
        return quotient, Polynomial(cleaned_remainder)
    
    def gcd(self, other: 'Polynomial') -> 'Polynomial':
        """Compute GCD of two polynomials using Euclidean algorithm."""
        a = self
        b = other

        while not b.is_zero():
            _, remainder = a / b
            a = b
            b = remainder
            print(f"a: {a}, b: {b}")
            print(f"degree a: {a.degree()}, degree b: {b.degree()}")

        # Make the GCD monic
        if not a.is_zero():
            leading_coeff = a.coeffs[0]
            monic_coeffs = [c / leading_coeff for c in a.coeffs]
            return Polynomial(monic_coeffs)
        return a
    
    def sturm_sequence(self) -> List['Polynomial']:
        """Generate the Sturm sequence for the polynomial."""
        sequence = []
        sequence.append(self)
        sequence.append(self.derivative())
        
        i = 1
        while sequence[i].degree() > 0:
            _, remainder = sequence[i-1] / sequence[i]
            # Negate the remainder according to Sturm's theorem
            next_poly = remainder * Polynomial([-1])
            sequence.append(next_poly)
            i += 1
        
        return sequence
    
    def sign_changes(self, x: float) -> int:
        """Count sign changes in the Sturm sequence evaluated at x."""
        sequence = self.sturm_sequence()
        values = [poly.evaluate(x) for poly in sequence]
        
        # Remove zeros by treating them as positive (but this can be improved)
        cleaned_values = []
        for v in values:
            if math.isclose(v, 0, abs_tol=1e-10):
                cleaned_values.append(0.0)
            else:
                cleaned_values.append(v)
        
        sign_changes = 0
        prev_sign = 0
        
        for v in cleaned_values:
            if v > 0:
                current_sign = 1
            elif v < 0:
                current_sign = -1
            else:
                current_sign = 0
            
            if (prev_sign != 0) and (current_sign != 0) and (current_sign != prev_sign):
                sign_changes += 1
            
            if current_sign != 0:
                prev_sign = current_sign
        
        return sign_changes
    
    def num_real_roots(self, a: float, b: float) -> int:
        """Count the number of distinct real roots in interval [a, b] using Sturm's theorem."""
        if a > b:
            a, b = b, a
        
        return self.sign_changes(a) - self.sign_changes(b)
    
    def find_real_roots(self, tolerance: float = 1e-6, max_iter: int = 100) -> List[float]:
        """Find all real roots of the polynomial using Sturm sequence and binary search."""
        if self.degree() == 0:
            return []
        
        # Find a bound for the roots
        max_coeff = max(abs(c) for c in self.coeffs[1:])
        bound = 1 + max_coeff / abs(self.coeffs[0])
        search_interval = [-bound, bound]
        
        # Find intervals with roots
        def find_root_intervals(poly, a, b, intervals):
            num_roots = poly.num_real_roots(a, b)
            
            if num_roots == 0:
                return
            
            if num_roots == 1:
                intervals.append((a, b))
                return
            
            # Multiple roots, subdivide the interval
            mid = (a + b) / 2
            find_root_intervals(poly, a, mid, intervals)
            find_root_intervals(poly, mid, b, intervals)
        
        intervals = []
        find_root_intervals(self, search_interval[0], search_interval[1], intervals)
        
        # Refine roots using binary search
        roots = []
        for a, b in intervals:
            if math.isclose(self.evaluate(a), 0, abs_tol=tolerance):
                roots.append(a)
                continue
            
            if math.isclose(self.evaluate(b), 0, abs_tol=tolerance):
                roots.append(b)
                continue
            
            left = a
            right = b
            
            for _ in range(max_iter):
                mid = (left + right) / 2
                val = self.evaluate(mid)
                
                if math.isclose(val, 0, abs_tol=tolerance):
                    roots.append(mid)
                    break
                
                if self.evaluate(left) * val < 0:
                    right = mid
                else:
                    left = mid
            else:
                roots.append(mid)
        
        # Remove duplicates (roots that are very close)
        unique_roots = []
        for root in sorted(roots):
            if not unique_roots or abs(root - unique_roots[-1]) > tolerance:
                unique_roots.append(root)
        
        return unique_roots
    
    def is_zero(self) -> bool:
        return all(math.isclose(c, 0.0, abs_tol=1e-10) for c in self.coeffs)
    
    def bezout_identity(self, other: 'Polynomial') -> Tuple['Polynomial', 'Polynomial', 'Polynomial']:
        """Find polynomials s and t such that s*self + t*other = gcd(self, other).
        
        Returns:
            A tuple (gcd, s, t) where gcd is the GCD of self and other,
            and s, t are the coefficients in Bezout's identity.
        """
        old_r, r = self, other
        old_s, s = Polynomial([1]), Polynomial([0])
        old_t, t = Polynomial([0]), Polynomial([1])
        
        while r.degree() >= 0 and not math.isclose(r.coeffs[0], 0, abs_tol=1e-10):
            quotient, remainder = old_r / r
            
            old_r, r = r, remainder
            old_s, s = s, old_s - (quotient * s)
            old_t, t = t, old_t - (quotient * t)
        
        # Make the GCD monic
        if old_r.degree() >= 0 and not math.isclose(old_r.coeffs[0], 0, abs_tol=1e-10):
            leading_coeff = old_r.coeffs[0]
            monic_gcd = Polynomial([c / leading_coeff for c in old_r.coeffs])
            monic_s = Polynomial([c / leading_coeff for c in old_s.coeffs])
            monic_t = Polynomial([c / leading_coeff for c in old_t.coeffs])
            return (monic_gcd, monic_s, monic_t)
        
        return (old_r, old_s, old_t)


# Example usage
if __name__ == "__main__":
    # Create some polynomials
    p1 = Polynomial([1, -3, 2])  # x^2 - 3x + 2
    p2 = Polynomial([1, -1])     # x - 1
    
    print(f"p1 = {p1}")
    print(f"p2 = {p2}")
    
    # Test addition
    p_add = p1 + p2
    print(f"\np1 + p2 = {p_add}")
    
    # Test multiplication
    p_mul = p1 * p2
    print(f"p1 * p2 = {p_mul}")
    
    # Test division
    quotient, remainder = p1 / p2
    print(f"\np1 / p2: quotient = {quotient}, remainder = {remainder}")
    
    # Test GCD
    p_gcd = p1.gcd(p2)
    print(f"\nGCD(p1, p2) = {p_gcd}")
    
    # Test Bezout's identity
    gcd, s, t = p1.bezout_identity(p2)
    print(f"\nBezout's identity for p1 and p2:")
    print(f"s = {s}")
    print(f"t = {t}")
    print(f"Verification: s*p1 + t*p2 = {s*p1 + t*p2}")
    print(f"GCD: {gcd}")
    
    # Test Sturm sequence
    print("\nSturm sequence for p1 = x^2 - 3x + 2:")
    sturm = p1.sturm_sequence()
    for i, poly in enumerate(sturm):
        print(f"S_{i} = {poly}")
    
    # Test root finding
    print("\nRoots of p1:")
    roots = p1.find_real_roots()
    for root in roots:
        print(f"Root: {root:.6f}, p1(root) = {p1.evaluate(root):.2e}")
    
    # Test with a cubic polynomial
    p3 = Polynomial([1, -6, 11, -6])  # x^3 - 6x^2 + 11x - 6
    print("\nTesting with p3 = x^3 - 6x^2 + 11x - 6")
    print("Sturm sequence:")
    sturm3 = p3.sturm_sequence()
    for i, poly in enumerate(sturm3):
        print(f"S_{i} = {poly}")
    
    print("\nRoots of p3:")
    roots3 = p3.find_real_roots()
    for root in roots3:
        print(f"Root: {root:.6f}, p3(root) = {p3.evaluate(root):.2e}")