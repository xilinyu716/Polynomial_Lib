#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <cstddef>
#if defined(_WIN32)
#  ifdef BUILDING_POLY
#    define POLY_API __declspec(dllexport)
#  else
#    define POLY_API __declspec(dllimport)
#  endif
#else
#  define POLY_API
#endif
extern "C" {
POLY_API void poly_add(const double* a, int len_a, const double* b, int len_b, double* out);
POLY_API void poly_sub(const double* a, int len_a, const double* b, int len_b, double* out);
POLY_API void poly_mul(const double* a, int len_a, const double* b, int len_b, double* out);
POLY_API void poly_mul_fft(const double* a, int len_a, const double* b, int len_b, double* out);
POLY_API int poly_div(const double* a, int len_a, const double* b, int len_b, double* q_out, int* len_q, double* r_out, int* len_r);
POLY_API void poly_derivative(const double* a, int len_a, double* out);
POLY_API double poly_eval(const double* a, int len_a, double x);
POLY_API int poly_gcd(const double* a, int len_a, const double* b, int len_b, double* out, int* len_out);
POLY_API int poly_sturm_sign_changes(const double* a, int len_a, double x);
POLY_API int poly_num_real_roots_interval(const double* a, int len_a, double left, double right);

    // Find all real roots
    // roots: pre-allocated buffer of size at least len_a - 1
    // Returns number of roots found
    POLY_API int poly_find_real_roots(const double* a, int len_a, double* roots);

    // Bezout's identity: s*a + t*b = gcd(a, b)
    // gcd, s, t: pre-allocated buffers. 
    // Recommended sizes:
    // gcd: max(len_a, len_b)
    // s: len_b + 1
    // t: len_a + 1
    // len_gcd, len_s, len_t: pointers to store the actual lengths
    POLY_API int poly_bezout(const double* a, int len_a, const double* b, int len_b,
                             double* gcd, int* len_gcd,
                             double* s, int* len_s,
                             double* t, int* len_t);
    POLY_API void poly_set_fft_threshold(int thr);
}
#endif // POLYNOMIAL_H
