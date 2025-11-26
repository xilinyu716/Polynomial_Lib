#include "../polynomial.h"
#include <cassert>
#include <vector>
#include <cmath>
int main() {
    double p1[] = {1.0, -3.0, 2.0};
    int l1 = 3;
    double p2[] = {1.0, -1.0};
    int l2 = 2;
    int add_len = std::max(l1, l2);
    std::vector<double> add_out(add_len);
    poly_add(p1, l1, p2, l2, add_out.data());
    assert(add_out.size() == 3);
    assert(std::fabs(add_out[0] - 1.0) < 1e-9);
    assert(std::fabs(add_out[1] + 2.0) < 1e-9);
    assert(std::fabs(add_out[2] - 1.0) < 1e-9);
    int mul_len = l1 + l2 - 1;
    std::vector<double> mul_out(mul_len);
    poly_mul(p1, l1, p2, l2, mul_out.data());
    assert(mul_out.size() == 4);
    assert(std::fabs(mul_out[0] - 1.0) < 1e-9);
    assert(std::fabs(mul_out[1] + 4.0) < 1e-9);
    assert(std::fabs(mul_out[2] - 5.0) < 1e-9);
    assert(std::fabs(mul_out[3] + 2.0) < 1e-9);
    std::vector<double> mul_fft_out(mul_len);
    poly_mul_fft(p1, l1, p2, l2, mul_fft_out.data());
    for(int i=0;i<mul_len;i++) assert(std::fabs(mul_fft_out[i] - mul_out[i]) < 1e-9);
    std::vector<double> q_out(std::max(l1 - l2 + 1, 1));
    std::vector<double> r_out(l1);
    int len_q = 0, len_r = 0;
    int rc = poly_div(p1, l1, p2, l2, q_out.data(), &len_q, r_out.data(), &len_r);
    assert(rc == 0);
    assert(len_q == 2);
    assert(std::fabs(q_out[0] - 1.0) < 1e-9);
    assert(std::fabs(q_out[1] + 2.0) < 1e-9);
    assert(len_r == 1);
    assert(std::fabs(r_out[0]) < 1e-9);
    std::vector<double> der_out(std::max(l1 - 1, 1));
    poly_derivative(p1, l1, der_out.data());
    assert(der_out.size() == 2);
    assert(std::fabs(der_out[0] - 2.0) < 1e-9);
    assert(std::fabs(der_out[1] + 3.0) < 1e-9);
    double v = poly_eval(p1, l1, 1.0);
    assert(std::fabs(v) < 1e-9);
    std::vector<double> gcd_out(std::max(l1, l2));
    int len_g=0; int rc_g = poly_gcd(p1, l1, p2, l2, gcd_out.data(), &len_g);
    assert(rc_g == 0);
    assert(len_g == 2);
    assert(std::fabs(gcd_out[0] - 1.0) < 1e-9);
    assert(std::fabs(gcd_out[1] + 1.0) < 1e-9);
    int changes_left = poly_sturm_sign_changes(p1, l1, -10.0);
    int changes_right = poly_sturm_sign_changes(p1, l1, 10.0);
    int num_roots = poly_num_real_roots_interval(p1, l1, -10.0, 10.0);
    assert(changes_left - changes_right == 2);
    assert(num_roots == 2);
    return 0;
}
