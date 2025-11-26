#include "polynomial.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
static size_t next_pow2(size_t x){ size_t n=1; while(n<x) n<<=1; return n; }
static void bit_reverse(std::vector<std::complex<double>>& a){ size_t n=a.size(); size_t j=0; for(size_t i=1;i<n;i++){ size_t bit=n>>1; for(; j & bit; bit>>=1) j^=bit; j^=bit; if(i<j) std::swap(a[i],a[j]); } }
static void fft(std::vector<std::complex<double>>& a, bool inv){ size_t n=a.size(); bit_reverse(a); for(size_t len=2; len<=n; len<<=1){ double ang=(inv?-1.0:1.0)*2.0*M_PI/len; std::complex<double> wlen(std::cos(ang), std::sin(ang)); for(size_t i=0;i<n;i+=len){ std::complex<double> w(1.0,0.0); for(size_t j=0;j<len/2;j++){ auto u=a[i+j]; auto v=a[i+j+len/2]*w; a[i+j]=u+v; a[i+j+len/2]=u-v; w*=wlen; } } } if(inv){ for(size_t i=0;i<n;i++) a[i]/=double(n); } }
static void clamp_eps(std::vector<double>& v){ for(auto &c: v){ if(std::fabs(c)<1e-10) c=0.0; } }
static int degree_len(int len) { return len - 1; }
static void trim(std::vector<double>& v) {
    int i = 0;
    while (i + 1 < (int)v.size() && std::fabs(v[i]) < 1e-12) i++;
    if (i > 0) v.erase(v.begin(), v.begin() + i);
}
static int g_fft_threshold = 64;
extern "C" {
void poly_add(const double* a, int len_a, const double* b, int len_b, double* out) {
    int len = std::max(len_a, len_b);
    int sa = len - len_a;
    int sb = len - len_b;
    for (int i = 0; i < len; ++i) {
        double s = 0.0;
        if (i >= sa) s += a[i - sa];
        if (i >= sb) s += b[i - sb];
        out[i] = s;
    }
}
void poly_sub(const double* a, int len_a, const double* b, int len_b, double* out) {
    int len = std::max(len_a, len_b);
    int sa = len - len_a;
    int sb = len - len_b;
    for (int i = 0; i < len; ++i) {
        double s = 0.0;
        if (i >= sa) s += a[i - sa];
        if (i >= sb) s -= b[i - sb];
        out[i] = s;
    }
}
void poly_mul(const double* a, int len_a, const double* b, int len_b, double* out) {
    int len = len_a + len_b - 1;
    if (len >= g_fft_threshold) { poly_mul_fft(a, len_a, b, len_b, out); return; }
    for (int i = 0; i < len; ++i) out[i] = 0.0;
    for (int i = 0; i < len_a; ++i) {
        for (int j = 0; j < len_b; ++j) {
            out[i + j] += a[i] * b[j];
        }
    }
}
static bool is_zero_vec(const std::vector<double>& v) {
    for (double c : v) if (std::fabs(c) > 1e-12) return false;
    return true;
}
int poly_div(const double* a, int len_a, const double* b, int len_b, double* q_out, int* len_q, double* r_out, int* len_r) {
    if (len_b == 1 && std::fabs(b[0]) < 1e-18) return 1;
    std::vector<double> rem(a, a + len_a);
    std::vector<double> div(b, b + len_b);
    trim(rem);
    trim(div);
    if ((int)rem.size() < (int)div.size()) {
        std::vector<double> q(1, 0.0);
        trim(q);
        *len_q = (int)q.size();
        for (int i = 0; i < *len_q; ++i) q_out[i] = q[i];
        std::vector<double> r = rem;
        if (r.empty()) r.push_back(0.0);
        *len_r = (int)r.size();
        for (int i = 0; i < *len_r; ++i) r_out[i] = r[i];
        return 0;
    }
    std::vector<double> q(degree_len((int)rem.size()) - degree_len((int)div.size()) + 1, 0.0);
    while ((int)rem.size() >= (int)div.size() && !is_zero_vec(rem)) {
        double lc = rem[0] / div[0];
        int p = (int)rem.size() - (int)div.size();
        std::vector<double> term(p + 1, 0.0);
        term[0] = lc;
        std::vector<double> prod(p + (int)div.size(), 0.0);
        for (int i = 0; i < (int)div.size(); ++i) prod[i] = div[i] * lc;
        for (int i = 0; i < (int)prod.size(); ++i) {
            if (i < (int)rem.size()) rem[i] -= prod[i];
        }
        trim(rem);
        int qi = (int)q.size() - 1 - p;
        q[qi] = lc;
    }
    trim(q);
    if (q.empty()) q.push_back(0.0);
    if (rem.empty()) rem.push_back(0.0);
    *len_q = (int)q.size();
    *len_r = (int)rem.size();
    for (int i = 0; i < *len_q; ++i) q_out[i] = q[i];
    for (int i = 0; i < *len_r; ++i) r_out[i] = rem[i];
    return 0;
}
void poly_derivative(const double* a, int len_a, double* out) {
    if (len_a <= 1) { out[0] = 0.0; return; }
    for (int i = 0; i < len_a - 1; ++i) {
        out[i] = a[i] * (len_a - 1 - i);
    }
}
double poly_eval(const double* a, int len_a, double x) {
    double r = 0.0;
    for (int i = 0; i < len_a; ++i) r = r * x + a[i];
    return r;
}
static int sign(double v){ const double eps=1e-12; if(std::fabs(v)<eps) return 0; return v>0?1:-1; }
int poly_gcd(const double* a, int len_a, const double* b, int len_b, double* out, int* len_out) {
    std::vector<double> A(a, a + len_a); trim(A);
    std::vector<double> B(b, b + len_b); trim(B);
    if (is_zero_vec(A) && is_zero_vec(B)) {
        out[0] = 0.0; *len_out = 1; return 0;
    }
    if (is_zero_vec(A)) { if (!is_zero_vec(B)) { double lc = B[0]; for (double &c : B) c /= lc; } if(B.empty()) B.push_back(0.0); *len_out = (int)B.size(); for(int i=0;i<*len_out;++i) out[i]=B[i]; return 0; }
    if (is_zero_vec(B)) { double lc = A[0]; for (double &c : A) c /= lc; if(A.empty()) A.push_back(0.0); *len_out = (int)A.size(); for(int i=0;i<*len_out;++i) out[i]=A[i]; return 0; }
    while (!is_zero_vec(B)) {
        std::vector<double> q(std::max((int)A.size() - (int)B.size() + 1, 1));
        std::vector<double> r(A.size());
        int len_q=0, len_r=0;
        int rc = poly_div(A.data(), (int)A.size(), B.data(), (int)B.size(), q.data(), &len_q, r.data(), &len_r);
        if (rc != 0) { return rc; }
        std::vector<double> R(r.begin(), r.begin() + len_r);
        trim(R);
        A.swap(B);
        B.swap(R);
    }
    trim(A);
    if (!is_zero_vec(A)) { double lc = A[0]; for (double &c : A) c /= lc; }
    if (A.empty()) A.push_back(0.0);
    *len_out = (int)A.size();
    for (int i = 0; i < *len_out; ++i) out[i] = A[i];
    return 0;
}
static std::vector<std::vector<double>> build_sturm_seq(const std::vector<double>& P) {
    std::vector<std::vector<double>> seq;
    seq.push_back(P);
    if (P.size() <= 1) return seq;
    std::vector<double> d(std::max((int)P.size() - 1, 1));
    poly_derivative(P.data(), (int)P.size(), d.data());
    trim(d);
    seq.push_back(d);
    int i = 1;
    while ((int)seq[i].size() > 1) {
        std::vector<double> q(std::max((int)seq[i-1].size() - (int)seq[i].size() + 1, 1));
        std::vector<double> r(seq[i-1].size());
        int len_q=0, len_r=0;
        int rc = poly_div(seq[i-1].data(), (int)seq[i-1].size(), seq[i].data(), (int)seq[i].size(), q.data(), &len_q, r.data(), &len_r);
        if (rc != 0) break; 
        std::vector<double> R(r.begin(), r.begin() + len_r);
        trim(R);
        for (double &c : R) c = -c;
        seq.push_back(R);
        i++;
        if (i > 200) break; // safety limit
    }
    return seq;
}
static int count_sign_changes(const std::vector<std::vector<double>>& seq, double x) {
    int changes = 0;
    int prev = 0;
    for (auto &poly : seq) {
        double val = poly_eval(poly.data(), (int)poly.size(), x);
        int s = sign(val);
        if (prev != 0 && s != 0 && s != prev) changes++;
        if (s != 0) prev = s;
    }
    return changes;
}

int poly_sturm_sign_changes(const double* a, int len_a, double x) {
    std::vector<double> P(a, a + len_a); trim(P);
    if (P.empty()) return 0;
    auto seq = build_sturm_seq(P);
    return count_sign_changes(seq, x);
}
int poly_num_real_roots_interval(const double* a, int len_a, double left, double right) {
    if (left > right) std::swap(left, right);
    std::vector<double> P(a, a + len_a); trim(P);
    if (P.empty()) return 0;
    auto seq = build_sturm_seq(P);
    return count_sign_changes(seq, left) - count_sign_changes(seq, right);
}

// Helper for finding roots
static void find_intervals(const std::vector<std::vector<double>>& seq, double a, double b, std::vector<std::pair<double, double>>& intervals) {
    int sa = count_sign_changes(seq, a);
    int sb = count_sign_changes(seq, b);
    int num_roots = sa - sb;
    
    if (num_roots == 0) return;
    if (num_roots == 1) {
        intervals.push_back({a, b});
        return;
    }
    double mid = (a + b) / 2.0;
    if (std::fabs(mid - a) < 1e-14) { // Avoid infinite recursion
        intervals.push_back({a, b});
        return;
    }
    find_intervals(seq, a, mid, intervals);
    find_intervals(seq, mid, b, intervals);
}

int poly_find_real_roots(const double* a, int len_a, double* roots) {
    std::vector<double> P(a, a + len_a); trim(P);
    if (P.size() <= 1) return 0; // Constant or empty
    
    // Bound
    double max_coeff = 0.0;
    for (size_t i = 1; i < P.size(); ++i) max_coeff = std::max(max_coeff, std::fabs(P[i]));
    double bound = 1.0 + max_coeff / std::fabs(P[0]);
    
    auto seq = build_sturm_seq(P);
    std::vector<std::pair<double, double>> intervals;
    find_intervals(seq, -bound, bound, intervals);
    
    auto push_unique = [&](std::vector<double>& vec, double r){
        for(double v : vec){ if(std::fabs(v - r) < 1e-6) return; }
        vec.push_back(r);
    };
    std::vector<double> out_roots;
    for (auto &interval : intervals) {
        double left = interval.first;
        double right = interval.second;
        
        // Check endpoints
        if (std::fabs(poly_eval(P.data(), (int)P.size(), left)) < 1e-9) { push_unique(out_roots, left); continue; }
        if (std::fabs(poly_eval(P.data(), (int)P.size(), right)) < 1e-9) { push_unique(out_roots, right); continue; }

        // Root refinement guided by Sturm counts
        for (int i = 0; i < 100; ++i) {
            double mid = (left + right) / 2.0;
            double val = poly_eval(P.data(), (int)P.size(), mid);
            if (std::fabs(val) < 1e-9) { left = mid; right = mid; break; }
            int sa = count_sign_changes(seq, left);
            int sm = count_sign_changes(seq, mid);
            int num_left = sa - sm;
            if (num_left > 0) right = mid; else left = mid;
        }
        double root = (left + right) / 2.0;
        push_unique(out_roots, root);
    }
    std::sort(out_roots.begin(), out_roots.end());
    // Final dedup (safety)
    std::vector<double> uniq;
    for(double r : out_roots){ if(uniq.empty() || std::fabs(uniq.back()-r) > 1e-6) uniq.push_back(r); }
    for(size_t i=0;i<uniq.size();++i) roots[i] = uniq[i];
    return (int)uniq.size();
}

int poly_bezout(const double* a, int len_a, const double* b, int len_b,
                double* gcd, int* len_gcd,
                double* s, int* len_s,
                double* t, int* len_t) {
    std::vector<double> r0(a, a + len_a); trim(r0);
    std::vector<double> r1(b, b + len_b); trim(r1);
    
    std::vector<double> s0(1, 1.0), s1(1, 0.0);
    std::vector<double> t0(1, 0.0), t1(1, 1.0);
    
    while (!is_zero_vec(r1)) {
        // q = r0 / r1
        std::vector<double> q(std::max((int)r0.size() - (int)r1.size() + 1, 1));
        std::vector<double> r_tmp(r0.size());
        int lq=0, lr=0;
        int rc = poly_div(r0.data(), (int)r0.size(), r1.data(), (int)r1.size(), q.data(), &lq, r_tmp.data(), &lr);
        if (rc != 0) return rc;
        trim(q);
        
        // r0 = r1, r1 = remainder (which is r0 - q*r1)
        // Actually poly_div returns remainder directly.
        // But we need to update s and t.
        // s_new = s0 - q*s1
        // t_new = t0 - q*t1
        
        // r_new is r_tmp (trimmed)
        std::vector<double> r_new(r_tmp.begin(), r_tmp.begin() + lr); trim(r_new);
        
        // Compute s_new
        // q * s1
        std::vector<double> qs1(q.size() + s1.size() - 1);
        poly_mul(q.data(), (int)q.size(), s1.data(), (int)s1.size(), qs1.data());
        // s0 - qs1
        std::vector<double> s_new(std::max(s0.size(), qs1.size()));
        poly_sub(s0.data(), (int)s0.size(), qs1.data(), (int)qs1.size(), s_new.data());
        trim(s_new);
        
        // Compute t_new
        std::vector<double> qt1(q.size() + t1.size() - 1);
        poly_mul(q.data(), (int)q.size(), t1.data(), (int)t1.size(), qt1.data());
        std::vector<double> t_new(std::max(t0.size(), qt1.size()));
        poly_sub(t0.data(), (int)t0.size(), qt1.data(), (int)qt1.size(), t_new.data());
        trim(t_new);
        
        r0 = r1; r1 = r_new;
        s0 = s1; s1 = s_new;
        t0 = t1; t1 = t_new;
    }
    
    // Normalize gcd to be monic
    if (!is_zero_vec(r0)) {
        double lc = r0[0];
        for (double &c : r0) c /= lc;
        for (double &c : s0) c /= lc;
        for (double &c : t0) c /= lc;
    } else {
        if(r0.empty()) r0.push_back(0.0);
    }
    
    if (r0.empty()) r0.push_back(0.0);
    if (s0.empty()) s0.push_back(0.0);
    if (t0.empty()) t0.push_back(0.0);
    
    *len_gcd = (int)r0.size();
    for(int i=0; i<*len_gcd; ++i) gcd[i] = r0[i];
    
    *len_s = (int)s0.size();
    for(int i=0; i<*len_s; ++i) s[i] = s0[i];
    
    *len_t = (int)t0.size();
    for(int i=0; i<*len_t; ++i) t[i] = t0[i];
    
    return 0;
}
}
void poly_mul_fft(const double* a, int len_a, const double* b, int len_b, double* out) {
    int len = len_a + len_b - 1;
    size_t n = next_pow2((size_t)len);
    std::vector<std::complex<double>> A(n), B(n);
    for (int i = 0; i < len_a; ++i) A[i] = std::complex<double>(a[len_a - 1 - i], 0.0);
    for (int i = len_a; i < (int)n; ++i) A[i] = std::complex<double>(0.0, 0.0);
    for (int i = 0; i < len_b; ++i) B[i] = std::complex<double>(b[len_b - 1 - i], 0.0);
    for (int i = len_b; i < (int)n; ++i) B[i] = std::complex<double>(0.0, 0.0);
    fft(A, false); fft(B, false);
    for (size_t i = 0; i < n; ++i) A[i] *= B[i];
    fft(A, true);
    std::vector<double> rev(len);
    for (int i = 0; i < len; ++i) rev[i] = A[i].real();
    clamp_eps(rev);
    for (int i = 0; i < len; ++i) out[i] = rev[len - 1 - i];
}
extern "C" {
void poly_set_fft_threshold(int thr) { if (thr < 0) thr = 0; g_fft_threshold = thr; }
}
