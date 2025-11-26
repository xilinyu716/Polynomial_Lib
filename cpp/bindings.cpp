#include "polynomial.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

static std::vector<double> wrap_add(const std::vector<double>& a, const std::vector<double>& b){
    int len = std::max((int)a.size(), (int)b.size());
    std::vector<double> out(len);
    poly_add(a.data(), (int)a.size(), b.data(), (int)b.size(), out.data());
    return out;
}
static std::vector<double> wrap_sub(const std::vector<double>& a, const std::vector<double>& b){
    int len = std::max((int)a.size(), (int)b.size());
    std::vector<double> out(len);
    poly_sub(a.data(), (int)a.size(), b.data(), (int)b.size(), out.data());
    return out;
}
static std::vector<double> wrap_mul(const std::vector<double>& a, const std::vector<double>& b){
    int len = (int)a.size() + (int)b.size() - 1;
    std::vector<double> out(std::max(len, 1));
    poly_mul(a.data(), (int)a.size(), b.data(), (int)b.size(), out.data());
    return out;
}
static py::tuple wrap_div(const std::vector<double>& a, const std::vector<double>& b){
    std::vector<double> q(std::max((int)a.size() - (int)b.size() + 1, 1));
    std::vector<double> r((int)a.size());
    int len_q=0, len_r=0;
    int rc = poly_div(a.data(), (int)a.size(), b.data(), (int)b.size(), q.data(), &len_q, r.data(), &len_r);
    if (rc != 0) throw std::runtime_error("division failed");
    std::vector<double> qv(q.begin(), q.begin()+len_q);
    std::vector<double> rv(r.begin(), r.begin()+len_r);
    return py::make_tuple(qv, rv);
}
static std::vector<double> wrap_derivative(const std::vector<double>& a){
    std::vector<double> out(std::max((int)a.size()-1, 1));
    poly_derivative(a.data(), (int)a.size(), out.data());
    return out;
}
static double wrap_eval(const std::vector<double>& a, double x){
    return poly_eval(a.data(), (int)a.size(), x);
}
static std::vector<double> wrap_gcd(const std::vector<double>& a, const std::vector<double>& b){
    std::vector<double> out(std::max((int)a.size(), (int)b.size()));
    int len_out=0; int rc = poly_gcd(a.data(), (int)a.size(), b.data(), (int)b.size(), out.data(), &len_out);
    if (rc != 0) throw std::runtime_error("gcd failed");
    return std::vector<double>(out.begin(), out.begin()+len_out);
}
static int wrap_sturm_sign_changes(const std::vector<double>& a, double x){
    return poly_sturm_sign_changes(a.data(), (int)a.size(), x);
}
static int wrap_num_real_roots_interval(const std::vector<double>& a, double left, double right){
    return poly_num_real_roots_interval(a.data(), (int)a.size(), left, right);
}
static std::vector<double> wrap_find_real_roots(const std::vector<double>& a){
    std::vector<double> roots(std::max((int)a.size(), 1));
    int count = poly_find_real_roots(a.data(), (int)a.size(), roots.data());
    return std::vector<double>(roots.begin(), roots.begin()+count);
}
static py::tuple wrap_bezout(const std::vector<double>& a, const std::vector<double>& b){
    std::vector<double> gcd(std::max((int)a.size(), (int)b.size()));
    std::vector<double> s((int)b.size()+1);
    std::vector<double> t((int)a.size()+1);
    int len_g=0, len_s=0, len_t=0;
    int rc = poly_bezout(a.data(), (int)a.size(), b.data(), (int)b.size(),
                         gcd.data(), &len_g, s.data(), &len_s, t.data(), &len_t);
    if (rc != 0) throw std::runtime_error("bezout failed");
    return py::make_tuple(
        std::vector<double>(gcd.begin(), gcd.begin()+len_g),
        std::vector<double>(s.begin(), s.begin()+len_s),
        std::vector<double>(t.begin(), t.begin()+len_t)
    );
}

PYBIND11_MODULE(polylib_ext, m){
    m.def("add", &wrap_add);
    m.def("sub", &wrap_sub);
    m.def("mul", &wrap_mul);
    m.def("div", &wrap_div);
    m.def("derivative", &wrap_derivative);
    m.def("eval", &wrap_eval);
    m.def("gcd", &wrap_gcd);
    m.def("sturm_sign_changes", &wrap_sturm_sign_changes);
    m.def("num_real_roots_interval", &wrap_num_real_roots_interval);
    m.def("find_real_roots", &wrap_find_real_roots);
    m.def("bezout", &wrap_bezout);
    m.def("set_fft_threshold", [](int thr){ poly_set_fft_threshold(thr); });
}

