#include "Polynomial.h"
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <random>
#include <complex>
#include <cmath>
#include <assert.h>
#include "NTRU.h"

typedef std::complex<double> cd;

Polynomial::Polynomial(unsigned int d, std::vector<int> c, unsigned int _q) {
    assert((unsigned int)c.size() == d);
    degree = d;
    q = _q;
    for (int i = 0; i < d; i++) {
        coeffs.push_back(c[i]);
    }
};
Polynomial::Polynomial(std::vector<int> c, unsigned int _q) :q(_q), degree(c.size()){
    for (int i = 0; i < c.size(); i++) {
        this->coeffs.push_back(c[i]);
    }
};

Polynomial::Polynomial(unsigned int d, unsigned int _q, unsigned int p, unsigned int a, bool if_Need_inversable): degree(d), q(_q), p(p), a(a), b(a) {

    if (if_Need_inversable) {
        coeffs = gen_Pol_inverse(d, a, b);
    }
    else {
        coeffs = gen_Pol(d, a, b);
    }
}

Polynomial::Polynomial(unsigned int d, unsigned int _q, unsigned int p, unsigned int a, unsigned int b, bool if_Need_inversable) : degree(d), q(_q), p(p), a(a), b(b) {

    if (if_Need_inversable) {
        coeffs = gen_Pol_inverse(d, a, b);
    }
    else {
        coeffs = gen_Pol(d, a, b);
    }
}


std::vector<int> Polynomial::gen_Pol(unsigned int d, unsigned a, unsigned b) {
    std::vector<int> init_coeffs(d, 0);
    for (int i = 0; i < (int)a; i++) init_coeffs[i] = 1;
    for (int i = a; i < (int)a + b; i++) init_coeffs[i] = -1;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(init_coeffs.begin(), init_coeffs.end(), gen);

    return init_coeffs;
}


std::vector<int> Polynomial::gen_Pol_inverse(unsigned int d, unsigned a, unsigned b) {
    std::vector<int> vec = gen_Pol(d, a, b);
    for (int i = 0; i < d; i++) std::cout << vec[i] << " ";
    std::cout << std::endl;
    while (!test_If_Inversable(vec, d)) {
        vec = gen_Pol(d, a, b);
    }
    return vec;
}

bool Polynomial::test_If_Inversable(const std::vector<int>& coeffs, unsigned int d) {
    // 构造多项式 x^d - 1
    std::vector<int> x_d_minus_1(d, 0);
    x_d_minus_1[0] = 1;
    x_d_minus_1[d - 1] = -1;

    // 计算 gcd(f(x), x^d - 1)
    std::vector<int> gcd = Polynomial::poly_gcd(coeffs, x_d_minus_1);

    // 判断 gcd 是否为常数 1
    return gcd.size() == 1 && gcd[0] == 1;
}


std::vector<int> Polynomial::poly_gcd(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> poly_a = a;
    std::vector<int> poly_b = b;

    while (!poly_b.empty() && !(poly_b.size() == 1 && poly_b[0] == 0)) {
        // 计算 a mod b
        std::vector<int> remainder = poly_mod(poly_a, poly_b.size() - 1);

        // 更新 poly_a 和 poly_b
        poly_a = poly_b;
        poly_b = remainder;
    }

    return remove_leading_zeros(poly_a); // 返回去除前导零的 gcd
}


std::vector<int> Polynomial::remove_leading_zeros(const std::vector<int>& poly) {
    auto it = std::find_if(poly.rbegin(), poly.rend(), [](int coef) { return coef != 0; });
    return std::vector<int>(poly.begin(), it.base());
}

std::vector<int> Polynomial::poly_mod(const std::vector<int>& poly, unsigned int d) {
    std::vector<int> result = poly;
    while (result.size() > d) {
        // 将最高次项减去
        int degree = result.size() - 1;
        int coef = result[degree];
        for (int i = 0; i < d; ++i) {
            result[degree - i] -= coef * (i == 0 ? 1 : 0); // 模掉最高次项
        }
        result.pop_back();
    }
    return remove_leading_zeros(result);
}




Polynomial Polynomial::operator+(const Polynomial& other) const {
    unsigned int deg = std::max(coeffs.size(), other.coeffs.size());
    unsigned int mod = NTRU::q;
    std::vector<int> res(deg, 0);
    int cnt = 0;
    while (cnt < coeffs.size() && cnt < other.coeffs.size()) {
        res[cnt] = (this->coeffs[cnt] + other.coeffs[cnt]) % mod;
        cnt++;
    }
    while (cnt < coeffs.size()) {
        res[cnt] = this->coeffs[cnt];
        cnt++;
    }
    while (cnt < other.coeffs.size()) {
        res[cnt] = other.coeffs[cnt];
        cnt++;
    }
    return Polynomial(deg, res, q);
}

Polynomial Polynomial::Pol_add(const Polynomial& other) const {
    unsigned int deg = std::max(coeffs.size(), other.coeffs.size());
    unsigned int mod = NTRU::q;
    std::vector<int> res(deg, 0);
    int cnt = 0;
    while (cnt < coeffs.size() && cnt < other.coeffs.size()) {
        res[cnt] = (this->coeffs[cnt] + other.coeffs[cnt]) % mod;
        cnt++;
    }
    while (cnt < coeffs.size()) {
        res[cnt] = this->coeffs[cnt];
        cnt++;
    }
    while (cnt < other.coeffs.size()) {
        res[cnt] = other.coeffs[cnt];
        cnt++;
    }
    return Polynomial(deg, res, q);
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
    unsigned int deg = std::max(coeffs.size(), other.coeffs.size());
    unsigned int mod = NTRU::q;
    std::vector<int> res(deg, 0);
    int cnt = 0;
    while (cnt < coeffs.size() && cnt < other.coeffs.size()) {
        res[cnt] = (this->coeffs[cnt] - other.coeffs[cnt] + mod) % mod;
        cnt++;
    }
    while (cnt < coeffs.size()) {
        res[cnt] = this->coeffs[cnt];
        cnt++;
    }
    while (cnt < other.coeffs.size()) {
        res[cnt] = other.coeffs[cnt];
        cnt++;
    }
    return Polynomial(deg, res, q);
}

Polynomial Polynomial::Pol_sub(const Polynomial& other) const {
    unsigned int deg = std::max(coeffs.size(), other.coeffs.size());
    unsigned int mod = NTRU::q;
    std::vector<int> res(deg, 0);
    int cnt = 0;
    while (cnt < coeffs.size() && cnt < other.coeffs.size()) {
        res[cnt] = (this->coeffs[cnt] - other.coeffs[cnt] + mod) % mod;
        cnt++;
    }
    while (cnt < coeffs.size()) {
        res[cnt] = this->coeffs[cnt];
        cnt++;
    }
    while (cnt < other.coeffs.size()) {
        res[cnt] = other.coeffs[cnt];
        cnt++;
    }
    return Polynomial(deg, res, q);
}

std::vector<int> Polynomial::get_coeffs() const {
    return coeffs;
}

void Polynomial::print_pol() {
    std::vector<int> tmp = this->get_coeffs();
    printf("The sk is : ");
    for (int i = 0; i < tmp.size(); i++) {
        printf("%d ", tmp[i]);
    }
    printf("\n");
}

void Polynomial::FFT(std::vector <cd>& a) {
    int n = a.size();
    if (n <= 1) return;

    std::vector<cd> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; i++) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    FFT(even);
    FFT(odd);

    double ang = 2 * std::acos(-1) / n;
    cd w(1), wn(cos(ang), sin(ang));
    for (int i = 0; i < n / 2; i++) {
        a[i] = even[i] + w * odd[i];
        a[i + n / 2] = even[i] - w * odd[i];
        w *= wn;
    }
}

void Polynomial::IFFT(std::vector <cd>& a) {
    int n = a.size();
    for (int i = 0; i < n; i++) a[i] = std::conj(a[i]);
    IFFT(a);
    for (int i = 0; i < n; i++) a[i] = std::conj(a[i]) / double(n);
}

Polynomial Polynomial::operator*(const Polynomial& other)  {
    std::vector<int> A = get_coeffs();
    std::vector<int> B = other.get_coeffs();
    int m = A.size();
    int n = B.size();
    std::vector<int> C(NTRU::N, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n && i + j < NTRU::N; j++) {
            C[i + j] += A[i] * B[j];
        }
        for (int j = std::max(0, int(NTRU::N - i)); j < n; j++) {
            C[i + j - NTRU::N] += A[i] * B[j];
        }
    }
    return Polynomial(C, NTRU::q);
}
Polynomial Polynomial::Pol_mul(const Polynomial& other) {
    std::vector<int> A = get_coeffs();
    std::vector<int> B = other.get_coeffs();
    int m = A.size();
    int n = B.size();
    std::vector<int> C(NTRU::N, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n && i + j < NTRU::N; j++) {
            C[i + j] += A[i] * B[j];
        }
        for (int j = std::max(0, int(NTRU::N - i)); j < n; j++) {
            C[i + j - NTRU::N] += A[i] * B[j];
        }
    }
    return Polynomial(C, NTRU::q);
}