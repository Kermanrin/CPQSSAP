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
typedef long long ll;

Polynomial::Polynomial(unsigned int d, unsigned int q, unsigned int p, unsigned int a, unsigned int b, bool if_Need_inversable) : degree(d), q(q), p(p), cnt_1(a), cnt_neg1(b) {

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

    while (!test_If_Inversable(vec, d)) {
        vec = gen_Pol(d, a, b);
    }
    return vec;
}

bool Polynomial::test_If_Inversable(const std::vector<int>& coeffs, unsigned int d) {
    return coeffs[0] != 0;
}

bool Polynomial::test_If_Inversable(const unsigned int d) {
    std::vector<int> vec = this->coeffs;
    return test_If_Inversable(vec, d);
}

std::vector<int> Polynomial::poly_mod_XN(std::vector<int>& poly, unsigned int d, unsigned int p) {
    int len = poly.size();
    for (int i = d; i < len; i++) {
        int j = i % d;
        if (poly[i] == 1) {
            poly[j] += 1;
        }
        else if (poly[i] == -1) {
            poly[j] += -1;
        }
    }
    poly.resize(d);
    return poly_mod_p(poly, p);
}

std::vector<int> Polynomial::poly_mod_p(std::vector<int>& poly, unsigned int p) {
    int len = poly.size();
    for (int i = 0; i < len; i++) {
        poly[i] %= 3;
        poly[i] = (poly[i] + 3) % 3;
        if (poly[i] == 2) poly[i] = -1;
    }
    return poly;
}

ll Polynomial::PowerMod(ll aa, ll nn, ll mod) {
    ll c = 1;
    while (nn) {
        if (nn & 1) {
            c = c * aa % mod;
        }
        aa = aa * aa % mod;
        nn >>= 1;
    }
    return c;
}

Polynomial Polynomial::Cal_Inverse_Of_Polynomial(unsigned int mod) {
    return Polynomial(degree, q, p, {0});
}

Polynomial Polynomial::Cal_Inverse_Of_Polynomial(std::vector<int> pol, unsigned int q, unsigned int degree) {
    return Polynomial(degree, q, p, { 0 });
    
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
    unsigned int deg = std::max(coeffs.size(), other.coeffs.size());
    unsigned int mod = q;
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
    return Polynomial(deg, q, p, res);
}

Polynomial Polynomial::Pol_add(const Polynomial& other) const {
    unsigned int deg = std::max(coeffs.size(), other.coeffs.size());
    unsigned int mod = q;
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
    return Polynomial(deg, q, p, res);
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
    unsigned int deg = std::max(coeffs.size(), other.coeffs.size());
    unsigned int mod = q;
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
    return Polynomial(deg, q, p, res);
}

Polynomial Polynomial::Pol_sub(const Polynomial& other) const {
    unsigned int deg = std::max(coeffs.size(), other.coeffs.size());
    unsigned int mod = q;
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
    return Polynomial(deg, q, p, res);
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
    std::vector<int> C(degree, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n && i + j < degree; j++) {
            C[i + j] += A[i] * B[j];
        }
        for (int j = std::max(0, int(degree - i)); j < n; j++) {
            C[i + j - degree] += A[i] * B[j];
        }
    }
    int deg = C.size();
    return Polynomial(deg, q, p, C);
}
Polynomial Polynomial::Pol_mul(const Polynomial& other) {
    /*std::vector<int> A = get_coeffs();
    std::vector<int> B = other.get_coeffs();
    int m = A.size();
    int n = B.size();
    std::cout << "m: " << m << " n : " << n << std::endl;
    std::vector<int> C(NTRU::N, 0);
    std::cout << "NTRUN= " << NTRU::N << std::endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n && i + j < NTRU::N; j++) {
            C[i + j] += A[i] * B[j];
        }
        for (int j = std::max(0, int(NTRU::N - i)); j < n; j++) {
            C[i + j - NTRU::N] += A[i] * B[j];
        }
    }
    int deg = C.size();
    return Polynomial(deg, this->q, this->p, -1, -1, C);*/
    int d = this->degree;
    int p = this->p;
    std::vector<int> A = get_coeffs();
    std::vector<int> B = other.get_coeffs();
    int m = A.size();
    int n = B.size();
    std::vector<int> C(degree * 2 - 1);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i + j] += A[i] * B[j];
        }
    }
    return Polynomial(d, q, p, poly_mod_XN(C, d, p));
}