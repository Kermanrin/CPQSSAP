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
        poly[j] = (poly[j] + poly[i]) % p;
    }
    for (int i = 0; i < d; i++) {
        poly[i] %= p;
    }
    poly.resize(d);
    return poly;
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
    const int n = degree;
    std::cout << "n:" << n << std::endl;
    std::vector<int> a(3 * n + 1, 0);
    for (int i = 0; i < n; i++) a[i] = coeffs[i];
    std::vector<int> b(3 * n + 1, 0);
    int g = 3, gi = 46;
    std::vector<int> rev(3 * n + 1, 0);
    std::vector<int> ta(3 * n + 1, 0);
    std::vector<int> tb(3 * n + 1, 0);

    NTT_solve(n, a, b, mod, rev, ta, tb, g, gi);
    //for (int i = 0; i < 3 * degree + 1; i++) std::cout << b[i] << " ";
    //std::cout << std::endl;
    b.resize(n);
    for (int i = 0; i < n; i++) b[i] -= mod / 2;

    return Polynomial(degree, q, p, b);
}

std::vector<int> Polynomial::CalVec_Inverse_Of_Polynomial(unsigned int mod) {
    const int n = degree;
    std::cout << "n:" << n << std::endl;
    std::vector<int> a(3 * n + 1, 0);
    for (int i = 0; i < n; i++) a[i] = coeffs[i];
    std::vector<int> b(3 * n + 1, 0);
    int g = 3, gi = 46;
    std::vector<int> rev(3 * n + 1, 0);
    std::vector<int> ta(3 * n + 1, 0);
    std::vector<int> tb(3 * n + 1, 0);

    NTT_solve(n, a, b, mod, rev, ta, tb, g, gi);
    //for (int i = 0; i < 3 * degree + 1; i++) std::cout << b[i] << " ";
    //std::cout << std::endl;
    b.resize(n);
    for (int i = 0; i < n; i++) b[i] -= mod / 2;

    return b;
}

Polynomial Polynomial::Cal_Inverse_Of_Polynomial(std::vector<int> pol, unsigned int q, unsigned int degree) {
    return Cal_Inverse_Of_Polynomial(q);
    
}

void Polynomial::NTT_solve(int len, std::vector<int>& a, std::vector<int>& b, int mod, std::vector<int> &rev, std::vector<int>& ta, std::vector<int>& tb,
    int g, int gi) {
    if (len == 1) {
        b[0] = Poly_qmi(a[0], mod - 2, mod);
        return;
    }
    NTT_solve((len + 1) >> 1, a, b, mod, rev, ta, tb, g, gi);
    int l = 1, bit = 0;
    while(l <= len + degree) {
        l <<= 1;
        ++bit;
    }
    NTT_pre(bit, rev);
    for (int i = 0; i < l; ++i) {
        ta[i] = a[i];
        tb[i] = (i < ((len + 1) >> 1) ? b[i] : 0);
    }
    NTT(ta, l, 1, rev, g, gi, mod);
    NTT(tb, l, 1, rev, g, gi, mod);
    for (int i = 0; i < l; ++i) {
        ta[i] = tb[i] * (((2 - ta[i] * tb[i]) % mod + mod) % mod) % mod;
    }
    NTT(ta, l, 0, rev, g, gi, mod);
    int inv = Poly_qmi(l, mod - 2, mod);
    for (int i = 0; i < len; ++i) {
        b[i] = ta[i] * inv % mod;
    }
}

int Polynomial::Poly_qmi(int a, int n, int mod) {
    int res = 1;
    while (n) {
        if (n & 1) {
            a = a * res % mod;
        }
        a = a * a % mod;
        n >>= 1;
    }
    return res;
}

void Polynomial::NTT_pre(int bit, std::vector<int>& rev) {
    for (int i = 0; i < (1 << bit); ++i) {
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
    }
}

void Polynomial::NTT(std::vector<int>& F, int len, int on, std::vector<int> & rev, int g, int gi, int mod) {
    for (int i = 0; i < len; ++i) {
        if (i < rev[i]) {
            std::swap(F[i], F[rev[i]]);
        }
    }
    for (int i = 2; i <= len; i <<= 1) {
        int gn = Poly_qmi(on ? g : gi, (mod - 1) / i, mod);
        for (int j = 0; j < len; j += i) {
            int gg = 1;
            for (int k = j; k < j + i / 2; k++) {
                int u = F[k];
                int v = gg * F[k + i / 2] % mod;
                F[k] = (u + v) % mod;
                F[k + i / 2] = (u - v + mod) % mod;
                gg = gg * gn % mod;
            }
        }
    }
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
    std::cout << tmp.size() << std::endl;
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
    int deg = A.size();
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
    return Polynomial(deg, q, p, C);
}
Polynomial Polynomial::Pol_mul(const Polynomial& other, int mod) {
    std::vector<int> A = this->coeffs;
    int d = A.size();
    std::vector<int> B = other.get_coeffs();
    int m = A.size();
    
    int n = B.size();
    std::cout << m << " " << n << std::endl;
    std::vector<int> C(251 * 2 - 1);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i + j] += A[i] * B[j];
        }
    }
    return Polynomial(251, q, p, poly_mod_XN(C, 251, mod));
}

std::vector<int> Polynomial::Pol_mul_Vec(const Polynomial& other, int mod) {
    std::vector<int> A = this->coeffs;
    int d = A.size();
    std::vector<int> B = other.get_coeffs();
    int m = A.size();

    int n = B.size();
    std::cout << m << " " << n << std::endl;
    std::vector<int> C(degree * 2 - 1);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i + j] += A[i] * B[j];
        }
    }
    return poly_mod_XN(C, d, mod);
}