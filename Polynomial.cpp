#include "Polynomial.h"
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
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