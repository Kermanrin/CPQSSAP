#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <iostream>
#include <complex>
#include <cmath>

typedef std::complex<double> cd;

class Polynomial {
private:
    std::vector<int> coeffs;//存储多项式的系数
    unsigned int degree;//存储多项式的阶数
    unsigned int q;//多项式模数 = NTRU static q
    //因为多项式需要mod q mod X^N - 1;

public:
    Polynomial() {}
    Polynomial(unsigned int d, unsigned int _q) : coeffs(degree + 1, 0), degree(d), q(_q) {};
    Polynomial(std::vector<int> c, unsigned int _q);
    Polynomial(unsigned int degree, std::vector<int> c, unsigned int _q);
    std::vector<int> Cal_Inverse_Of_Polynomial(std::vector<int> pol, unsigned int q);

    Polynomial operator+(const Polynomial& other) const; //多项式求和，重载方法 //前一个const代表不会更改参数第二个const代表不会更改调用函数的对象
    Polynomial Pol_add(const Polynomial& p) const;  //多项式求和

    Polynomial operator-(const Polynomial& other) const;//多项式求差
    Polynomial Pol_sub(const Polynomial& p) const;//多项式求差

    void FFT(std::vector <cd>& a) ;//傅里叶变换求逆和卷积
    void IFFT(std::vector <cd>& a) ;//但是还有些问题，现在还是在实数域上面实现的卷积

    Polynomial operator*(const Polynomial& other) ;//多项式卷积
    Polynomial Pol_mul(const Polynomial& p) ;//但是还有些问题，现在还是在实数域上面实现的卷积
    
    std::vector<int> get_coeffs() const;
    void print_pol();

    template<typename T>
    T Cal_Inverse_Of_Polynomial(const T& obj);
};

#endif