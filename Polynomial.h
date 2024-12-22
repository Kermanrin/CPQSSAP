#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <iostream>
#include <complex>
#include <cmath>

typedef std::complex<double> cd;
typedef long long ll;

class Polynomial {
private:
    std::vector<int> coeffs;//存储多项式的系数
    unsigned int degree;//存储多项式的阶数 = NTRU static N
    unsigned int q;//多项式模数 = NTRU static q
    unsigned int p;
    unsigned int cnt_1;//代表多项式有a个1
    unsigned int cnt_neg1;//代表多项式有b个-1
    //因为多项式需要mod q mod X^N - 1;
    ll PowerMod(ll aa, ll n, ll mod);

public:
    Polynomial() {};
    Polynomial(unsigned int d, unsigned int _q, unsigned int p, unsigned int a, unsigned int b, std::vector<int> vec): degree(d), q(_q), p(p), cnt_1(a), cnt_neg1(b), coeffs(vec) {};
    Polynomial(unsigned int d, unsigned int _q, unsigned int p) : coeffs(gen_Pol(d, 1, 1)), degree(d), q(_q), p(p) {};
    Polynomial(unsigned int d, unsigned int _q, unsigned int p, unsigned int a, unsigned int b, bool if_Need_inversable = false) ;
    Polynomial(unsigned int d, unsigned int _q, unsigned int p, unsigned int a, bool if_Need_inversable = false);

    std::vector<int> gen_Pol(unsigned int d, unsigned a, unsigned b);
    std::vector<int> gen_Pol_inverse(unsigned int d, unsigned a, unsigned b);
    bool test_If_Inversable(const std::vector<int>& coeffs, unsigned int d);
    // 计算多项式的模运算 (对 x^d - 1 取模)
    std::vector<int> poly_mod_XN(std::vector<int>& poly, unsigned int d, unsigned int p);
    std::vector<int> poly_mod_p(std::vector<int>& poly, unsigned int p);
    // 计算两个多项式的最大公因数
    std::vector<int> poly_gcd(const std::vector<int>& a, const std::vector<int>& b);
    // 辅助函数：去除多项式前导零
    std::vector<int> remove_leading_zeros(const std::vector<int>& poly);

    Polynomial(std::vector<int> c, unsigned int _q);//仅用于测试，多项式作为密钥应该是随机生成的
    Polynomial(unsigned int degree, std::vector<int> c, unsigned int _q);//仅用于测试，多项式作为密钥应该是随机生成的

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
    std::vector<int> Cal_Inverse_Of_Polynomial(std::vector<int> pol, unsigned int q, unsigned int degree);
    Polynomial Cal_Inverse_Of_Polynomial(Polynomial poly, unsigned q, unsigned int degree);

    template<typename T>
    T Cal_Inverse_Of_Polynomial(const T& obj, unsigned q, unsigned int degree);
};

#endif