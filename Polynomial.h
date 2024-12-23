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
    std::vector<int> coeffs;//�洢����ʽ��ϵ��
    unsigned int degree;//�洢����ʽ�Ľ��� = NTRU static N
    unsigned int q;//����ʽģ�� = NTRU static q
    unsigned int p;
    unsigned int cnt_1 = -1;//�������ʽ��cnt_1��1,���ڸ�һ����δ֪
    unsigned int cnt_neg1 = -1;//�������ʽ��cnt_neg1��-1
    //��Ϊ����ʽ��Ҫmod q mod X^N - 1;
    ll PowerMod(ll aa, ll n, ll mod);

public:
    Polynomial() {};
    Polynomial(unsigned int d, unsigned int q, unsigned int p, unsigned int a, unsigned int b, bool if_Need_inversable = false) ;
    Polynomial(unsigned int d, unsigned int q, unsigned int p, std::vector<int> vec) : degree(d), q(q), p(p),coeffs(vec) {};

    std::vector<int> gen_Pol(unsigned int d, unsigned a, unsigned b);
    std::vector<int> gen_Pol_inverse(unsigned int d, unsigned a, unsigned b);
    bool test_If_Inversable(unsigned int mod);
    bool test_If_Inversable(const std::vector<int>& coeffs, unsigned int mod);
    // �������ʽ��ģ���� (�� x^d - 1 ȡģ)
    std::vector<int> poly_mod_XN(std::vector<int>& poly, unsigned int d, unsigned int p);
    std::vector<int> poly_mod_p(std::vector<int>& poly, unsigned int p);
    // ������������ʽ���������
    std::vector<int> poly_gcd(const std::vector<int>& a, const std::vector<int>& b);
    // ����������ȥ������ʽǰ����
    std::vector<int> remove_leading_zeros(const std::vector<int>& poly);


    Polynomial operator+(const Polynomial& other) const; //����ʽ��ͣ����ط��� //ǰһ��const��������Ĳ����ڶ���const��������ĵ��ú����Ķ���
    Polynomial Pol_add(const Polynomial& p) const;  //����ʽ���

    Polynomial operator-(const Polynomial& other) const;//����ʽ���
    Polynomial Pol_sub(const Polynomial& p) const;//����ʽ���

    void FFT(std::vector <cd>& a) ;//����Ҷ�任����;��
    void IFFT(std::vector <cd>& a) ;//���ǻ���Щ���⣬���ڻ�����ʵ��������ʵ�ֵľ��

    Polynomial operator*(const Polynomial& other) ;//����ʽ���
    Polynomial Pol_mul(const Polynomial& p) ;//���ǻ���Щ���⣬���ڻ�����ʵ��������ʵ�ֵľ��
    
    std::vector<int> get_coeffs() const;
    void print_pol();
    Polynomial Cal_Inverse_Of_Polynomial(std::vector<int> pol, unsigned int q, unsigned int degree);
    Polynomial Cal_Inverse_Of_Polynomial(unsigned int mod);

    template<typename T>
    T Cal_Inverse_Of_Polynomial(const T& obj, unsigned q, unsigned int degree);
};

#endif