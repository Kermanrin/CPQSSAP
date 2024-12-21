#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <iostream>
#include <complex>
#include <cmath>

typedef std::complex<double> cd;

class Polynomial {
private:
    std::vector<int> coeffs;//�洢����ʽ��ϵ��
    unsigned int degree;//�洢����ʽ�Ľ��� = NTRU static N
    unsigned int q;//����ʽģ�� = NTRU static q
    unsigned int p;
    unsigned int a;//�������ʽ��a��1
    unsigned int b;//�������ʽ��b��-1
    //��Ϊ����ʽ��Ҫmod q mod X^N - 1;

public:
    Polynomial() {};
    Polynomial(unsigned int d, unsigned int _q, unsigned int p) : coeffs(gen_Pol(d, 1, 1)), degree(d), q(_q), p(p) {};
    Polynomial(unsigned int d, unsigned int _q, unsigned int p, unsigned int a, unsigned int b, bool if_Need_inversable = false) ;
    Polynomial(unsigned int d, unsigned int _q, unsigned int p, unsigned int a, bool if_Need_inversable = false);

    bool test_If_Inversable(const std::vector<int>& coeffs, unsigned int d);
    std::vector<int> gen_Pol(unsigned int d, unsigned a, unsigned b);
    std::vector<int> gen_Pol_inverse(unsigned int d, unsigned a, unsigned b);
    // �������ʽ��ģ���� (�� x^d - 1 ȡģ)
    std::vector<int> poly_mod(const std::vector<int>& poly, unsigned int d);
    // ������������ʽ���������
    std::vector<int> poly_gcd(const std::vector<int>& a, const std::vector<int>& b);
    // ����������ȥ������ʽǰ����
    std::vector<int> remove_leading_zeros(const std::vector<int>& poly);

    Polynomial(std::vector<int> c, unsigned int _q);//�����ڲ��ԣ�����ʽ��Ϊ��ԿӦ����������ɵ�
    Polynomial(unsigned int degree, std::vector<int> c, unsigned int _q);//�����ڲ��ԣ�����ʽ��Ϊ��ԿӦ����������ɵ�

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
    std::vector<int> Cal_Inverse_Of_Polynomial(std::vector<int> pol, unsigned int q);

    template<typename T>
    T Cal_Inverse_Of_Polynomial(const T& obj);
};

#endif