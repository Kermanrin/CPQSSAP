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
    unsigned int degree;//�洢����ʽ�Ľ���
    unsigned int q;//����ʽģ�� = NTRU static q
    //��Ϊ����ʽ��Ҫmod q mod X^N - 1;

public:
    Polynomial() {}
    Polynomial(unsigned int d, unsigned int _q) : coeffs(degree + 1, 0), degree(d), q(_q) {};
    Polynomial(std::vector<int> c, unsigned int _q);
    Polynomial(unsigned int degree, std::vector<int> c, unsigned int _q);
    std::vector<int> Cal_Inverse_Of_Polynomial(std::vector<int> pol, unsigned int q);

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

    template<typename T>
    T Cal_Inverse_Of_Polynomial(const T& obj);
};

#endif