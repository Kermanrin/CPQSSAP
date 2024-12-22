#pragma once
#ifndef NTRU_H
#define NTRU_H

#include "Polynomial.h"

class NTRU {
public:
    NTRU(const unsigned int N , const unsigned int q, const unsigned int p) : f(Polynomial(N, q, p)), g(Polynomial(N, q, p)) {};
    NTRU(const unsigned int N, const unsigned int q, const unsigned int p, const unsigned int a) : f(Polynomial(N, q, p, a, a - 1, true)), g(Polynomial(N, q, p, a, a, false)) {};
    NTRU(const unsigned int N, const unsigned int q, const unsigned int p, const unsigned int a, const unsigned int b)
        : f(Polynomial(N, q, p, a, b, true)), g(Polynomial(N, q, p, a, b, false)) {};
    NTRU(Polynomial f, Polynomial g) : f(f), g(g) {};
    NTRU(std::vector<int> f, std::vector<int> g);
    unsigned get_q() { return q; };
    Polynomial get_sk_f() { return f; };
    Polynomial get_sk_g() { return g; };
    Polynomial get_pk_h() { return h; };


private:
    friend class Polynomial;
    const static unsigned int const N = 3;
    const static unsigned int const p = 3;//ͨ��Ϊ3,��������ĳ����Ա�ı����������������ʼ��
    const static unsigned int const q = 100;//��������ĳ����Ա��statec����
    const static unsigned int const a = 1;
    const static unsigned int const b = 1;
    Polynomial f;//˽Կ
    Polynomial g;//˽Կ
    Polynomial h;//��Կ
    Polynomial f_q;//f mod q����
    Polynomial f_p;// f mod p ����
    
};


#endif

