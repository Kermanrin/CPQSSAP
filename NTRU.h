#pragma once
#ifndef NTRU_H
#define NTRU_H

#include "Polynomial.h"

class NTRU {
public:
    NTRU(Polynomial f, Polynomial g);
    NTRU(std::vector<int> f, std::vector<int> g);
    unsigned get_q() { return q; };
    Polynomial get_sk_f() { 
        return f; 
    };
    Polynomial get_sk_g() { return g; };
    Polynomial get_pk_h() { return h; };


private:
    friend class Polynomial;
    const static unsigned int const N = 3;
    const static unsigned int const p = 3;//ͨ��Ϊ3,��������ĳ����Ա�ı����������������ʼ��
    const static unsigned int const q = 100;//��������ĳ����Ա��statec����
    Polynomial f;//˽Կ
    Polynomial g;//˽Կ
    Polynomial h;//��Կ
    Polynomial f_q;//f mod q����
    Polynomial f_p;// f mod p ����
};


#endif

