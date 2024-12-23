#pragma once
#ifndef NTRU_H
#define NTRU_H

#include "Polynomial.h"
#include <iostream>
#include <unordered_map>

class NTRU {
public:
    
    NTRU(int sl) : N(security_level_table[sl % 4][0]), p(security_level_table[sl % 4][1]), q(security_level_table[sl % 4][2]), cntf_1(security_level_table[sl % 4][3])
        , cntf_neg1(security_level_table[sl % 4][4]), cntg_1(security_level_table[sl % 4][5]), cntg_neg1(security_level_table[sl % 4][6]), f(KeyGen_f()), g(KeyGen_g()) {};

    Polynomial KeyGen_f();//����˽Կf
    Polynomial KeyGen_g();//����˽Կg
    Polynomial KeyGen_h();//���ɹ�Կh,��Ҫ����Polynomial��������ʽ�溯������ûʵ��
    Polynomial KeyGen_f_mod(int mod);//�������ʽģmod���棬��Ҫ����Polynomial��������ʽ�溯������ûʵ��

    Polynomial get_sk_f() { return f; };
    Polynomial get_sk_g() { return g; };
    Polynomial get_pk_h() { return h; };


private:
    friend class Polynomial;
    int security_level_table[4][7] = { {3, 3, 100, 2, 0, 1, 1},{512, 3, 10000001, 10000, 10000, 10000, 9999},{768, 3, 10000001, 10000, 9999}, {1024, 3, 10000001, 10000, 9999} };
    int security_level;
    unsigned int const N;
    unsigned int const p;//ͨ��Ϊ3,��������ĳ����Ա�ı����������������ʼ��
    unsigned int const q;//��������ĳ����Ա��statec����
    unsigned int const cntf_1;
    unsigned int const cntf_neg1;
    unsigned int const cntg_1;
    unsigned int const cntg_neg1;
    Polynomial f;//˽Կ
    Polynomial g;//˽Կ
    Polynomial h;//��Կ
    Polynomial f_q;//f mod q����
    Polynomial f_p;// f mod p ����
    
};


#endif

