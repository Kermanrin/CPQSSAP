#pragma once
#ifndef NTRU_H
#define NTRU_H

#include "Polynomial.h"
#include <iostream>
#include <unordered_map>

class NTRU {
public:
    
    NTRU(std::vector<int> vec, int sl) : N(security_level_table[sl % 4][0]), p(security_level_table[sl % 4][1]), q(security_level_table[sl % 4][2]), cntf_1(security_level_table[sl % 4][3])
        , cntf_neg1(security_level_table[sl % 4][4]), cntg_1(security_level_table[sl % 4][5]), cntg_neg1(security_level_table[sl % 4][6]), h(KeyGen_h()), f(KeyGen_f()), g(KeyGen_g())
        , f_q(KeyGen_f_q()), f_p(KeyGen_f_p()) {
        int len = vec.size();
        for (int i = 0; i < len; i++) {
            vec[i] = vec[i] % p - 1;
        }
        m = Polynomial(N, q, p, vec);
    }
    Polynomial KeyGen_f();//����˽Կf
    Polynomial KeyGen_g();//����˽Կg
    Polynomial KeyGen_h();//���ɹ�Կh,��Ҫ����Polynomial��������ʽ�溯��
    Polynomial KeyGen_f_q();
    Polynomial KeyGen_f_p();
    Polynomial KeyGen_f_mod(int mod);//�������ʽģmod���棬��Ҫ����Polynomial��������ʽ�溯������ûʵ��

    Polynomial get_sk_f() { 
        printf("The sk_f is : ");
        return f; 
    };
    Polynomial get_sk_g() {
        printf("The sk_g is : ");
        return g; 
    };
    Polynomial get_pk_h() {
        printf("The pk_h is : ");
        return h; 
    };
    Polynomial get_sk_f_q() {
        printf("The sk_f_q is : ");
        return f_q; 
    };
    Polynomial get_sk_f_p() {
        printf("The sk_f_p is : ");
        return f_p; 
    };

    Polynomial get_message() {
        printf("The plaintext is : ");
        return m;
    };


private:
    friend class Polynomial;
    int security_level_table[4][7] = { {3, 3, 100, 2, 0, 1, 1},{251, 3, 257, 72, 71, 72, 72},{768, 3, 10000001, 10000, 9999}, {1024, 3, 10000001, 10000, 9999} };
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
    Polynomial f_q;//f mod q����
    Polynomial f_p;// f mod p ����
    Polynomial h;//��Կ
    Polynomial m;
    
};


#endif

