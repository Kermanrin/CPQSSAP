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
    const static unsigned int const p = 3;//通常为3,不特属于某个成员的变量，必须在类外初始化
    const static unsigned int const q = 100;//不特属于某个成员的statec变量
    Polynomial f;//私钥
    Polynomial g;//私钥
    Polynomial h;//公钥
    Polynomial f_q;//f mod q的逆
    Polynomial f_p;// f mod p 的逆
};


#endif

