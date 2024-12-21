#include "NTRU.h"
#include "Polynomial.h"



NTRU::NTRU(std::vector<int> f, std::vector<int> g) {
    Polynomial p1(f, NTRU::q);
    Polynomial p2(g, NTRU::q);
    //这里不要写成 NTRU(N, f, g)会递归构调用造函数自身，这是不允许的。
    this->f = p1;
    this->g = p2;
};