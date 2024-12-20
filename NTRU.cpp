#include "NTRU.h"
#include "Polynomial.h"


NTRU::NTRU( Polynomial _f, Polynomial _g) : g(_g), f(_f) {


};

NTRU::NTRU(std::vector<int> f, std::vector<int> g) {
    Polynomial p1(f, NTRU::q);
    Polynomial p2(g, NTRU::q);
    //���ﲻҪд�� NTRU(N, f, g)��ݹ鹹�����캯���������ǲ�����ġ�
    this->f = p1;
    this->g = p2;
};