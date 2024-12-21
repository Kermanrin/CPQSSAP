#include <iostream>
#include "Polynomial.h"
#include "NTRU.h"

using namespace std;



int main() {
    vector<int> a = { 1, 1, 1 };
    vector<int> b = { 1, -1, -1 };
    for (int i = 0; i < a.size(); i++) cout << a[i] << " ";
    cout << endl;
    NTRU ntr(3, 10, 3, 1, 1);
    //cout << ntr.get_q() << endl;
    //Polynomial pp(a, 3);
    //pp.print_pol();
    cout << endl; 
    ntr.get_sk_f().print_pol();
    cout << endl;
    ntr.get_sk_g().print_pol();
    cout << endl;
    Polynomial p = ntr.get_sk_f().Pol_mul(  ntr.get_sk_g()); //²âÊÔ¾í»ý x^2 + x + 1 con x^ - x - 1 = -x^2 - x - 1 mod x^3 - 1
    p.print_pol();

    return 0;
}