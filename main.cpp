#include <iostream>
#include "Polynomial.h"
#include "NTRU.h"

using namespace std;



int main() {
    vector<int> a = { 1, 1, 1 };
    vector<int> b = { 1, -1, -1 };
    for (int i = 0; i < a.size(); i++) cout << a[i] << " ";
    cout << endl;
    int security_level;
    printf("请输入0作为测试安全参数\n");
    cin >> security_level;
    NTRU ntr(security_level);
    //cout << ntr.get_q() << endl;
    //Polynomial pp(a, 3);
    //pp.print_pol();
    cout << endl; 
    ntr.get_sk_f().print_pol();
    cout << endl;
    ntr.get_sk_g().print_pol();
    cout << endl;
    Polynomial p = ntr.get_sk_f().Pol_mul(ntr.get_sk_g()); //测试卷积 x^2 + x + 1 con x^ - x - 1 = -x^2 - x - 1 mod x^3 - 1
    p.print_pol();

    return 0;
}