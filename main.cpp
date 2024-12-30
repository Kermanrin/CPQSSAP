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
    vector<int> message(251, 1);

    printf("请输入不超过251位的消息,终止请输入114514\n");

    int input;
    int idx = 0;
    while (cin >> input, input != 114514) {
        message[idx++] = input;
    }

    printf("请输入0作为测试安全参数或1作为模拟安全参数\n");
    cin >> security_level;
    
    NTRU ntr(message, security_level);
    //cout << ntr.get_q() << endl;
    //Polynomial pp(a, 3);
    //pp.print_pol();
    cout << endl; 
    ntr.get_sk_f().print_pol();
    cout << endl;
    ntr.get_sk_g().print_pol();
    cout << endl;
    ntr.get_sk_f_p().print_pol();
    cout << endl;
    ntr.get_sk_f_q().print_pol();
    cout << endl;
    ntr.get_pk_h().print_pol();
    cout << endl;
    ntr.get_message().print_pol();
    cout << endl;
    //Polynomial p = ntr.get_sk_f().Cal_Inverse_Of_Polynomial(257); //测试卷积 x^2 + x + 1 con x^ - x - 1 = -x^2 - x - 1 mod x^3 - 1
    //p.print_pol();
    ntr.get_sk_f().Pol_mul(ntr.get_sk_f_p(), 3).print_pol();

    return 0;
}