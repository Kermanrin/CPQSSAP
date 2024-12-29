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
    printf("������0��Ϊ���԰�ȫ������1��Ϊģ�ⰲȫ����\n");
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
    Polynomial p = ntr.get_sk_f().Cal_Inverse_Of_Polynomial(257); //���Ծ�� x^2 + x + 1 con x^ - x - 1 = -x^2 - x - 1 mod x^3 - 1
    p.print_pol();

    return 0;
}