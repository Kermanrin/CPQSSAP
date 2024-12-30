#include "NTRU.h"
#include "Polynomial.h"



Polynomial NTRU::KeyGen_f(){
	std::cout << "f generated!" << std::endl;
	return Polynomial(N, q, p, cntf_1, cntf_neg1, true);
}

Polynomial NTRU::KeyGen_g(){
	std::cout << "g generated!" << std::endl;
	return Polynomial(N, q, p, cntg_1, cntg_neg1, false);

}



Polynomial NTRU::KeyGen_f_q() {
	std::vector<int> vec = f.CalVec_Inverse_Of_Polynomial(q);
	std::cout << "f_q generated!" << std::endl;
	std::cout << vec.size() << std::endl;
	int d = vec.size();

	return Polynomial(N, q, p, vec);
}


Polynomial NTRU::KeyGen_f_p() {
	std::cout << "f_p generated!" << std::endl;
	std::vector<int> vec = f.CalVec_Inverse_Of_Polynomial(p);
	int d = vec.size();

	return Polynomial(N, q, p, vec);
}

Polynomial NTRU::KeyGen_h(){
	std::vector<int> vec =  f_q.Pol_mul_Vec(g, q);
	std::cout << "h generated!" << std::endl;
	std::cout << std::endl;
	return Polynomial(N, q, p, vec);
}
//Polynomial NTRU::KeyGen_f_mod(int mod) {
//	return f.Cal_Inverse_Of_Polynomial(mod);
//}
