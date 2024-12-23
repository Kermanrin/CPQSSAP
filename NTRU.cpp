#include "NTRU.h"
#include "Polynomial.h"



Polynomial NTRU::KeyGen_f(){
	return Polynomial(N, q, p, cntf_1, cntf_neg1, true);
}

Polynomial NTRU::KeyGen_g(){
	return Polynomial(N, q, p, cntg_1, cntg_neg1, true);

}

Polynomial NTRU::KeyGen_h(){
	return f;
}

Polynomial NTRU::KeyGen_f_mod(int mod) {
	return f.Cal_Inverse_Of_Polynomial(mod);
}
