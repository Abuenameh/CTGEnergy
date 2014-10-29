#include <complex>
#include <iostream>

using namespace std;

template<class T>
complex<T> operator~(const complex<T> a) {
	return conj(a);
}

#define L 50
#define nmax 5
#define dim (nmax+1)

inline int mod(int i) {
    return (i + L) % L;
}

/*
complex<double> energyi(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyi.incl"
    return Ec;
}
complex<double> denergyi0(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyi0.incl"
    return dEc;
}
complex<double> denergyi1(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyi1.incl"
    return dEc;
}
complex<double> denergyi2(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyi2.incl"
    return dEc;
}
complex<double> denergyi3(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyi3.incl"
    return dEc;
}
complex<double> denergyi4(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyi4.incl"
    return dEc;
}
complex<double> denergyi5(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyi5.incl"
    return dEc;
}

complex<double> energyj1(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj1.incl"
    return Ec;
}
complex<double> denergyj10(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyj10.incl"
    return dEc;
}
complex<double> denergyj11(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyj11.incl"
    return dEc;
}
complex<double> denergyj12(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyj12.incl"
    return dEc;
}
complex<double> denergyj3(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyj13.incl"
    return dEc;
}
complex<double> denergyj4(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyj14.incl"
    return dEc;
}
complex<double> denergyj5(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> dEc;
#include "denergyj15.incl"
    return dEc;
}

complex<double> energyj2(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj2.incl"
    return Ec;
}

complex<double> energyj1j2_1(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj1j2_1.incl"
    return Ec;
}
complex<double> energyj1j2_2(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj1j2_2.incl"
    return Ec;
}
complex<double> energyj1j2_3(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj1j2_3.incl"
    return Ec;
}
complex<double> energyj1j2_4(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj1j2_4.incl"
    return Ec;
}
complex<double> energyj1k1_1(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj1k1_1.incl"
    return Ec;
}
complex<double> energyj1k1_2(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj1k1_2.incl"
    return Ec;
}
complex<double> energyj1k1_3(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj1k1_3.incl"
    return Ec;
}
complex<double> energyj1k1_4(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj1k1_4.incl"
    return Ec;
}
complex<double> energyj2k2_1(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj2k2_1.incl"
    return Ec;
}
complex<double> energyj2k2_2(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj2k2_2.incl"
    return Ec;
}
complex<double> energyj2k2_3(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj2k2_3.incl"
    return Ec;
}
complex<double> energyj2k2_4(int k1, int j1, int i, int j2, int k2, complex<double>** f, double* U, double U0, double* dU, double* J, double mu, complex<double> expth, complex<double> expmth, complex<double> exp2th, complex<double> expm2th) {
    complex<double> Ec;
#include "energyj2k2_4.incl"
    return Ec;
}
*/

#include "energyi.incl"
#include "energyj1.incl"
#include "energyj2.incl"
#include "energyj1j2_1.incl"
#include "energyj1j2_2_1.incl"
#include "energyj1j2_2_2.incl"
#include "energyj1j2_3_1.incl"
#include "energyj1j2_3_2.incl"
#include "energyj1j2_4.incl"
#include "energyj1k1_1.incl"
#include "energyj1k1_2.incl"
#include "energyj1k1_3.incl"
#include "energyj1k1_4.incl"
#include "energyj2k2_1.incl"
#include "energyj2k2_2.incl"
#include "energyj2k2_3.incl"
#include "energyj2k2_4.incl"

#include "denergyi0.incl"
#include "denergyi1.incl"
#include "denergyi2.incl"
#include "denergyi3.incl"
#include "denergyi4.incl"
#include "denergyi5.incl"

#include "denergyj10.incl"
#include "denergyj11.incl"
#include "denergyj12.incl"
#include "denergyj13.incl"
#include "denergyj14.incl"
#include "denergyj15.incl"

#include "denergyj20.incl"
#include "denergyj21.incl"
#include "denergyj22.incl"
#include "denergyj23.incl"
#include "denergyj24.incl"
#include "denergyj25.incl"

#include "denergyj1j2_10.incl"
#include "denergyj1j2_11.incl"
#include "denergyj1j2_12.incl"
#include "denergyj1j2_13.incl"
#include "denergyj1j2_14.incl"
#include "denergyj1j2_15.incl"

//#include "denergyj1j2_20.incl"
//#include "denergyj1j2_21.incl"
//#include "denergyj1j2_22.incl"
//#include "denergyj1j2_23.incl"
//#include "denergyj1j2_24.incl"
//#include "denergyj1j2_25.incl"
#include "denergyj1j2_2_10.incl"
#include "denergyj1j2_2_11.incl"
#include "denergyj1j2_2_12.incl"
#include "denergyj1j2_2_13.incl"
#include "denergyj1j2_2_14.incl"
#include "denergyj1j2_2_15.incl"
#include "denergyj1j2_2_20.incl"
#include "denergyj1j2_2_21.incl"
#include "denergyj1j2_2_22.incl"
#include "denergyj1j2_2_23.incl"
#include "denergyj1j2_2_24.incl"
#include "denergyj1j2_2_25.incl"

//#include "denergyj1j2_30.incl"
//#include "denergyj1j2_31.incl"
//#include "denergyj1j2_32.incl"
//#include "denergyj1j2_33.incl"
//#include "denergyj1j2_34.incl"
//#include "denergyj1j2_35.incl"
#include "denergyj1j2_3_10.incl"
#include "denergyj1j2_3_11.incl"
#include "denergyj1j2_3_12.incl"
#include "denergyj1j2_3_13.incl"
#include "denergyj1j2_3_14.incl"
#include "denergyj1j2_3_15.incl"
#include "denergyj1j2_3_20.incl"
#include "denergyj1j2_3_21.incl"
#include "denergyj1j2_3_22.incl"
#include "denergyj1j2_3_23.incl"
#include "denergyj1j2_3_24.incl"
#include "denergyj1j2_3_25.incl"

#include "denergyj1j2_40.incl"
#include "denergyj1j2_41.incl"
#include "denergyj1j2_42.incl"
#include "denergyj1j2_43.incl"
#include "denergyj1j2_44.incl"
#include "denergyj1j2_45.incl"
//#include "denergyj1j2_4_10.incl"
//#include "denergyj1j2_4_11.incl"
//#include "denergyj1j2_4_12.incl"
//#include "denergyj1j2_4_13.incl"
//#include "denergyj1j2_4_14.incl"
//#include "denergyj1j2_4_15.incl"
//#include "denergyj1j2_4_20.incl"
//#include "denergyj1j2_4_21.incl"
//#include "denergyj1j2_4_22.incl"
//#include "denergyj1j2_4_23.incl"
//#include "denergyj1j2_4_24.incl"
//#include "denergyj1j2_4_25.incl"

#include "denergyj1k1_10.incl"
#include "denergyj1k1_11.incl"
#include "denergyj1k1_12.incl"
#include "denergyj1k1_13.incl"
#include "denergyj1k1_14.incl"
#include "denergyj1k1_15.incl"

#include "denergyj1k1_20.incl"
#include "denergyj1k1_21.incl"
#include "denergyj1k1_22.incl"
#include "denergyj1k1_23.incl"
#include "denergyj1k1_24.incl"
#include "denergyj1k1_25.incl"

#include "denergyj1k1_30.incl"
#include "denergyj1k1_31.incl"
#include "denergyj1k1_32.incl"
#include "denergyj1k1_33.incl"
#include "denergyj1k1_34.incl"
#include "denergyj1k1_35.incl"

#include "denergyj1k1_40.incl"
#include "denergyj1k1_41.incl"
#include "denergyj1k1_42.incl"
#include "denergyj1k1_43.incl"
#include "denergyj1k1_44.incl"
#include "denergyj1k1_45.incl"

#include "denergyj2k2_10.incl"
#include "denergyj2k2_11.incl"
#include "denergyj2k2_12.incl"
#include "denergyj2k2_13.incl"
#include "denergyj2k2_14.incl"
#include "denergyj2k2_15.incl"

#include "denergyj2k2_20.incl"
#include "denergyj2k2_21.incl"
#include "denergyj2k2_22.incl"
#include "denergyj2k2_23.incl"
#include "denergyj2k2_24.incl"
#include "denergyj2k2_25.incl"

#include "denergyj2k2_30.incl"
#include "denergyj2k2_31.incl"
#include "denergyj2k2_32.incl"
#include "denergyj2k2_33.incl"
#include "denergyj2k2_34.incl"
#include "denergyj2k2_35.incl"

#include "denergyj2k2_40.incl"
#include "denergyj2k2_41.incl"
#include "denergyj2k2_42.incl"
#include "denergyj2k2_43.incl"
#include "denergyj2k2_44.incl"
#include "denergyj2k2_45.incl"

#define ENERGY(ei) Ec += energy##ei(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);

double energy(double* x, double* U, double U0, double* dU, double* J, double mu, double theta) {
    complex<double> expth = exp(complex<double>(0, 1)*theta);
    complex<double> expmth = conj(expth);
    complex<double> exp2th = expth*expth;
    complex<double> expm2th = conj(exp2th);
    
    complex<double> * f[L];
    for (int i = 0; i < L; i++) {
        f[i] = reinterpret_cast<complex<double>*> (&x[2 * i * dim]);
    }

    complex<double> Ec = 0;
    
    for (int i = 0; i < L; i++) {
        int k1 = mod(i - 2);
        int j1 = mod(i - 1);
        int j2 = mod(i + 1);
        int k2 = mod(i + 2);
        
        ENERGY(i)
        ENERGY(j1)
        ENERGY(j2)
        ENERGY(j1j2_1)
        ENERGY(j1j2_2_1)
        ENERGY(j1j2_2_2)
        ENERGY(j1j2_3_1)
        ENERGY(j1j2_3_2)
        ENERGY(j1j2_4)
        ENERGY(j1k1_1)
        ENERGY(j1k1_2)
        ENERGY(j1k1_3)
        ENERGY(j1k1_4)
        ENERGY(j2k2_1)
        ENERGY(j2k2_2)
        ENERGY(j2k2_3)
        ENERGY(j2k2_4)
//        Ec += energyi(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj1(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj2(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj1j2_1(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj1j2_2(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj1j2_3(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj1j2_4(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj1k1_1(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj1k1_2(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj1k1_3(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj1k1_4(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj2k2_1(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj2k2_2(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj2k2_3(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        Ec += energyj2k2_4(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);

    }
    
//    cout << Ec << endl;
    
    return real(Ec);
}

//#define DENERGY(Ei, n) cout << "About to dE" #Ei #n << endl; \
//dEc = denergy##Ei##n(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th); \
//        grad[2*(i*dim+n)] += 2*dEc.real(); \
//        grad[2*(i*dim+n)+1] += 2*dEc.imag(); \
//        cout << "Did dE" #Ei #n << endl;

#define DENERGY(Ei, n) dEc = denergy##Ei##n(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th); \
        grad[2*(i*dim+n)] += 2*dEc.real(); \
        grad[2*(i*dim+n)+1] += 2*dEc.imag();

void energygrad(double* x, double* U, double U0, double* dU, double* J, double mu, double theta, double* grad) {
//    cout << "energygrad " << grad << endl;
    complex<double> expth = exp(complex<double>(0, 1) * theta);
    complex<double> expmth = conj(expth);
    complex<double> exp2th = expth*expth;
    complex<double> expm2th = conj(exp2th);

    complex<double> * f[L];
    for (int i = 0; i < L; i++) {
        f[i] = reinterpret_cast<complex<double>*> (&x[2 * i * dim]);
    }

    complex<double> Ec = 0;
//    vector<vector<complex<double> > > dEc(L, vector<complex<double> >(dim, 0));
    complex<double> dEc = 0;
    
    fill(grad, grad+2*L*dim, 0);

    for (int i = 0; i < L; i++) {
        int k1 = mod(i - 2);
        int j1 = mod(i - 1);
        int j2 = mod(i + 1);
        int k2 = mod(i + 2);

        DENERGY(i, 0);
        DENERGY(i, 1);
        DENERGY(i, 2);
        DENERGY(i, 3);
        DENERGY(i, 4);
        DENERGY(i, 5);

        DENERGY(j1, 0);
        DENERGY(j1, 1);
        DENERGY(j1, 2);
        DENERGY(j1, 3);
        DENERGY(j1, 4);
        DENERGY(j1, 5);
        
        DENERGY(j2, 0);
        DENERGY(j2, 1);
        DENERGY(j2, 2);
        DENERGY(j2, 3);
        DENERGY(j2, 4);
        DENERGY(j2, 5);

        DENERGY(j1j2_1, 0);
        DENERGY(j1j2_1, 1);
        DENERGY(j1j2_1, 2);
        DENERGY(j1j2_1, 3);
        DENERGY(j1j2_1, 4);
        DENERGY(j1j2_1, 5);

        DENERGY(j1j2_2_1, 0);
        DENERGY(j1j2_2_2, 0);
        DENERGY(j1j2_2_1, 1);
        DENERGY(j1j2_2_2, 1);
        DENERGY(j1j2_2_1, 2);
        DENERGY(j1j2_2_2, 2);
        DENERGY(j1j2_2_1, 3);
        DENERGY(j1j2_2_2, 3);
        DENERGY(j1j2_2_1, 4);
        DENERGY(j1j2_2_2, 4);
        DENERGY(j1j2_2_1, 5);
        DENERGY(j1j2_2_2, 5);
        
        DENERGY(j1j2_3_1, 0);
        DENERGY(j1j2_3_2, 0);
        DENERGY(j1j2_3_1, 1);
        DENERGY(j1j2_3_2, 1);
        DENERGY(j1j2_3_1, 2);
        DENERGY(j1j2_3_2, 2);
        DENERGY(j1j2_3_1, 3);
        DENERGY(j1j2_3_2, 3);
        DENERGY(j1j2_3_1, 4);
        DENERGY(j1j2_3_2, 4);
        DENERGY(j1j2_3_1, 5);
        DENERGY(j1j2_3_2, 5);
        
//        DENERGY(j1j2_2, 0);
//        DENERGY(j1j2_2, 1);
//        DENERGY(j1j2_2, 2);
//        DENERGY(j1j2_2, 3);
//        DENERGY(j1j2_2, 4);
//        DENERGY(j1j2_2, 5);

//        DENERGY(j1j2_3, 0);
//        DENERGY(j1j2_3, 1);
//        DENERGY(j1j2_3, 2);
//        DENERGY(j1j2_3, 3);
//        DENERGY(j1j2_3, 4);
//        DENERGY(j1j2_3, 5);

        DENERGY(j1j2_4, 0);
        DENERGY(j1j2_4, 1);
        DENERGY(j1j2_4, 2);
        DENERGY(j1j2_4, 3);
        DENERGY(j1j2_4, 4);
        DENERGY(j1j2_4, 5);

        DENERGY(j1k1_1, 0);
        DENERGY(j1k1_1, 1);
        DENERGY(j1k1_1, 2);
        DENERGY(j1k1_1, 3);
        DENERGY(j1k1_1, 4);
        DENERGY(j1k1_1, 5);

        DENERGY(j1k1_2, 0);
        DENERGY(j1k1_2, 1);
        DENERGY(j1k1_2, 2);
        DENERGY(j1k1_2, 3);
        DENERGY(j1k1_2, 4);
        DENERGY(j1k1_2, 5);

        DENERGY(j1k1_3, 0);
        DENERGY(j1k1_3, 1);
        DENERGY(j1k1_3, 2);
        DENERGY(j1k1_3, 3);
        DENERGY(j1k1_3, 4);
        DENERGY(j1k1_3, 5);

        DENERGY(j1k1_4, 0);
        DENERGY(j1k1_4, 1);
        DENERGY(j1k1_4, 2);
        DENERGY(j1k1_4, 3);
        DENERGY(j1k1_4, 4);
        DENERGY(j1k1_4, 5);

        DENERGY(j2k2_1, 0);
        DENERGY(j2k2_1, 1);
        DENERGY(j2k2_1, 2);
        DENERGY(j2k2_1, 3);
        DENERGY(j2k2_1, 4);
        DENERGY(j2k2_1, 5);

        DENERGY(j2k2_2, 0);
        DENERGY(j2k2_2, 1);
        DENERGY(j2k2_2, 2);
        DENERGY(j2k2_2, 3);
        DENERGY(j2k2_2, 4);
        DENERGY(j2k2_2, 5);

        DENERGY(j2k2_3, 0);
        DENERGY(j2k2_3, 1);
        DENERGY(j2k2_3, 2);
        DENERGY(j2k2_3, 3);
        DENERGY(j2k2_3, 4);
        DENERGY(j2k2_3, 5);

        DENERGY(j2k2_4, 0);
        DENERGY(j2k2_4, 1);
        DENERGY(j2k2_4, 2);
        DENERGY(j2k2_4, 3);
        DENERGY(j2k2_4, 4);
        DENERGY(j2k2_4, 5);

        //        dEc = denergy0(k1, j1, i, j2, k2, f, U, U0, dU, J, mu, expth, expmth, exp2th, expm2th);
//        grad[2*(i*dim+0)] += 2*dEc.real();
//        grad[2*(i*dim+0)+1] += 2*dEc.imag();

//        dEc = denergyi0(i, f, U, U0, dU, J, mu, theta);
//        grad[2*(i*dim+0)] = 2*dEc.real();
//        grad[2*(i*dim+0)+1] = 2*dEc.imag();
//        dEc = denergy1(i, x, U, U0, dU, J, mu, theta);
//        grad[2*(i*dim+1)] = 2*dEc.real();
//        grad[2*(i*dim+1)+1] = 2*dEc.imag();
//        dEc = denergy2(i, x, U, U0, dU, J, mu, theta);
//        grad[2*(i*dim+2)] = 2*dEc.real();
//        grad[2*(i*dim+2)+1] = 2*dEc.imag();
//        dEc = denergy3(i, x, U, U0, dU, J, mu, theta);
//        grad[2*(i*dim+3)] = 2*dEc.real();
//        grad[2*(i*dim+3)+1] = 2*dEc.imag();
//        dEc = denergy4(i, x, U, U0, dU, J, mu, theta);
//        grad[2*(i*dim+4)] = 2*dEc.real();
//        grad[2*(i*dim+4)+1] = 2*dEc.imag();
//        dEc = denergy5(i, x, U, U0, dU, J, mu, theta);
//        grad[2*(i*dim+5)] = 2*dEc.real();
//        grad[2*(i*dim+5)+1] = 2*dEc.imag();

//        dEc[i][0] = denergy0(i, x, U, U0, dU, J, mu, theta);
//        dEc[i][1] = denergy1(i, x, U, U0, dU, J, mu, theta);
//        dEc[i][2] = denergy2(i, x, U, U0, dU, J, mu, theta);
//        dEc[i][3] = denergy3(i, x, U, U0, dU, J, mu, theta);
//        dEc[i][4] = denergy4(i, x, U, U0, dU, J, mu, theta);
//        dEc[i][5] = denergy5(i, x, U, U0, dU, J, mu, theta);
//
    }

//    return real(Ec);
}
