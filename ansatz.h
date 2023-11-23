#ifndef ANSATZ_H_INCLUDED
#define ANSATZ_H_INCLUDED

/*
    Header for classes of ansatz functions for steady states
*/

#include <vector>
#include <string>
#include <iostream>
#include <limits> //To use numeric-limits (like double_epsilon)
#include <cmath>

using namespace std;

class Ansatz { //Polymorphic Parent Class for ansatz function
    public:
        double C; //Constant/Factor for ansatz function
        double L0; //Parameter L_0 of stst
        double ell; //Polytropic exponent ell

        virtual double fct(double eta, double L) = 0; //Phi(eta) * (L-L_0)^\ell
        virtual double prime(double eta, double L) = 0; //Phi'(eta) * (L-L_0)^\ell
        virtual double rho(double r, double y) = 0; //rho(r,y)
        void setConstant(double newC) {C = newC;}

        virtual double g(double z) = 0; //Function g(z) relating rho_0 and U_0
        virtual double gprime(double z) = 0; //g'(z)
        virtual double gprimeprime(double z) = 0; //g''(z)
};


class IsoAnsatz: public Ansatz { //class for isotropic ansatz functions
    public:
        double rho(double r, double y); //rho(r,y)
        double g(double z); //g(z)
        double gprime(double z); //g'(z)
        double gprimeprime(double z); //g''(z)
};


class IsoPolyAnsatz: public IsoAnsatz { //class for isotropic polytropic ansatz function, i.e., \Phi(\eta)=C*\eta_+^k
    public:
        double k; //exponent
        double fct(double eta, double L); //Phi(eta)
        double prime(double eta, double L); //Phi'(eta)

        IsoPolyAnsatz(double ki, double Ci=1.)
        { //constructor
            L0 = 0.;
            ell = 0.;

            C = Ci;
            k = ki;

            if(k < -numeric_limits<double>::epsilon())
                cout << endl << "WARNING! Polytropes with k < 0 might not work properly..." << endl << endl;

            if(C < numeric_limits<double>::epsilon())
                cout << endl << "WARNING! Constant C <= 0 does not make sense?!" << endl << endl;
        }
};

class IsoKingAnsatz: public IsoAnsatz { //class for King ansatz function, i.e., \Phi(\eta)= C * (e^\eta-1.)_+
    public:
        double fct(double eta, double L); //Phi(eta)
        double prime(double eta, double L); //Phi'(eta)

        IsoKingAnsatz(double Ci=1.)
        {
            L0 = 0.;
            ell = 0.;

            C = Ci;

            if(C < numeric_limits<double>::epsilon())
                cout << endl << "WARNING! Constant C <= 0 does not make sense?!" << endl << endl;
        }
};

class PolyAnsatz: public Ansatz { //class for polytropic ansatz function, i.e., \varphi(E,L)= C * (E_0-E)_+^k * (L-L_0)_+^l
    public:
        double k; //exponent
        double ckl; //constant

        double fct(double eta, double L);
        double prime(double eta, double L);
        double rho(double r, double y); //rho(r,y)
        double g(double z); //g(z)
        double gprime(double z); //g'(z)
        double gprimeprime(double z); //g''(z)

        PolyAnsatz(double ki, double elli, double L0i, double Ci=1.)
        { //constructor
            k = ki;
            ell = elli;
            L0 = L0i;
            C = Ci;

            double cl = pow(2., ell+1.5)*M_PI*sqrt(M_PI)*tgamma(ell+1.)/tgamma(ell+1.5);
            ckl = cl*tgamma(k+1.)*tgamma(ell+1.5)/tgamma(k+ell+2.5);

            if(k < -numeric_limits<double>::epsilon())
                cout << endl << "WARNING! Polytropes with k < 0 might not work properly..." << endl << endl;

            if(C < numeric_limits<double>::epsilon())
                cout << endl << "WARNING! Constant C <= 0 does not make sense?!" << endl << endl;

            if(ell+.5 < numeric_limits<double>::epsilon())
                cout << endl << "WARNING! Polytropes with ell <= -.5 might not work properly..." << endl << endl;

            if(L0 < -numeric_limits<double>::epsilon())
                cout << endl << "WARNING! Parameter L0 < 0 does not make sense?!" << endl << endl;
        }
};


#endif // ANSATZ_H_INCLUDED
