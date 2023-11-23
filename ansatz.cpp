/*
    Source file for classes of ansatz functions
*/

#include "ansatz.h"
#include <cmath>
#include <iostream>
#include <fstream>

double sqrtp(double x) //Helper-Function
{
    if(x <= 0.)
        return 0.;
    else
        return sqrt(x);
}

double oneoversqrtp(double x) //Helper-Function
{
    if(x < numeric_limits<double>::epsilon())
        return 0.;
    else
        return 1./sqrt(x);
}

double IsoAnsatz::rho(double r, double y)
{
    return g(y);
}

//Computes function g
//In isotropic case: g(z) = 2^(5/2)*Pi*int_0^z Phi(eta)*sqrt(z-eta) d(eta)
//Integral is computed using simpsons rule
double IsoAnsatz::g(double z)
{
    if(z <= 0.)
        return 0.;

    int N = 500;
    double h = z/((double)N); //stepsize
    double s = 0.; //integral
    for(int i = 0; i < N; i++)
        s += fct(((double)i)*h, 0.)*sqrtp(z-((double)i)*h) + 4.*fct(((double)i)*h + .5*h, 0.)*sqrtp(z - ((double)i)*h - .5*h) + fct(((double)i)*h + h, 0.)*sqrtp(z - ((double)i)*h - h);
    s *= 4.*sqrt(2.)*M_PI*h/6.;
    return s;
}

//Computes g'(z)
double IsoAnsatz::gprime(double z)
{
    if(z <= 0.)
        return 0.;

    int N = 500;
    double h = z/((double)N); //stepsize
    double s = 0.; //integral
    for(int i = 0; i < N; i++)
        s += fct(((double)i)*h, 0.)*oneoversqrtp(z-((double)i)*h) + 4.*fct(((double)i)*h + .5*h, 0.)*oneoversqrtp(z - ((double)i)*h - .5*h) + fct(((double)i)*h + h, 0.)*oneoversqrtp(z - ((double)i)*h - h);
    s *= 2.*sqrt(2.)*M_PI*h/6.;
    return s;
}

//Computes g''(z)
double IsoAnsatz::gprimeprime(double z)
{
    if(z <= 0.)
        return 0.;

    int N = 500;
    double h = z/((double)N); //stepsize
    double s = 0.; //integral
    for(int i = 0; i < N; i++)
        s += prime(((double)i)*h, 0.)*oneoversqrtp(z-((double)i)*h) + 4.*prime(((double)i)*h + .5*h, 0.)*oneoversqrtp(z - ((double)i)*h - .5*h) + prime(((double)i)*h + h, 0.)*oneoversqrtp(z - ((double)i)*h - h);
    s *= 2.*sqrt(2.)*M_PI*h/6.;
    return s;
}

double IsoPolyAnsatz::fct(double eta, double L)
{ //Phi(eta)
    if(eta > 0.)
        return C*pow(eta, k);
    else
        return 0.;
}

double IsoPolyAnsatz::prime(double eta, double L)
{ //Phi'(eta)
    if(eta > 0.)
        return C*k*pow(eta, k-1.);
    else
        return 0.;
}

double IsoKingAnsatz::fct(double eta, double L)
{ //Phi(eta)
    if(eta > 0.)
        return C*(exp(eta)-1.);
    else
        return 0.;
}

double IsoKingAnsatz::prime(double eta, double L)
{ //Phi'(eta)
    if(eta > 0.)
        return C*exp(eta);
    else
        return 0.;
}

double PolyAnsatz::fct(double eta, double L)
{ //Phi(eta)
    if(eta > 0. && L > L0)
        return C*pow(eta, k)*pow(L-L0, ell);
    else
        return 0.;
}
double PolyAnsatz::prime(double eta, double L)
{ //Phi'(eta)
    if(eta > 0. && L > L0)
        return C*k*pow(eta, k-1.)*pow(L-L0, ell);
    else
        return 0.;
}

double PolyAnsatz::rho(double r, double y)
{
    if(r < numeric_limits<double>::epsilon())
    {
        if(L0 > numeric_limits<double>::epsilon() || fabs(ell) > numeric_limits<double>::epsilon())
            return 0.;
        else
            return g(y);
    }

    return pow(r, 2.*ell)*g(y-L0/(2.*r*r));
}

double PolyAnsatz::g(double z)
{
    if(z < numeric_limits<double>::epsilon())
        return 0.;

    return C*ckl*pow(z, k+ell+1.5);
}

double PolyAnsatz::gprime(double z)
{
    if(k+ell+.5 < numeric_limits<double>::epsilon())
        cout << " WARNING! g' computation not save for polytropic ansatz of k=" << k << ", l=" << ell << ". ";

    if(z < numeric_limits<double>::epsilon())
        return 0.;

    return C*ckl*(k+ell+1.5)*pow(z, k+ell+.5);
}

double PolyAnsatz::gprimeprime(double z)
{
    if(k+ell-.5 < numeric_limits<double>::epsilon())
        cout << " WARNING! g'' computation not save for polytropic ansatz of k=" << k << ", l=" << ell << ". ";

    if(z < numeric_limits<double>::epsilon())
        return 0.;

    return C*ckl*(k+ell+1.5)*(k+ell+.5)*pow(z, k+ell-.5);
}
