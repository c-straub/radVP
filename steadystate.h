#ifndef STEADYSTATE_H_INCLUDED
#define STEADYSTATE_H_INCLUDED

/////////////////////////////////////////////////
/// Class for Steady State Stuff
/////////////////////////////////////////////////

#include "rad_fct.h"
#include "ansatz.h"
#include "XV.h"

#include <string>
#include <vector>


class EL{ //class for pairs (E,L)
    public:
        double E;
        double L;
        EL(double Ei, double Li) {E=Ei; L=Li;} //constructor
};

class SteadyState { //class for a steady state
    protected:
        double mtot; //total mass
        double rmax; //outer radius; more precisely, first radius of r-grid where rho vanishes after support
        double rmin; //inner radius; more precisely, in case of shell, last radius where rho still vanishes
        double e0; //cut off energy
        RadFct pot; //potential
        RadFct mass; //local mass
        RadFct rho; //massdensity
        Ansatz* phi; //ansatzfunction

        double L0; //minimal occurring L in support
        double Lmax; //maximal occurring L in support

    public:
        double getmtot() {return mtot;}
        double getrmax() {return rmax;}
        double getrmin() {return rmin;}
        double gete0() {return e0;}
        RadFct getpot() {return pot;}
        double getL0() {return L0;}
        double getLmax() {return Lmax;}
        RadFct getmass() {return mass;}
        RadFct getrho() {return rho;}
        double getdr() {return pot.getdr();}
        double getkappa() {return (e0-pot_atr(0.));}

        double potati(int ind) {return pot.atind(ind);}
        double massati(int ind) {return mass.atind(ind);}
        double rhoati(int ind) {return rho.atind(ind);}

        SteadyState(double kappa, double dr, Ansatz* phii); //computation of steady states using basic input

        double energy(double r, double w, double L); //Evaluates particle energy at (r,w,L)
        double atrwl(double r, double w, double L); //Evaluates Steady State at (r,w,L)
        double phiprimeatEL(double E, double L) //varphi'(E,L)
            {return -phi->prime(e0-E, L);}
        double phiprimeatrwl(double r, double w, double L); //Calculates varphi'(E(r,w,L), L)

        double pot_atr(double r); //Evaluates stationary potential at radius
        double mass_atr(double r); //Evaluates stationary local mass at radius
        double rho_atr(double r); //Evaluates stationary rho at radius
        double potprime_atr(double r); //Evaluates radial derivative of stationary potential at radius

        double rhoprime_atr(double r); //Evaluates stationary rho' at radius
        double rhoprimeprime_atr(double r); //Evaluates stationary rho'' at radius

        double potprimeoverr_atr(double r); //Evaluates U_0'(r)/r
        double potprimeoverrmax(int Nr); //Computes Maximum of above function

        void plotrad(string fname, int radsteps); //creates Plot file for radial data with a maximum of maxsteps radial evaluations
             //FORMAT: 1) radius |2) potential |3) rho |4) mass |5) rhofalse
        void plotrad(string fname, int radsteps, double rplotmin, double rplotmax); //Same as above, but with prescribed min and max radius

        double effpot(double L, double r); //Computes Psi_L(r)
        double effpotprime(double L, double r); //Computes Psi_L'(r)
        double effpotprimeprime(double L, double r); //Computes Psi_L''(r)
        double effpotprime3(double L, double r); //Computes Psi_L'''(r)
        double effpotprime4(double L, double r); //Computes Psi_L''''(r)

        void ploteffpot(string fname, int amountL, int radsteps, double rplotmin, double rplotmax); //Plot effpot and derivatives in file for some L-values
            //FORMAT: 1) radius |2) Psi_L(r) |3) Psi_L'(r) |4) Psi_L''(r)

        double computerminus(EL el); //Computes r_-(E,L)
        double computerplus(EL el); //Computes r_+(E,L);
        double computerl(double L); //Computes r_L
        double Emin(double L); //Computes minimum of effective potential, i.e., Psi_L(r_L)

        bool inSupport(EL el); //returns whether a given (E,L) pair is inside the steady state's support
        bool inSupport(EL el, double EminL); //Same as above, but with given EminL
        bool atSupportboundary(EL el); //return whether a given (E,L) lies at boundary of steady state's support

        double setL0(double dr, double kappa); //Calculates minimal occurring L value insed steady state support
        double setLmax(); //Calculates maximal occurring L value inside steady state support
        double getLmax(double E); //Returns maximal L value for some E

        XV charSysRHS(XV xv); //Returns right-hand side of characteristic system (of stst) evaluated at xv

        double computeT(EL el, double dt); //Computes T(E,L) for given pair (E,L)
        double computeTCorner(double dt); //Computes T(E_0,L_0)
        vector<double> computeT(vector<EL> els, double dt); //Computes T(E,L) for every element in given vector; parallelised

        pair<double,double> trajectoriesolver(double rinit, double winit, double L, double dt,
                                              bool stopatradius, double stoprad, bool stopaftertime, double stoptime, bool stopafterperiod,
                                              bool compInt, double *integral, RadFct *Uprime); //function to solve the characteristic system of the steady state

        pair<double,double> trajectoriesolver(double rinit, double winit, double L, double dt, double T); //Shortcut for above function; stop after time T and dont compute integral
        pair<double,double> trajectoriesolver(double rinit, double winit, double L, double dt, double T,
                                              double *integral, RadFct *Uprime); //Shortcut for above function; stop after time and compute integral

        vector<EL> createELtriangle(int Esteps, int Lsteps); //Returns vector of (E,L) values which create the (E,L)-triangle, i.e., the (E,L)-support of the solution
        vector<EL> createELtriangle(int Esteps, int Lsteps, int minELs); //Same as above, but makes sure that we get at least minELs EL-pairs by iteratively increasing Esteps and Lsteps

        pair<double,double> computeTminTmax(double dt, int Esteps, int Lsteps); //Search for Min and Max of Period Function on stst support on given grid
        pair<double,double> computeTminTmax(double dt, int Esteps, int Lsteps, int minELs); //Same as above but takes at least minELs (E,L)-pairs
        pair<double,double> computeTminTmax(double dt, vector<double> Tvals); //Same as above but with already given T-vals (however, add some extreme values)

        bool TmaxatCorner(double dt, double Tmax, double EPS); //Returns whether Tmax value is attained at corner (E_0,L_0) (up to error of EPS)
        bool Tmaxatboundary(vector<EL> els, vector<double> Ts, double EPS); //Return whether Maximum of T is attained at boundary points of (E,L)-triangle

        bool TincE(vector<double> Ts, vector<EL> ELs, double EPS); //Returns whether T is monotone in E (up to error of EPS)
        bool TdecL(vector<double> Ts, vector<EL> ELs, double EPS); //Returns whether T is monotone in E (up to error of EPS)

        vector<bool> Tmaxmons(double dt, int Esteps, int Lsteps, double EPS); //Does all of three functions together and returns results
            //FORMAT: 1) Tmax at Corner? |2) Tmax at Boundary? |3) T increasing in E? |4) T decreasing in L?
        vector<bool> Tmaxmons(double dt, int Esteps, int Lsteps, int minELs, double EPS); //Same as above with minimum number of (E,L)s
        vector<bool> Tmaxmons(double dt, vector<EL> ELs, double EPS); //Same as above with prescribed (E,L)'s
        vector<bool> Tmaxmons(double dt, vector<EL> ELs, vector<double> Ts, double EPS); //Same as above with prescribed Tvals

        void plotTs(string fname, vector<EL> els, double dt); //Computes T-values and plot them in file
        void plotTs(string fname, vector<EL> els, vector<double> Ts, double dt); //Same as above but T-values already computed
        void plotTs(string fname, int Esteps, int Lsteps, double dt); //Computes T-values and plot them in file
        void plotTs(string fname, int Esteps, int Lsteps, int minELs, double dt); //Same as above but takes at least minELs (E,L)-pairs

        void plotTs_fixedL(string fname, int amountL, int Esteps, double dt); //Computes T-values for plots of T( . ,L) for several L's
        void plotTs_fixedE(string fname, int amountE, int Lsteps, double dt); //Computes T-values for plots of T(E, . ) for several E's

        double GL(double L, double r); //Computes G_L(r)
        double GL(double L, double r, double EminL); //Computes G_L(r), but with given EminL
        double HL(double L, double r); //Computes H_L(r)
        double HL(double L, double r, double EminL); //Computes H_L(r), but with given EminL
        double radreflection(double L, double r); //For given r and L, return the radius s\neq r such that Psi_L(r)=Psi_L(s) ; this function is sometimes called zeta_L(r), a radial reflection mapping
        double GLref(double L, double r); //Computes G_L^{ref}(r)
        double GLref(double L, double r, double rL, double EminL); //Computes G_L^{ref}(r), but with given rL and EminL
        double HLref(double L, double r); //Computes H_L^{ref}(r)
        double HLref(double L, double r, double rL, double EminL); //Computes H_L^{ref}(r), but with given rL and EminL

        bool GLrefHLrefPos(int amountL, double dr, double drEPS, double valEPS); //Returns whether functions GLref and HLref are always positive; parallelised

        void plotGLHL(string fname, int amountL, int radsteps); //Prints plotable data for GL, HL, GLref, and HLref for different L's into file
            //FORMAT: 1) radius |2) G_L(r) |3) H_L(r) |4) G_Lref(r) |5) H_Lref(r)

        double dETnearcirc(double L); //Returns limit of \partial_E T(E,LL) as (E,LL) -> (Emin L, L)
        double dLTnearcirc(double L); //Returns limit of \partial_L T(E,LL) as (E,LL) -> (Emin L, L)

        double kunzerhs(double Tmax); //Returns right-hand side of Kunze Crit, i.e., 4*M_PI*M_PI/(Tmax*Tmax)

        SteadyState rescaleConstRadius(); //Returns Rescale steady state by multiplying phi by a constant s.t. new stst has Rmax=1
};

SteadyState ststfromfile(string fname); //Returns steady state based on input file
    // FORMAT: k ell L0 kappa
    //         drst
    // If k<0, we do a King Model. If kappa<=0 (and k>=0, L0=0), we choose kappa s.t. Rmax=1

void ststfamily_kappa(string fname, Ansatz* phi, vector<double> kappas, double dr); //Compute whole stst family for fixed ansatz and variable kappa; save global data in file
    //FORMAT: 1) kappa |2) M |3) Rmin |4) Rmax |5) E_0 |6) L_0 |7) L_max |8) kappa
void ststfamily_ks(string fname, vector<double> ks, double dr); //Compute iso poly stst family with different ks and kappa=1; save global data in file
    //FORMAT: 1) k |2) M |3) Rmin |4) Rmax |5) E_0 |6) L_0 |7) L_max |8) kappa
void ststfamily_ks_radnorm(string fname, vector<double> ks, double dr); //Same as above but choose kappa s.t. Rmax=1
    //FORMAT: 1) k |2) M |3) Rmin |4) Rmax |5) E_0 |6) L_0 |7) L_max |8) kappa

SteadyState rescaleConstRadius(); //Returns Rescale steady state by multiplying phi by a constant s.t. new stst has Rmax=1
double kappa4r1(double k, double ell, double dr); //Returns the kappa value which we need for Rmax=1 in case of polytropic ansatz with L0=0
vector<double> kappas4r1(vector<double> ks, double ell, double dr); //Same as above, but parallelised for whole family of iso polys
SteadyState polykappascaled(double k, double ell, double dr); //Return Poly stst with kappa chosen s.t. Rmax=1

void ststfamily_Ts_kappa(string fname, Ansatz* phi, vector<double> kappas, double dr,
                         int Esteps, int Lsteps, int minELs, double dt, bool seperatePlotFiles,
                         double EPST,
                         bool withGHref, int amountL, double drrefs, double drrefsEPS, double valrefsEPS);
                         //Compute stst family for fixed ansatz and calculate Tmin Tmax values & analyse monotonicity; eventually containing GLHLref analysis
    //FORMAT: 1) kappa |2) Tmin |3) Tmax |4) Tmax/Tmin
    //       |5) Tmax at (E_0,L_0)? |6) Tmax on Boundary? |7) T increasing in E? |8) T decreasing in L? |9) GLref>0 and HLref>0 ?
    //       |10) Kunze Crit LHS |11) Kunze Crit RHS |12) Kunze Crit holds?
void ststfamily_Ts_ks(string fname, vector<double> ks, double dr,
                      int Esteps, int Lsteps, int minELs, double dt, bool seperatePlotFiles,
                      double EPST,
                      bool withGHref, int amountL, double drrefs, double drrefsEPS, double valrefsEPS);
                      //Computes iso poly stst family with different ks and kappa=1 and calculate same as above
    //FORMAT: Similar as above, but |1) k |...
void ststfamily_Ts_ks_radnorm(string fname, vector<double> ks, double dr,
                      int Esteps, int Lsteps, int minELs, double dt, bool seperatePlotFiles,
                      double EPST,
                      bool withGHref, int amountL, double drrefs, double drrefsEPS, double valrefsEPS);
                      //Same as above but choose kappa s.t. Rmax=1

void clearOutputFile(string fname, string sprefix);
void appendSSToOutputFile(string fname, string sprefix, stringstream *ss);

#endif // STEADYSTATE_H_INCLUDED
