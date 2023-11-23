#define _USE_MATH_DEFINES //To use PI
#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream> //For the output in txt-files via ofstream etc
#include <vector> //Modern way of handling arrays in c++
#include <sstream> //For string streams
#include <limits> //To use numeric-limits (like double_epsilon)
#include <chrono> //For time measurement
#include <algorithm> //For find
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string>
#include <functional>

#include "steadystate.h"




SteadyState::SteadyState(double kappa, double dr, Ansatz* phii)
{
    phi = phii;

    int kr = 0;
    double y = kappa;
    double m = .0;
    double rmaxst = 100000.; //Safes us from infinite loops
    vector<double> yarr;
    vector<double> rhoarr;
    vector<double> marr;

    if(kappa < numeric_limits<double>::epsilon())
        exit(1);

    double r = ((double)kr)*dr;
    do //Determine inner radius on dr grid
    {
        kr++;
        r = ((double)kr)*dr;

        if(r > rmaxst)
        {
            cout << endl << endl << "Fatal error; cant find matter for steady state..." << endl << endl;
            exit(6);
        }
    }
    while(phi->rho(r, kappa) < numeric_limits<double>::epsilon());

    kr--;
    rmin = ((double)kr)*dr;

    r = rmin;

    double rhoprev = phi->rho(r, y);
    double yprev = y;
    double mprev = 0.;

    double yhalf, rhalf, rhohalf, mhalf;
    double rhocurr = 0.;

    bool reachedvacuum = false;

    while(!reachedvacuum)
    {
        yprev = y;
        mprev = m;

        // First step of midpoint method
        rhoprev = phi->rho(r, y);
        if(r > numeric_limits<double>::epsilon())
            yhalf = y - .5*(m/r)*(dr/r);
        else
            yhalf = y;
        rhalf = r + .5*dr;
        mhalf = m + 4.*M_PI*r*r*.5*dr*rhoprev;

        // Second step of midpoint method
        rhohalf  = phi->rho(rhalf, yhalf);
        if(r > numeric_limits<double>::epsilon())
            y -= (mhalf/r)*(dr/r); //Compute y using mass at half grid point
        m += 4.*M_PI*rhalf*rhalf*dr*rhohalf;

        //Save data in arrays, for radial point r (before update!)
        yarr.push_back(yprev);
        rhoarr.push_back(rhoprev);
        marr.push_back(mprev);

        //Go to next radius
        kr++;
        r = ((double)kr)*dr;

        rhocurr = phi->rho(r, y);

        //Check whether (old) radius is in vacuum
        if(rhoprev < numeric_limits<double>::epsilon() && rhocurr < numeric_limits<double>::epsilon())
            reachedvacuum = true;

        if(r > rmaxst)
        {
            cout << endl << endl
                 << "Fatal Error during steady state computation; reached maximal radius r = " << r << " ..."
                 << endl << endl;
            exit(3);
        }
    }

    r = ((double)(kr-1))*dr;
    mtot = m;
    rmax = r;
    e0 = y-m/r;

    vector<double> potarr;
    for(int i = 0; i < yarr.size(); i++)
        potarr.push_back(e0 - yarr.at(i));

    pot = RadFct(dr, yarr.size(), rmin, potarr);
    mass = RadFct(dr, yarr.size(), rmin, marr);
    rho = RadFct(dr, yarr.size(), rmin, rhoarr);

    L0 = phi->L0; //setL0(dr, kappa);
    setLmax();
}

SteadyState ststfromfile(string fname)
{
    bool success = true;
    double k, ell, L0, kappa;
    double drst;
    ifstream parameterin(fname);

    if(!(parameterin >> k >> ell >> L0 >> kappa))
        success = false;
    if(!(parameterin >> drst))
        success = false;

    if(!success)
    {
        cout << endl << "Illegal input in Steady State Input File " << fname << " ..." << endl;
        exit(570401);
    }

    Ansatz *phi;
    if(k < -numeric_limits<double>::epsilon())
        phi = new IsoKingAnsatz();
    else
        phi = new PolyAnsatz(k, ell, L0);

    if(kappa > numeric_limits<double>::epsilon())
        return SteadyState(kappa, drst, phi);
    else
    {
        if(k > -numeric_limits<double>::epsilon() && fabs(L0) < numeric_limits<double>::epsilon())
            return polykappascaled(k, ell, drst);
        else
        {
            cout << endl << "Illegal input in Steady State Input File " << fname << " ..." << endl;
            exit(570402);
        }
    }
}


double SteadyState::setL0(double dr, double kappa)
{//Assumes that rmin is already computed!
    if(rmin < numeric_limits<double>::epsilon())
    {
        L0 = 0.;
        return 0.;
    }
    else
    {//Compute rmin more accurately to determine L0
        double rminl = rmin; //Left point
        double rminr = rmin+dr; //Right point

        if(phi->rho(rminr, kappa) < numeric_limits<double>::epsilon())
        {
            cout << endl << endl << "Unexpected error when computing L0...?!" << endl << endl;
            exit(4);
        }

        int counterhelper = 0;
        while(counterhelper < 10)
        {
            double rminhalf = .5*(rminl+rminr);
            if(phi->rho(rminhalf, kappa) < numeric_limits<double>::epsilon())
                rminl = rminhalf;
            else
                rminr = rminhalf;
            counterhelper++;
        }

        double rminprec = .5*(rminl+rminr); //More precise rmin

        L0 = 2.*kappa*rminprec*rminprec;
        return L0;
    }
}


double SteadyState::setLmax()
{
    double Lmaxbound = 2.*rmax*rmax*(e0-potati(0));
    if(Lmaxbound - L0 < numeric_limits<double>::epsilon())
    {
        cout << endl << endl << "Unexpected error when computing Lmax...?!" << endl << endl;
        exit(5);
    }

    double L = .5*(Lmaxbound+L0);
    double Lleft = L0; //Asumes that L0 is already set!
    double Lright = Lmaxbound;
    int Lk = 1;
    while(Lk <= 20)
    {
        if(Emin(L) < e0)
        {
            Lleft = L;
            L = .5*(Lleft + Lright);
        }
        else
        {
            Lright = L;
            L = .5*(Lleft + Lright);
        }

        Lk++;
    }

    Lmax = Lleft;
    return Lleft;
}

double SteadyState::getLmax(double E)
{
    double Lmaxbound = 2.*rmax*rmax*(E-potati(0));
    if(Lmaxbound - 0. < numeric_limits<double>::epsilon())
    {
        cout << endl << endl << "Unexpected error when computing Lmax for E=" << E << endl << endl;
        exit(5);
    }

    double L = .5*(Lmaxbound+0.);
    double Lleft = 0.;
    double Lright = Lmaxbound;
    int Lk = 1;
    while(Lk <= 20)
    {
        if(Emin(L) < E)
        {
            Lleft = L;
            L = .5*(Lleft + Lright);
        }
        else
        {
            Lright = L;
            L = .5*(Lleft + Lright);
        }

        Lk++;
    }

    return Lleft;
}



double SteadyState::energy(double r, double w, double L)
{
    return effpot(L,r) + .5*w*w;
}


double SteadyState::atrwl(double r, double w, double L)
{
    return phi->fct(e0 - energy(r,w,L), L);
}

double SteadyState::phiprimeatrwl(double r, double w, double L)
{
    return phiprimeatEL(energy(r,w,L), L);
}

double SteadyState::pot_atr(double r)
{
    if(r - pot.getrmin() < numeric_limits<double>::epsilon())
        return pot.atind(0); //potential
    else if(pot.inrange(r))
        return pot.atr(r);
    else
        return -mtot/r;
}

double SteadyState::mass_atr(double r)
{
    if(r - mass.getrmin() < numeric_limits<double>::epsilon())
        return 0.;
    else if(mass.inrange(r))
        return mass.atr(r);
    else
        return mtot;
}

double SteadyState::rho_atr(double r)
{
    if(r - rho.getrmin() < numeric_limits<double>::epsilon())
        return rho.atind(0);
    else if(rho.inrange(r))
        return rho.atr(r);
    else
        return 0.;
}

double SteadyState::potprime_atr(double r)
{
    if(r - mass.getrmin() < numeric_limits<double>::epsilon() || r < numeric_limits<double>::epsilon())
        return 0.;
    else
        return (mass_atr(r)/(r*r));
}

double SteadyState::potprimeoverr_atr(double r)
{
    if(r < numeric_limits<double>::epsilon())
    {
        if(phi->L0 > numeric_limits<double>::epsilon())
            return 0.;
        else
        {
            if(phi->ell > -numeric_limits<double>::epsilon())
                return (4.*M_PI/3.)*rho_atr(0.);
            else
                return numeric_limits<double>::max();
        }
    }

    return potprime_atr(r)/r;
}

double SteadyState::potprimeoverrmax(int Nr)
{
    double returner = numeric_limits<double>::lowest();
    for(int i = 0; i <= Nr; i++)
    {
        double valueloc = potprimeoverr_atr(rmin + ((double)i)*(rmax-rmin)/((double)Nr));
        if(valueloc > returner)
            returner = valueloc;
    }
    return returner;
}

double SteadyState::effpot(double L, double r)
{
    if(L < numeric_limits<double>::epsilon())
        return pot_atr(r);

    if(r < numeric_limits<double>::epsilon())
    {
        cout << endl << endl << "Fatal error; trying to evaluate Psi_L(r) at r<=0..." << endl << endl;
        exit(3);
    }
    else
        return (pot_atr(r) + .5*L/(r*r));
}

double SteadyState::effpotprime(double L, double r)
{
    if(L < numeric_limits<double>::epsilon())
    {
        if(r < numeric_limits<double>::epsilon())
            return 0.;
        else
            return mass_atr(r)/(r*r);
    }

    if(r < -numeric_limits<double>::epsilon())
        return numeric_limits<double>::lowest();

    if(r < numeric_limits<double>::epsilon())
    {
        cout << endl << endl << "Fatal error; trying to evaluate Psi_L'(r) at r<=0..." << endl;
        cout << setprecision(12) << "L=" << L << ", r=" << r << "." << endl;
        cout << "Steady State has Rmax=" << rmax << ", L0=" << L0 << ", Lmax=" << Lmax << endl << endl;
        exit(4);
    }

    return (mass_atr(r)/(r*r) - L/(r*r*r));
}

double SteadyState::effpotprimeprime(double L, double r)
{
    if(L < numeric_limits<double>::epsilon())
    {
        if(r < numeric_limits<double>::epsilon())
            return (4.*M_PI/3.)*rhoati(0);
        else
            return 4.*M_PI*rho_atr(r) - 2.*mass_atr(r)/(r*r*r);
    }

    if(r < numeric_limits<double>::epsilon())
    {
        cout << endl << endl << "Fatal error; trying to evaluate Psi_L''(r) at r<=0..." << endl << endl;
        exit(5);
    }

    return (4.*M_PI*rho_atr(r) - 2.*mass_atr(r)/(r*r*r) + 3.*L/(r*r*r*r));
}

double SteadyState::rhoprime_atr(double r)
{
    double psil = effpot(L0,r);
    double psilp = effpotprime(L0,r);

    if(fabs(phi->ell) < numeric_limits<double>::epsilon())
        return (-1.)*phi->gprime(e0-psil)*psilp;
    else
        return 2.*phi->ell*pow(r, 2.*phi->ell-1.)*phi->g(e0-psil) - pow(r, 2.*phi->ell)*phi->gprime(e0-psil)*psilp;
}

double SteadyState::rhoprimeprime_atr(double r)
{
    double psil = effpot(L0,r);
    double psilp = effpotprime(L0,r);
    double psilpp = effpotprimeprime(L0,r);

    if(fabs(phi->ell) < numeric_limits<double>::epsilon())
        return phi->gprimeprime(e0-psil)*psilp*psilp - phi->gprime(e0-psil)*psilpp;
    else
        return 2.*phi->ell*(2.*phi->ell-1.)*pow(r, 2.*phi->ell-2.)*phi->g(e0-psil)
                - 4.*phi->ell*pow(r, 2.*phi->ell-1.)*phi->gprime(e0-psil)*psilp
                + pow(r, 2.*phi->ell)*phi->gprimeprime(e0-psil)*psilp*psilp
                - pow(r, 2.*phi->ell)*phi->gprime(e0-psil)*psilpp;
}

double SteadyState::effpotprime3(double L, double r)
{
    if(r < numeric_limits<double>::epsilon())
    {
        if(L < numeric_limits<double>::epsilon())
            return 0.;
        else
            return numeric_limits<double>::lowest();
    }

    return 4.*M_PI*rhoprime_atr(r) - 8.*M_PI*rho_atr(r)/r + 6.*mass_atr(r)/(r*r*r*r) - 12.*L/(r*r*r*r*r);
}

double SteadyState::effpotprime4(double L, double r)
{
    if(r < numeric_limits<double>::epsilon())
    {
        if(L < numeric_limits<double>::epsilon())
            return 0.; //Might not be the correct extension; but shouldn't be relevant
        else
            return numeric_limits<double>::max();
    }

    return 4.*M_PI*rhoprimeprime_atr(r) - 8.*M_PI*rhoprime_atr(r)/r + 32.*M_PI*rho_atr(r)/(r*r)
            - 24.*mass_atr(r)/(r*r*r*r*r) - 60.*L/(r*r*r*r*r*r);
}

double SteadyState::computerl(double L)
{
    if(L < numeric_limits<double>::epsilon())
        return 0.;

    double dr = getdr();

    int kr = 1;
    double rloc = dr;

    while(true)
    {
        //rloc = ((double)kr)*dr;

        if(kr <= 0)
        {
            cout << endl << endl << "Fatal Error! Ueberlauf!!" << endl << endl;
            exit(4051);
        }

        if(effpotprime(L, rloc) > 0.)
        {
            if(kr == 1)
                return dr;

            double psilpleft = effpotprime(L,rloc-dr);
            double psilpright = effpotprime(L,rloc);

            if(psilpleft > 0.) //Shouldnt be the case
                return rloc-dr;

            double convcoeff = (0. - psilpleft)/(psilpright-psilpleft); //Interpolate effpotprime linearly

            return rloc - convcoeff*dr;
        }

        kr++;
        rloc += dr;

        if(rloc > 1000000.)
        {
            cout << endl << endl
                 << "Fatal Error when computing r_L for L = " << L << " , reached maximal radius r = " << rloc << " ..."
                 << endl << endl;
            exit(3);
        }
    }
}

double SteadyState::Emin(double L)
{
    return effpot(L, computerl(L));
}



bool SteadyState::inSupport(EL el)
{
    return inSupport(el, Emin(el.L));
}

bool SteadyState::inSupport(EL el, double EminL)
{
    double E = el.E;
    double L = el.L;

    if((E-e0) > numeric_limits<double>::epsilon() || (L-L0) < -numeric_limits<double>::epsilon()  || (L-Lmax) > numeric_limits<double>::epsilon())
        return false;

    if((E-EminL) < -numeric_limits<double>::epsilon())
        return false;

    return true;
}

bool SteadyState::atSupportboundary(EL el)
{
    double E = el.E;
    double L = el.L;

    //Check bottom part L=L0
    if(fabs(L-L0) < numeric_limits<double>::epsilon() && E-pot_atr(0.) > -numeric_limits<double>::epsilon() && E-e0 < numeric_limits<double>::epsilon())
        return true;
    //Check right part E=E0
    if(fabs(E-e0) < numeric_limits<double>::epsilon() && L-L0 > -numeric_limits<double>::epsilon() && L-Lmax < numeric_limits<double>::epsilon())
        return true;
    //Check near circular curve
    if(E-pot_atr(0.) > -numeric_limits<double>::epsilon() && E-e0 < numeric_limits<double>::epsilon()
       && L-L0 > -numeric_limits<double>::epsilon() && L-Lmax < numeric_limits<double>::epsilon())
        if(fabs(E-Emin(L)) < numeric_limits<double>::epsilon())
            return true;

    return false;
}

void SteadyState::plotrad(string fname, int radsteps)
{
    plotrad(fname, radsteps, rmin, rmax);
}


void SteadyState::plotrad(string fname, int radsteps, double rplotmin, double rplotmax)
{
    stringstream outer;

    //First some global stuff
    outer << setprecision(7) << "# Rmin = " << getrmin() << "\n";
    outer << "# Rmax = " << getrmax() << "\n";
    outer << "# Total Mass = " << getmtot() << "\n";
    outer << "# U_0(0) = " << pot_atr(0.) << "\n";
    outer << "# Cut Off Energy = " << gete0() << "\n";
    outer << "# Kappa = " << getkappa() << "\n";

    if(radsteps > 1)
    {
        outer << "#\n# Now radial data in Format |1) radius |2) potential |3) rho |4) mass |5) rhofalse \n#\n";

        double drp = (rplotmax-rplotmin)/((double)(radsteps-1));

        for(int i = 0; i < radsteps; i++)
        {
            double rloc = rplotmin + ((double)i)*drp;
            outer << setprecision(12) << rloc << " " << pot_atr(rloc) << " " << rho_atr(rloc) << " "
                  << mass_atr(rloc) << " " << 4.*M_PI*rloc*rloc*rho_atr(rloc) << "\n";
        }
    }

    //Write data in file
    string prefixer = "ststr";
    clearOutputFile(fname, prefixer);
    appendSSToOutputFile(fname, prefixer, &outer);
}


void SteadyState::ploteffpot(string fname, int amountL, int radsteps, double rplotmin, double rplotmax)
{
    stringstream outer;

    //First some global stuff
    outer << setprecision(7) << "# Rmin = " << getrmin() << "\n";
    outer << "# Rmax = " << getrmax() << "\n";
    outer << "# Total Mass = " << getmtot() << "\n";
    outer << "# U_0(0) = " << pot_atr(0.) << "\n";
    outer << "# Cut Off Energy = " << gete0() << "\n";
    outer << "# Kappa = " << getkappa() << "\n";

    if(radsteps > 1 && amountL > 0)
    {
        double drp = (rplotmax-rplotmin)/((double)(radsteps-1));

        for(int j = 0; j < amountL; j++)
        {
            double Lloc = L0 + (((double)j)/((double)(amountL-1)))*(Lmax-L0);

            outer << setprecision(3) << "L=" << Lloc << "\n"; //Columnheader
            for(int i = 0; i < radsteps; i++)
            {
                double rloc = rplotmin + ((double)i)*drp;

                if(Lloc > numeric_limits<double>::epsilon() && rloc < numeric_limits<double>::epsilon())
                    outer << setprecision(12) << rloc << " " << numeric_limits<double>::max() << " "
                          << numeric_limits<double>::lowest() << " " << numeric_limits<double>::max() << "\n";
                else
                    outer << setprecision(12) << rloc << " " << effpot(Lloc, rloc) << " "
                          << effpotprime(Lloc, rloc) << " " << effpotprimeprime(Lloc, rloc) << "\n";
            }

            outer << "\n\n";
        }
    }

    //Write data in file
    string prefixer = "effpotr";
    clearOutputFile(fname, prefixer);
    appendSSToOutputFile(fname, prefixer, &outer);
}



struct StStGlobData { //Wrapper Struct for global data of stst
    double mtot;
    double rmin;
    double rmax;
    double e0;
    double L0;
    double Lmax;
    double kappa;

    StStGlobData() //default constructor
    {
        mtot = -1.; rmin = -1.; rmax = -1.;
        e0 = -1.; L0 = -1.; Lmax = -1.;
        kappa = -1.;
    }

    StStGlobData(double mtoti, double rmini, double rmaxi, double e0i, double L0i, double Lmaxi, double kappai)
    {
        mtot = mtoti; rmin = rmini; rmax = rmaxi;
        e0 = e0i; L0 = L0i; Lmax = Lmaxi;
        kappa = kappai;
    }

    string tostring()
    {
        stringstream returnerss;
        returnerss << mtot << " " << rmin << " " << rmax << " "
                   << e0 << " " << L0 << " " << Lmax << " " << kappa;
        return returnerss.str();
    }
};

void ststfamilyparahelper(Ansatz* phi, double kappa, double dr, StStGlobData *data)
{//Helper Function which "returns" just essential data of steady state
    SteadyState stst(kappa, dr, phi);
    *data = StStGlobData(stst.getmtot(), stst.getrmin(), stst.getrmax(),
                         stst.gete0(), stst.getL0(), stst.getLmax(),
                         kappa);
    return;
}

void ststfamilygeneral(string fname, vector<Ansatz*> phis, vector<double> kappas, vector<string> adddata, double dr)
{
    int nrststs = kappas.size();
    int batchsize = 20;

    // Calculate the number of batches
    int numBatches = (nrststs + batchsize - 1)/batchsize;

    stringstream outer;
    string prefixer = "ststfam";
    clearOutputFile(fname, prefixer);

    vector<thread> threads;

    for(int i = 0; i < numBatches; i++)
    {
        threads.clear();

        int nrststs_batch = min(batchsize, nrststs - i*batchsize);
        vector<StStGlobData> data_batch(nrststs_batch);

        for(int j = 0; j < nrststs_batch; j++)
        {
            int l = i*batchsize + j;
            threads.emplace_back(thread(ststfamilyparahelper, phis.at(l), kappas.at(l), dr, &data_batch.at(j)));
        }

        for (auto &thread : threads)
            thread.join();

        for(int j = 0; j < nrststs_batch; j++)
        {
            int l = i*batchsize + j;
            outer << adddata.at(l) << " " << data_batch.at(j).tostring() << "\n";
        }

        //Write data immediately into file and clear sstream
        appendSSToOutputFile(fname, prefixer, &outer);
        outer.str(string()); //Clear string
    }

    return;
}

void ststfamily_kappa(string fname, Ansatz* phi, vector<double> kappas, double dr)
{
    vector<string> adddata;
    vector<Ansatz*> phis;
    for(int i = 0; i < kappas.size(); i++)
    {
        stringstream helper;
        helper << kappas.at(i);
        adddata.push_back(helper.str());
        phis.push_back(phi);
    }

    ststfamilygeneral(fname, phis, kappas, adddata, dr);
}


void ststfamily_ks(string fname, vector<double> ks, double dr)
{
    vector<string> adddata;
    vector<Ansatz*> phis;
    vector<double> kappas;
    for(int i = 0; i < ks.size(); i++)
    {
        stringstream helper;
        helper << ks.at(i);
        adddata.push_back(helper.str());
        Ansatz* phi = new PolyAnsatz(ks.at(i), 0., 0.);
        phis.push_back(phi);
        kappas.push_back(1.);
    }

    ststfamilygeneral(fname, phis, kappas, adddata, dr);
}


void kappa4r1_parahelper(double k, double ell, double dr, double *result)
{
    Ansatz* phi = new PolyAnsatz(k, ell, 0.);
    SteadyState stst(1., dr, phi);

    *result = pow(stst.getrmax(), (4.*ell+4.)/(2.*k+2.*ell+1.));
}

double kappa4r1(double k, double ell, double dr)
{
    double returner;
    kappa4r1_parahelper(k, ell, dr, &returner);
    return returner;
}

vector<double> kappas4r1(vector<double> ks, double ell, double dr)
{
    vector<double> kappas(ks.size(), -1.);

    int nrststs = ks.size();
    int batchsize = 50;

    // Calculate the number of batches
    int numBatches = (nrststs + batchsize - 1)/batchsize;

    vector<thread> threads;
    for(int i = 0; i < numBatches; i++)
    {
        threads.clear();

        int nrststs_batch = min(batchsize, nrststs - i*batchsize);

        for(int j = 0; j < nrststs_batch; j++)
        {
            int l = i*batchsize + j;
            threads.emplace_back(thread(kappa4r1_parahelper, ks.at(l), ell, dr, &kappas.at(l)));
        }

        for (auto &thread : threads)
            thread.join();
    }
    threads.clear();

    return kappas;
}

SteadyState polykappascaled(double k, double ell, double dr)
{
    double kappa  = kappa4r1(k, ell, dr);
    Ansatz* phi;
    if(fabs(ell) < numeric_limits<double>::epsilon())
        phi = new IsoPolyAnsatz(k);
    else
        phi = new PolyAnsatz(k, ell, 0.);

    return SteadyState(kappa, dr, phi);
}



void ststfamily_ks_radnorm(string fname, vector<double> ks, double dr)
{
    vector<string> adddata;
    vector<Ansatz*> phis;
    vector<double> kappas = kappas4r1(ks, 0., dr);

    for(int i = 0; i < ks.size(); i++)
    {
        stringstream helper;
        helper << ks.at(i);
        adddata.push_back(helper.str());
        Ansatz* phi = new PolyAnsatz(ks.at(i), 0., 0.);
        phis.push_back(phi);
    }

    ststfamilygeneral(fname, phis, kappas, adddata, dr);
}





SteadyState SteadyState::rescaleConstRadius()
{
    phi->setConstant(rmax*rmax*phi->C);
    return SteadyState(e0-potati(0), pot.getdr(), phi);
}




double SteadyState::computerminus(EL el)
{
    double rloc = 0.;
    if(el.L < numeric_limits<double>::epsilon())
        return 0.;
    if(el.E > -numeric_limits<double>::epsilon())
        return numeric_limits<double>::max();

    if(el.E-Emin(el.L) < numeric_limits<double>::epsilon())
        return computerl(el.L);

    double dr = getdr();
    int i = 1;
    while(rloc < 1000000.)
    {
        rloc = ((double)i)*dr;
        if(effpot(el.L, rloc) < el.E)
        {
            if(i == 1)
                return dr;

            double psilleft = effpot(el.L,rloc-dr);
            double psilright = effpot(el.L,rloc);

            if(psilleft < el.E) //Shouldnt be the case
                return rloc-dr;

            double convcoeff = (el.E - psilleft)/(psilright-psilleft); //Interpolate effpot linearly

            return rloc - convcoeff*dr;
        }

        i++;
    }

    return numeric_limits<double>::max();
}


double SteadyState::computerplus(EL el)
{
    if(el.L < -numeric_limits<double>::epsilon())
        return 0.;
    if(el.E > -numeric_limits<double>::epsilon())
        return numeric_limits<double>::max();

    if(el.E-Emin(el.L) < numeric_limits<double>::epsilon())
        return computerl(el.L);

    double rl = computerl(el.L);
    double rloc = rl;
    double dr = getdr();
    int i = 0;
    while(rloc < 1000000.)
    {
        rloc = rl + ((double)i)*dr;
        if(effpot(el.L, rloc) > el.E)
        {
            if(i == 1)
                return dr;

            double psilleft = effpot(el.L,rloc-dr);
            double psilright = effpot(el.L,rloc);

            if(psilleft > el.E) //Shouldnt be the case
                return rloc-dr;

            double convcoeff = (el.E - psilleft)/(psilright-psilleft); //Interpolate effpot linearly

            return rloc - convcoeff*dr;
        }
        i++;
    }

    return numeric_limits<double>::max();
}


double SteadyState::radreflection(double L, double r)
{
    if(L < numeric_limits<double>::epsilon())
        return r;
    if(effpot(L,r) > -numeric_limits<double>::epsilon())
        return numeric_limits<double>::max();

    double rl = computerl(L);
    if(r-rl > -numeric_limits<double>::epsilon())
        return r;

    return computerplus(EL(effpot(L,r),L));
}

double SteadyState::GL(double L, double r)
{
    return GL(L, r, Emin(L));
}

double SteadyState::GL(double L, double r, double EminL)
{
    if(r < -numeric_limits<double>::epsilon())
        return 0.;
    if(r < numeric_limits<double>::epsilon() && L > numeric_limits<double>::epsilon())
        return 0.;

    double psil = effpot(L,r)-EminL;
    double psilp = effpotprime(L,r);
    double psilpp = effpotprimeprime(L,r);

    if(fabs(psilp) < numeric_limits<double>::epsilon())
        return 0.;

    return ( 1. - 2.*(psil*psilpp)/(psilp*psilp) );
}


double SteadyState::HL(double L, double r)
{
    return HL(L, r, Emin(L));
}

double SteadyState::HL(double L, double r, double EminL)
{
    if(r < -numeric_limits<double>::epsilon())
        return 0.;
    if(r < numeric_limits<double>::epsilon() && L < numeric_limits<double>::epsilon())
        return -2.;
    if(r < numeric_limits<double>::epsilon() && L > numeric_limits<double>::epsilon())
        return 0.;

    double psil = effpot(L,r)-EminL;
    double psilp = effpotprime(L,r);
    double psilpp = effpotprimeprime(L,r);

    if(fabs(psilp) < numeric_limits<double>::epsilon())
        return 0.;

    return ( 1. - 2.*(psil*(psilpp+(2./r)*psilp))/(psilp*psilp) );
}

double SteadyState::GLref(double L, double r)
{
    return GLref(L, r, computerl(L), Emin(L));
}

double SteadyState::GLref(double L, double r, double rL, double EminL)
{
    if(r < -numeric_limits<double>::epsilon())
        return 0.;
    if(r < numeric_limits<double>::epsilon() && L > numeric_limits<double>::epsilon())
        return 0.;
    if(L < numeric_limits<double>::epsilon())
        return GL(L, r, EminL);
    if(rL-r < numeric_limits<double>::epsilon())
        return 0.;

    double zetaL = radreflection(L,r);

    return (GL(L,r,EminL) - GL(L,zetaL,EminL)*effpotprime(L,r)/effpotprime(L,zetaL));
}


double SteadyState::HLref(double L, double r)
{
    return HLref(L, r, computerl(L), Emin(L));
}

double SteadyState::HLref(double L, double r, double rL, double EminL)
{
    if(r < -numeric_limits<double>::epsilon())
        return 0.;
    if(r < numeric_limits<double>::epsilon() && L > numeric_limits<double>::epsilon())
        return 0.;
    if(r < numeric_limits<double>::epsilon() && L < numeric_limits<double>::epsilon())
        return numeric_limits<double>::lowest();
    if(L < numeric_limits<double>::epsilon())
        return (HL(L, r, EminL)/(r*r));
    if(rL-r < numeric_limits<double>::epsilon())
        return 0.;

    double zetaL = radreflection(L,r);

    return (HL(L,r,EminL)/(r*r) - (HL(L,zetaL,EminL)/(zetaL*zetaL))*effpotprime(L,r)/effpotprime(L,zetaL));
}


void GLrefHLrefPosParahelper(double L, double dr, double drEPS, double valEPS, SteadyState* stst, int* resulter)
{ //Helper function for function below; Check Positivity at fixed L
    double rloc = 0.;
    double rminus = stst->computerminus(EL(stst->gete0(),L));
    double rl = stst->computerl(L);
    double eminloc = stst->Emin(L);

    rloc = rl-drEPS; //Maximal r
    while(rloc > rminus)
    {
        double GLrefloc = stst->GLref(L, rloc, rl, eminloc);
        double HLrefloc = stst->HLref(L, rloc, rl, eminloc);

        if(GLrefloc+valEPS < 0. || HLrefloc+valEPS < 0.)
        {
            *resulter = 0;
            cout << "Violated GLref&HLref Positivity at L=" << L << ", r=" << rloc << ": "
                 << "GLref=" << GLrefloc << ", HLref=" << HLrefloc << "."
                 << "Note r_L=" << rl << "." << endl;
            return;
        }
        rloc -= dr;
    }
    *resulter = 1;
}

bool SteadyState::GLrefHLrefPos(int amountL, double dr, double drEPS, double valEPS)
{
    vector<thread> threads;

    vector<int> resulters(amountL, 1);

    for(int i = 0; i < amountL; i++)
    {
        double Lloc = L0 + (((double)i)/((double)(amountL)))*(Lmax-L0);

        threads.emplace_back(thread(GLrefHLrefPosParahelper, Lloc, dr, drEPS, valEPS, this, &resulters.at(i)));
    }

    for(auto &thread : threads)
        thread.join();

    for(int i = 0; i < resulters.size(); i++)
        if(resulters.at(i) == 0)
            return false;

    return true;
}

void SteadyState::plotGLHL(string fname, int amountL, int radsteps)
{
    stringstream ss;

    for(int i = 0; i < amountL; i++)
    {
        double Lloc = L0 + (((double)i)/((double)(amountL)))*(Lmax-L0);

        double rminus = computerminus(EL(e0,Lloc));
        double rplus = computerplus(EL(e0,Lloc));
        double eminl = Emin(Lloc);
        double rl = computerl(Lloc);

        ss << "# L = " << Lloc << "\n";
        ss << "# r_-(E_0,L) = " << rminus << "\n";
        ss << "# r_+(E_0,L) = " << rplus << "\n";
        ss << "# E_min(L) = " << eminl << "\n";
        ss << "# r_L = " << rl << "\n";

        ss << setprecision(3) << "L=" << Lloc << "\n"; //Columnheader

        for(int j = 0; j < radsteps; j++)
        {
            double rloc = rminus + (((double)j)/((double)(radsteps-1)))*(rplus-rminus);
            ss << setprecision(12) << rloc << " " << GL(Lloc, rloc, eminl) << " " << HL(Lloc, rloc, eminl) << " "
               << GLref(Lloc, rloc, rl, eminl) << " " << HLref(Lloc, rloc, rl, eminl) << "\n";
        }

        ss << "\n\n";
    }

    string fnameprefix = "GLHL";

    clearOutputFile(fname, fnameprefix);
    appendSSToOutputFile(fname, fnameprefix, &ss);

    return;
}


double SteadyState::dETnearcirc(double L)
{
    if(L < numeric_limits<double>::epsilon())
        return 0.; //Not save for L=0

    double rl = computerl(L);
    double psilpp = effpotprimeprime(L, rl);
    double psilppp = effpotprime3(L, rl);
    double psilpppp = effpotprime4(L, rl);

    return (5.*M_PI/12.)*psilppp*psilppp/pow(psilpp, 3.5) - .25*M_PI*psilpppp/pow(psilpp, 2.5);
}

double SteadyState::dLTnearcirc(double L)
{
    if(L < numeric_limits<double>::epsilon())
        return 0.; //Not save for L=0

    double rl = computerl(L);
    double psilpp = effpotprimeprime(L, rl);
    double psilppp = effpotprime3(L, rl);
    double psilpppp = effpotprime4(L, rl);

    return (-M_PI/(rl*rl*rl*rl*pow(psilpp, 3.5))) * (3.*psilpp*psilpp + rl*psilpp*psilppp
                                                     - rl*rl*psilpp*psilpppp/8. + rl*rl*psilppp*psilppp*5./24.);
}

double SteadyState::kunzerhs(double Tmax)
{
    return 4.*M_PI*M_PI/(Tmax*Tmax);
}

XV SteadyState::charSysRHS(XV xv)
{
    double rloc = xv.comp_r();
    double Uprimeloc = potprime_atr(rloc);

    XV returner;
    if(rloc < numeric_limits<double>::epsilon())
        returner = XV(xv.v1, xv.v2, xv.v3, 0., 0., 0.);
    else
        returner = XV(xv.v1, xv.v2, xv.v3, -(xv.x1/rloc)*Uprimeloc, -(xv.x2/rloc)*Uprimeloc, -(xv.x3/rloc)*Uprimeloc);
    return returner;
}




pair<double,double> SteadyState::trajectoriesolver(double rinit, double winit, double L, double dt,
                                              bool stopatradius, double stoprad, bool stopaftertime, double stoptime, bool stopafterperiod,
                                              bool compInt, double *integral, RadFct *Uprime)
{
    if(stopafterperiod && L > numeric_limits<double>::epsilon() && fabs(winit) > numeric_limits<double>::epsilon())
    {
        cout << endl << endl << "Error when solving char System of stst to compute Period..." << endl << endl;
        exit(10);
    }
    if(stopaftertime && stoptime < numeric_limits<double>::epsilon())
    {
        cout << endl << endl << "Error when evolving char System over time; time t=" << stoptime << " too small" << endl << endl;
        exit(11);
    }
    if(stopatradius && rinit > stoprad)
    {
        pair<double, double> returner {0.,0.};
        return returner;
    }

    XV xv(rinit, winit, L);

    double t = 0.;
    int tcounter = 0;

    while(tcounter < numeric_limits<int>::max()-2)
    {
        tcounter++;
        t = ((double)tcounter)*dt;

        //Compute RK4 Terms
        XV k1 = charSysRHS(xv);

        XV k2 = charSysRHS(addWithFactor(xv, .5*dt, k1));

        XV k3 = charSysRHS(addWithFactor(xv, .5*dt, k2));

        XV k4 = charSysRHS(addWithFactor(xv, dt, k3));

        if(compInt)
        { //Calculate Integral over w*Uprime; integrand is only evaluated at full dt steps...
            double rloc = xv.comp_r();

            if(tcounter == 1)
                *integral += .5*xv.comp_w(rloc)*Uprime->atr(rloc)*dt;
            else
                *integral += xv.comp_w(rloc)*Uprime->atr(rloc)*dt;
        }

        xv.rk4_update(k1, k2, k3, k4, dt);

        double rloc = xv.comp_r();

        if(stopatradius && rloc > stoprad)
        {
            pair<double, double> returner {t,t};
            return returner;
        }

        if(stopaftertime && (stoptime-t) < numeric_limits<double>::epsilon())
        {
            double wloc = xv.comp_w(rloc);

            if(compInt)
                *integral += .5*wloc*Uprime->atr(rloc)*dt;

            pair<double, double> returner {rloc, wloc};
            return returner;
        }

        if(stopafterperiod && xv.comp_w() > 0.)
        {
            pair<double, double> returner {2.*t, 2.*t};
            return returner;
        }
    }

    exit(570403);
}

pair<double,double> SteadyState::trajectoriesolver(double rinit, double winit, double L, double dt, double T)
{
    double inthelper = 0.;
    RadFct radfcthelper;

    return trajectoriesolver(rinit, winit, L, dt, false, 0.,
                             true, T, false, false, &inthelper, &radfcthelper);
}

pair<double,double> SteadyState::trajectoriesolver(double rinit, double winit, double L, double dt, double T,
                                                   double *integral, RadFct *Uprime)
{
    return trajectoriesolver(rinit, winit, L, dt, false, 0.,
                             true, T, false, true, integral, Uprime);
}


double SteadyState::computeT(EL el, double dt)
{
    if((el.E - Emin(el.L)) < numeric_limits<double>::epsilon())
    {
        double psilpp = effpotprimeprime(el.L, computerl(el.L));
        if(psilpp < numeric_limits<double>::epsilon())
            return numeric_limits<double>::max();
        else
        {
            if(el.L < numeric_limits<double>::epsilon())
                return M_PI/sqrt(psilpp);
            else
                return 2.*M_PI/sqrt(psilpp);
        }
    }


    double inthelper = 0.;
    RadFct radfcthelper;

    return trajectoriesolver(computerplus(el), 0., el.L, dt,
                             false, 0., false, 0., true,
                             false, &inthelper, &radfcthelper).first;
}


double SteadyState::computeTCorner(double dt)
{
    return computeT(EL(e0, L0), dt);
}


void computeT_parahelper(EL el, double dt, SteadyState* stst, double *res)
{
    *res = stst->computeT(el, dt);
    return;
}

vector<double> SteadyState::computeT(vector<EL> els, double dt)
{
    int nrststs = els.size();
    int batchsize = 500;

    // Calculate the number of batches
    int numBatches = (nrststs + batchsize - 1)/batchsize;

    vector<double> returner(els.size(), 0.);

    vector<thread> threads;

    for(int i = 0; i < numBatches; i++)
    {
        threads.clear();

        int nrststs_batch = min(batchsize, nrststs - i*batchsize);

        for(int j = 0; j < nrststs_batch; j++)
        {
            int l = i*batchsize + j;
            threads.emplace_back(thread(computeT_parahelper, els.at(l), dt, this, &returner.at(l)));
        }

        for (auto &thread : threads)
            thread.join();
    }

    return returner;
}



vector<EL> SteadyState::createELtriangle(int Esteps, int Lsteps)
{
    double Emintot = Emin(L0);
    vector<EL> returner;
    for(int i = 0; i < Lsteps; i++)
    {
        double Lloc = L0 + (((double)i)/((double)Lsteps))*(Lmax-L0);
        double eminlloc = Emin(Lloc);
        if(i > 0) //Make sure that minimal energy values are also in triangle
            returner.push_back(EL(eminlloc, Lloc));

        for(int j = 0; j <= Esteps; j++)
        {
            double Eloc = Emintot + (((double)j)/((double)Esteps))*(e0-Emintot);
            EL elloc(Eloc, Lloc);
            if(inSupport(elloc, eminlloc))
                returner.push_back(elloc);
        }
    }
    returner.push_back(EL(e0, Lmax));
    return returner;
}

vector<EL> SteadyState::createELtriangle(int Esteps, int Lsteps, int minELs)
{
    vector<EL> returner = createELtriangle(Esteps, Lsteps);
    //cout << "Triangle with min " << minELs << ", size is " << returner.size() << endl;

    if(minELs <= returner.size())
        return returner;
    else
    {
        int Estepsold = Esteps;
        int Estepsnew = Estepsold + 10;
        int Lstepsnew = Lsteps + (Lsteps*Estepsnew-Lsteps*Estepsold)/Estepsold;
        return createELtriangle(Estepsnew, Lstepsnew, minELs);
    }
}


pair<double,double> SteadyState::computeTminTmax(double dt, int Esteps, int Lsteps)
{
    return computeTminTmax(dt, computeT( createELtriangle(Esteps, Lsteps), dt) );
}

pair<double,double> SteadyState::computeTminTmax(double dt, int Esteps, int Lsteps, int minELs)
{
    return computeTminTmax(dt, computeT( createELtriangle(Esteps, Lsteps, minELs), dt) );
}

pair<double,double> SteadyState::computeTminTmax(double dt, vector<double> Tvals)
{
    double Tmin = numeric_limits<double>::max();
    double Tmax = numeric_limits<double>::lowest();

    for(int i = 0; i < Tvals.size(); i++)
    {
        if(Tvals.at(i) < Tmin)
            Tmin = Tvals.at(i);
        if(Tvals.at(i) > Tmax)
            Tmax = Tvals.at(i);
    }

//    cout << endl << endl << "Away from Minimal Energies: Tmin=" << Tmin << ", Tmax=" << Tmax
//         << " (" << Tvals.size() << " values)."  << endl;

    //Search on Minimal energy curve and at (E_0,L_0) -> The former is (computationally) fast due to explicit formula
    vector<double> Tadditional;
    Tadditional.push_back(computeTCorner(dt));

    int nrLs = 100;
    for(int i = 0; i <= nrLs; i++)
    {
        double Lloc = L0 + (((double)i)/((double)nrLs))*(Lmax-L0);
        Tadditional.push_back(computeT(EL(Emin(Lloc),Lloc), dt));
    }

    double Tminadd = numeric_limits<double>::max();
    double Tmaxadd = numeric_limits<double>::lowest();
    for(int i = 0; i < Tadditional.size(); i++)
    {
        if(Tadditional.at(i) < Tminadd)
            Tminadd = Tadditional.at(i);
        if(Tadditional.at(i) > Tmaxadd)
            Tmaxadd = Tadditional.at(i);
    }

//    cout << "At E=Emin and Corner: Tmin=" << Tminadd << ", Tmax=" << Tmaxadd << "." << endl;

    Tmin = min(Tmin, Tminadd);
    Tmax = max(Tmax, Tmaxadd);

    pair<double,double> returner {Tmin, Tmax};
    return returner;
}


bool SteadyState::TmaxatCorner(double dt, double Tmax, double EPS)
{
    return ((computeTCorner(dt)+EPS-Tmax) > 0.);
}

bool SteadyState::Tmaxatboundary(vector<EL> els, vector<double> Ts, double EPS)
{
    double Tmaxloc = numeric_limits<double>::lowest();
    EL tmaxer(e0,L0);

    for(int i = 0; i < els.size(); i++)
    {
        if(Ts.at(i) > Tmaxloc)
        {
            Tmaxloc = Ts.at(i);
            tmaxer = els.at(i);
        }
    }

    if(atSupportboundary(tmaxer))
        return true;

    for(int i = 0; i < els.size(); i++)
    {
        if(Ts.at(i) > Tmaxloc-EPS)
            if(atSupportboundary(els.at(i)))
                return true;
    }

    return false;
}

bool SteadyState::TincE(vector<double> Ts, vector<EL> ELs, double EPS)
{
    //First determine different L-values and associated indices
    vector<double> Ls;
    vector<vector<int>> inds;

    for(int i = 0; i < ELs.size(); i++)
    {
        bool newL = true;
        int Lind = -1;
        EL ELloc = ELs.at(i);
        double Lloc = ELloc.L;
        for(int j = 0; j < Ls.size(); j++)
        {
            if(abs(Ls.at(j)-Lloc) < numeric_limits<double>::epsilon())
            {
                newL = false;
                Lind = j;
            }
        }

        if(newL)
        {
            Ls.push_back(Lloc);
            vector<int> newLind{i};
            inds.push_back(newLind);
        }
        else
            inds.at(Lind).push_back(i);
    }

    /*
    cout << "We got " << Ls.size() << " different Ls" << endl;
    for(int j = 0; j < Ls.size(); j++)
        cout << j <<  "th L is " << Ls.at(j) << ", amount of Es is " << inds.at(j).size() << "." << endl;
    */

    //For every fixed L-value go through E's and check T-monotonicity
    for(int j = 0; j < Ls.size(); j++)
    {
        for(int i = 0; i < inds.at(j).size(); i++)
        {
            double Eloc = ELs.at(inds.at(j).at(i)).E;
            double Tloc = Ts.at(inds.at(j).at(i));

            for(int i2 = 0; i2 < inds.at(j).size(); i2++)
            {
                double Eloc2 = ELs.at(inds.at(j).at(i2)).E;
                double Tloc2 = Ts.at(inds.at(j).at(i2));

                if(i2 != i)
                {
                    if(Eloc2 < Eloc && Tloc2-Tloc-EPS > 0.)
                    {
                        cout << "Violated Monotonicity: L=" << Ls.at(j) << ". T(" << Eloc << ",L)=" << Tloc
                             << ", T(" << Eloc2 << ",L)=" << Tloc2 << ". Note EminL=" << Emin(Ls.at(j)) << "." << endl;
                        return false;
                    }
                    if(Eloc2 > Eloc && Tloc2-Tloc+EPS < 0.)
                    {
                        cout << "Violated Monotonicity: L=" << Ls.at(j) << ". T(" << Eloc << ",L)=" << Tloc
                             << ", T(" << Eloc2 << ",L)=" << Tloc2 << ". Note EminL=" << Emin(Ls.at(j)) << "." << endl;
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

bool SteadyState::TdecL(vector<double> Ts, vector<EL> ELs, double EPS)
{ //Be clever and just call above function after switching roles of E and L and multiplying T with (-1)
    vector<EL> ELstrick;
    for(int i = 0; i < ELs.size(); i++)
        ELstrick.push_back(EL(ELs.at(i).L, ELs.at(i).E));
    vector<double> Tstrick;
    for(int i = 0; i < Ts.size(); i++)
        Tstrick.push_back((-1.)*Ts.at(i));

    return TincE(Tstrick, ELstrick, EPS);
}


vector<bool> SteadyState::Tmaxmons(double dt, int Esteps, int Lsteps, double EPS)
{
    return Tmaxmons(dt, createELtriangle(Esteps, Lsteps), EPS);
}

vector<bool> SteadyState::Tmaxmons(double dt, int Esteps, int Lsteps, int minELs, double EPS)
{
    return Tmaxmons(dt, createELtriangle(Esteps, Lsteps, minELs), EPS);
}

vector<bool> SteadyState::Tmaxmons(double dt, vector<EL> ELs, double EPS)
{
    return Tmaxmons(dt, ELs, computeT(ELs, dt), EPS);
}

vector<bool> SteadyState::Tmaxmons(double dt, vector<EL> ELs, vector<double> Tvals, double EPS)
{
    vector<bool> returner;
    returner.push_back(TmaxatCorner(dt, computeTminTmax(dt, Tvals).second, EPS));
    returner.push_back(Tmaxatboundary(ELs, Tvals, EPS));
    returner.push_back(TincE(Tvals, ELs, EPS));
    returner.push_back(TdecL(Tvals, ELs, EPS));

    return returner;
}

void SteadyState::plotTs(string fname, int Esteps, int Lsteps, double dt)
{
    plotTs(fname, createELtriangle(Esteps, Lsteps), dt);
}

void SteadyState::plotTs(string fname, int Esteps, int Lsteps, int minELs, double dt)
{
    plotTs(fname, createELtriangle(Esteps, Lsteps, minELs), dt);
}

void SteadyState::plotTs(string fname, vector<EL> els, double dt)
{
    plotTs(fname, els, computeT(els, dt), dt);
}



void SteadyState::plotTs(string fname, vector<EL> els, vector<double> Ts, double dt)
{
    pair<double, double> Tminmax = computeTminTmax(dt, Ts);

    stringstream ss;
    ss << "# Tmin = " << Tminmax.first << "\n";
    ss << "# Tmax = " << Tminmax.second << "\n";
    ss << "# T-Values used: " << Ts.size() << "\n";

    ss << "#\n"; //Now general global data
    ss << "# Rmin = " << getrmin() << "\n";
    ss << "# Rmax = " << getrmax() << "\n";
    ss << "# Total Mass = " << getmtot() << "\n";
    ss << "# U_0(0) = " << pot_atr(0.) << "\n";
    ss << "# Cut Off Energy = " << gete0() << "\n";
    ss << "#\n";

    for(int i = 0; i < Ts.size(); i++)
        ss << els.at(i).E << " " << els.at(i).L << " " << Ts.at(i) << "\n";

    string prefixer = "Ts";

    clearOutputFile(fname, prefixer);
    appendSSToOutputFile(fname, prefixer, &ss);
}


void SteadyState::plotTs_fixedL(string fname, int amountL, int Esteps, double dt)
{
    stringstream ss;
    string prefixer = "TfixedL";

    clearOutputFile(fname, prefixer);

    for(int i = 0; i < amountL; i++)
    {
        double Lloc = L0 + (((double)i)/((double)amountL))*(Lmax-L0);
        vector<EL> elsloc;
        double eminloc = Emin(Lloc);
        for(int j = 0; j < Esteps; j++)
        {
            double Eloc = eminloc + (((double)j)/((double)(Esteps-1)))*(e0-eminloc);

            elsloc.push_back(EL(Eloc, Lloc));
        }

        vector<double> Ts = computeT(elsloc, dt);

        ss << setprecision(3) << "L=" << Lloc << "\n";

        for(int j = 0; j < elsloc.size(); j++)
        {
            ss << setprecision(12) << elsloc.at(j).E << " " << Ts.at(j) << "\n";
        }

        ss << "\n\n";

        //Write data immediately into file and clear sstream
        appendSSToOutputFile(fname, prefixer, &ss);
        ss.str(string()); //Clear string
    }

    return;
}


void SteadyState::plotTs_fixedE(string fname, int amountE, int Lsteps, double dt)
{
    stringstream ss;
    string prefixer = "TfixedE";

    clearOutputFile(fname, prefixer);

    double eminall = Emin(L0);
    for(int i = 1; i <= amountE; i++)
    {
        double Eloc = eminall + (((double)i)/((double)amountE))*(e0-eminall);

        vector<EL> elsloc;
        double lmaxloc = getLmax(Eloc);
        for(int j = 0; j < Lsteps; j++)
        {
            double Lloc = L0 + (((double)j)/((double)(Lsteps-1)))*(lmaxloc-L0);

            elsloc.push_back(EL(Eloc, Lloc));
        }

        vector<double> Ts = computeT(elsloc, dt);

        ss << setprecision(3) << "E=" << Eloc << "\n";

        for(int j = 0; j < elsloc.size(); j++)
        {
            ss << setprecision(12) << elsloc.at(j).L << " " << Ts.at(j) << "\n";
        }

        ss << "\n\n";

        //Write data immediately into file and clear sstream
        appendSSToOutputFile(fname, prefixer, &ss);
        ss.str(string()); //Clear string
    }

    return;
}



//Helper function for general Tmin/Tmax Computation of family of steady states
void ststfamily_Ts_general(string fname, vector<Ansatz*> phis, vector<double> kappas, double dr, vector<string> adddata,
                           int Esteps, int Lsteps, int minELs, double dt, bool seperatePlotFiles,
                           double EPST,
                           bool withGHref, int amountL, double drrefs, double drrefsEPS, double valrefsEPS)
{
    stringstream ss;
    string prefixer = "Tfam";

    clearOutputFile(fname, prefixer);

    for(int i = 0; i < kappas.size(); i++)
    {
        cout << adddata.at(i) << ": ";

        SteadyState stst(kappas.at(i), dr, phis.at(i));
        cout << endl << "DONE with stst computation " << endl;

        vector<EL> DEL = stst.createELtriangle(Esteps, Lsteps, minELs);
        cout << endl << "Ended comp of triangle; " << DEL.size() << " values." << endl << endl;
        vector<double> Ts = stst.computeT(DEL, dt);

        cout << endl << "Ended comp of Ts" << endl << endl;

        pair<double,double> Tminmax = stst.computeTminTmax(dt, Ts);

        vector<bool> Tmaxmons = stst.Tmaxmons(dt, DEL, Ts, EPST);

        if(seperatePlotFiles)
        {
            string fnamesep = fname;
            fnamesep += adddata.at(i);
            stst.plotTs(fnamesep, DEL, Ts, dt);
            stst.plotTs_fixedL(fnamesep, 6, Esteps, dt);
            stst.plotTs_fixedE(fnamesep, 7, Lsteps, dt);
        }

        ss << adddata.at(i) << " " << Tminmax.first << " " << Tminmax.second << " " << Tminmax.second/Tminmax.first << " "
           << Tmaxmons.at(0) << " " << Tmaxmons.at(1) << " " << Tmaxmons.at(2) << " " << Tmaxmons.at(3);

        if(withGHref)
        {
            bool GLrefHLrefpos = stst.GLrefHLrefPos(amountL, drrefs, drrefsEPS, valrefsEPS);

            ss << " " << GLrefHLrefpos;
        }
        else
            ss << " " << -1;

        //Kunze Crit
        int NrKunze = (stst.getrmax()-stst.getrmin())/drrefs;
        double kunzelhs = stst.potprimeoverrmax(NrKunze);
        double kunzerhs = stst.kunzerhs(Tminmax.second);
        ss << " " << kunzelhs << " " << kunzerhs << " " << (kunzelhs<kunzerhs);

        ss << "\n";

        //Write data immediately into file and clear sstream
        appendSSToOutputFile(fname, prefixer, &ss);
        ss.str(string()); //Clear string

        cout << " DONE." << endl;
    }

    return;
}


void ststfamily_Ts_kappa(string fname, Ansatz* phi, vector<double> kappas, double dr,
                         int Esteps, int Lsteps, int minELs, double dt, bool seperatePlotFiles,
                         double EPST,
                         bool withGHref, int amountL, double drrefs, double drrefsEPS, double valrefsEPS)
{
    vector<string> adddata;
    vector<Ansatz*> phis;
    for(int i = 0; i < kappas.size(); i++)
    {
        stringstream helper;
        helper << kappas.at(i);
        adddata.push_back(helper.str());
        phis.push_back(phi);
    }

    ststfamily_Ts_general(fname, phis,  kappas, dr, adddata, Esteps, Lsteps, minELs, dt, seperatePlotFiles, EPST,
                          withGHref, amountL, drrefs, drrefsEPS, valrefsEPS);
}

void ststfamily_Ts_ks(string fname, vector<double> ks, double dr,
                      int Esteps, int Lsteps, int minELs, double dt, bool seperatePlotFiles,
                      double EPST,
                      bool withGHref, int amountL, double drrefs, double drrefsEPS, double valrefsEPS)
{
    vector<string> adddata;
    vector<Ansatz*> phis;
    vector<double> kappas(ks.size(), 1.);

    for(int i = 0; i < ks.size(); i++)
    {
        stringstream helper;
        helper << ks.at(i);
        adddata.push_back(helper.str());
        Ansatz* phi = new PolyAnsatz(ks.at(i), 0., 0.); //Use PolyAnsatz instead of IsoPolyAnsatz for faster computation
        phis.push_back(phi);
    }

    ststfamily_Ts_general(fname, phis,  kappas, dr, adddata, Esteps, Lsteps, minELs, dt, seperatePlotFiles, EPST,
                          withGHref, amountL, drrefs, drrefsEPS, valrefsEPS);
}

void ststfamily_Ts_ks_radnorm(string fname, vector<double> ks, double dr,
                      int Esteps, int Lsteps, int minELs, double dt, bool seperatePlotFiles,
                      double EPST,
                      bool withGHref, int amountL, double drrefs, double drrefsEPS, double valrefsEPS)
{
    vector<string> adddata;
    vector<Ansatz*> phis;
    vector<double> kappas = kappas4r1(ks, 0., dr);

    for(int i = 0; i < ks.size(); i++)
    {
        stringstream helper;
        helper << ks.at(i);
        adddata.push_back(helper.str());
        Ansatz* phi = new PolyAnsatz(ks.at(i),  0., 0.); //Use PolyAnsatz instead of IsoPolyAnsatz for faster computation
        phis.push_back(phi);
    }

    ststfamily_Ts_general(fname, phis,  kappas, dr, adddata, Esteps, Lsteps, minELs, dt, seperatePlotFiles, EPST,
                          withGHref, amountL, drrefs, drrefsEPS, valrefsEPS);
}






void clearOutputFile(string fname, string sprefix)
{
    string outputhelper = "";
    string fnametxt = sprefix + "_" + fname + ".txt";
    ofstream dataout(fnametxt.c_str());
    dataout << outputhelper;
    dataout.close();
}

void appendSSToOutputFile(string fname, string sprefix, stringstream *ss)
{
    string fnametxt = sprefix + "_" + fname + ".txt";
    ofstream dataout;
    dataout.open(fnametxt, ios_base::app);
    string output = (*ss).str();
    dataout << output;
    dataout.close();

    ss->str(string()); //Clear stringstream
}




