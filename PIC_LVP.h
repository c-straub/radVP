#ifndef PIC_LVP_H_INCLUDED
#define PIC_LVP_H_INCLUDED

/////////////////////////////////////////////////
/// File for LVP PIC
/////////////////////////////////////////////////

#include "PIC_general.h"
#include "steadystate.h"

using namespace std;

class LinPert { //Parent Class for linear perturbation, i.e., the initial data for the LVP
    public:
        virtual double atrwL(double r, double w, double L, SteadyState *stst) = 0; //returns value of perturbation; input r,w,L variables and Steady State (latter for additional weight)
}; //cf. below for different types of perturbations

class FieldsLVP : public Fields { //Class for Fields of LVP
    public:
        using Fields::Fields; //Inherit constructor of parent class; possible thanks to C++11

        double getEkinfree_seq(const vector<ParPos> *parpos, SteadyState *stst); //Computes kinetic part of free energy sequentially
        double getEpotlin_seq(const vector<ParPos> *parpos, SteadyState *stst); //Computes potential part of linearised energy sequentially
};

class LocsLVP : public Locs { //Class for local quantities used for parallelisation as well as RK4 Data
    public:
        vector<double> ekinfreeloc; //local parts of ekinfree
        vector<double> epotlinloc; //local parts of epotlin

        //RK4 Data for computing integral:
        vector<double> int_k1;
        vector<double> int_k2;
        vector<double> int_k3;

        LocsLVP(int am_parts, RadFct *rho);

        void compEkinfree_parahelper(int tid, const vector<ParPos> *parpos, SteadyState *stst); //computes part of ekinfree
        double get_Ekinfree_all(); //Returns sum of all ekinfreelocs

        void compEpotlin_parahelper(int tid, const vector<ParPos> *parpos, SteadyState *stst); //computes part of epotlin
        double get_Epotlin_all(); //Returns sum of all ekinfreelocs
};

class ParaDataLVP : public ParaData { //Class for data used in parallelised part of LVP PIC
    public:
        double efree; //(Total) free energy
        double efreei; //Initial (Total) free energy
        double ekinfree; //Kinetic part of free energy
        double epotfree; //Potential part of free energy

        double elin; //(Total) linearised energy
        double elini; //Initial linearised energy
        double ekinlin; //Kinetic part of linearised energy
        double epotlin; //Potential part of linearised energy

        double mstst; //Total Mass of underlying Steady State

        ParaDataLVP(double mtotin, double ekinfreein, double epotfreein, double ekinlinin, double epotlinin,
                    double mststin, double pot0in,
                    double dt, int tsteps); //Constructor, based on initial values of all quantities

        double get_efreeerrorrel(); //Returns relative efree error
        double get_elinerrorrel(); //Returns relative elin error

        void update_vars(int stepscur, double mtotcur, double pot0cur,
                         double ekinfreecur, double epotfreecur, double ekinlincur, double epotlincur,
                         double rmincur, double rmaxcur); //Function to update time dependent variables to current values
};


class OutersLVP : public Outers { //Class for Outputs of LVP
    public:
        OutersLVP(string fnamei);

        void data2console(ParaDataLVP *pada); //Write Current Data to Console

        void initwrite_globout(PICParams pp); //Writes initial information into globout stream
        void write_globout(ParaDataLVP *pada); //Writes global data into stream
};

string conTabularLVPHead(); //Returns head of console tabular
string conTabularLVPMidHead(); //Returns mid-head of console tabular
string conTabularLVPSeparator(); //Returns separator row of console tabular

/////////////////////////////////////////////////
/// Now come the key functions :
/////////////////////////////////////////////////

void PIC_LVP_fromFiles(); //Reads in Parameters for LVP PIC from input files and runs PIC

vector<ParPos> init_LVP(InitParams params, LinPert *pert, SteadyState *stst, Outers *outs); //Initialises particles

void PIC_LVP(InitParams ip, LinPert* pert, SteadyState* stst, PICParams pp, bool with_response); //Starts PIC for LVP, including initialisation of particles
void PIC_LVP(vector<ParPos> *parpos, SteadyState *stst, PICParams pp, OutersLVP *outs, bool with_response); //Same as above, but particles already initialised

void PIC_LVP_parahelper(const int tid, PICParams pp, bool withresponse,
                        SteadyState *stst,
                        vector<ParPos> *parpos, FieldsLVP *fields, LocsLVP *locz,
                        ParaDataLVP *pada,
                        OutersLVP *outs); //Helper Function for parallelised part of LVP PIC

void addSingleParticle_num(const ParPos *parpos, double *res, int fpow, bool overphiprime, bool factorU0, SteadyState *stst); //Adds influence of single particle to number res
    //with optional power of f, denominator |\varphi'|, and factor U_0(x)

///Different Types of Perturbations:

class LPwpowerphi: public LinPert { //Nr. 0, linear perturbation of form w^j*\varphi(E,L), where j is an integer
    public:
        int j; //exponent for w
        double atrwL(double r, double w, double L, SteadyState *stst);

        LPwpowerphi(int ji) {j=ji;}
};

class LPwpowerphiprime: public LinPert { //Nr. 1, linear perturbation of form w^j*\varphi'(E,L), where j is an integer
    public:
        int j; //exponent for w
        double atrwL(double r, double w, double L, SteadyState *stst);

        LPwpowerphiprime(int ji) {j=ji;}
};

class LPsqrtrwphiprime: public LinPert { //Nr. 2, linear perturbation of form sqrt(|r^j+w^j|)*\varphi'(E,L), where j is an integer
    public:
        int j; //exponent for r and w
        double atrwL(double r, double w, double L, SteadyState *stst);

        LPsqrtrwphiprime(int ji) {j=ji;}
};

class LPsqrtrwphiprimeawaywell: public LinPert { //Nr. 3, linear perturbation of form sqrt(|r^j+w^j|)*\varphi'(E)*chi(E,L), where j is an integer and chi(E,L) smoothly vanishes close to E=EminL
    public:
        int j; //exponent for w ; j = 0 means no w
        double convanish; //chi(E,L) vanishes for E <= EminL + convanish*(E_0-U_0(0)), chi(E,L) = 1 elsewhere
        double atrwL(double r, double w, double L, SteadyState *stst);

        LPsqrtrwphiprimeawaywell(int ji, double iconvanish) {j=ji; convanish=iconvanish;}
};

LinPert* pertfromfile(string parameterfilename); //read in perturbation from file
int pertfromfile_nr(string parameterfilename); //Returns perturbation number from file

void atrua_parahelper(double r, double u, double a, LinPert *pert, SteadyState *stst, double *resval); //Helper function to parallelise evaluation of perturbation at different (r,u,alpha)


#endif // PIC_LVP_H_INCLUDED
