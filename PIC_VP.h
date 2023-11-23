#ifndef PIC_VP_H_INCLUDED
#define PIC_VP_H_INCLUDED

#include "PIC_general.h"
#include "steadystate.h"


class FieldsVP : public Fields { //Class for Fields of VP
    public:
        using Fields::Fields; //Inherit constructor of parent class; possible thanks to C++11

        XV RHScharsys(double r, double w, double L); //Returns Right-Hand Side of Char Sys of these fields at given (r,w,L)
        XV RHScharsys(XV *xv);
        XV RHScharsys(ParPos *papo); //Returns Right-Hand Side of Char Sys of these fields, input is papo
};

class LocsVP : public Locs { //Class for local quantities used for parallelisation as well as RK4 Data
    public:
        //RK4 Data for computing new particle position:
        vector<XV> par_k1;
        vector<XV> par_k2;
        vector<XV> par_k3;

        LocsVP(int am_parts, RadFct *rho);
};

class ParaDataVP : public ParaData { //Class for data used in parallelised part of VP PIC
    public:
        double etot; //(Total) energy
        double etoti; //Initial (Total)  energy
        double ekin; //Kinetic part of energy
        double epot; //Potential part of energy
        double ekini; //Initial kinetic part of energy

        double rminplot; //Minimal Radius for Plotting macroscopic functions
        double rmaxplot; //Maximal Radius for Plotting macroscopic functions
        int nrdelparts; //Number of deleted particles

        ParaDataVP(double mtotin, double ekinin, double epotin, double pot0in,
                    double dt, int tsteps,
                    double rminplotin, double rmaxplotin); //Constructor, based on initial values of all quantities
        ParaDataVP() {;}; //Helper default constructur; just used to compute rmin and rmax later

        double get_etoterrorrel(); //Returns relative elin error

        void update_vars(int stepscur, double mtotcur, double pot0cur,
                         double ekincur, double epotcur,
                         double rmincur, double rmaxcur, double wmincur, double wmaxcur,
                         int nrdelpartscur); //Function to update time dependent variables to current values
};

class OutersVP : public Outers { //Class for Outputs of VP
    public:
        OutersVP(string fnamei);

        void data2console(ParaDataVP *pada); //Write Current Data to Console

        void initwrite_globout(PICParams pp); //Writes initial information into globout stream
        void write_globout(ParaDataVP *pada); //Writes global data into stream
};

string conTabularVPHead(); //Returns head of console tabular
string conTabularVPMidHead(); //Returns mid-head of console tabular
string conTabularVPSeparator(); //Returns separator row of console tabular

/////////////////////////////////////////////////
/// Now come the key functions :
/////////////////////////////////////////////////

void PIC_VP_fromFiles(); //Reads in Parameters for VP PIC from input files and runs PIC

vector<ParPos> init_VP(InitParams params, double pertamp, SteadyState *stst, Outers *outs); //Initialises particles

void PIC_VP(InitParams ip, double pertamp, SteadyState* stst, PICParams pp); //Starts PIC for VP, including initialisation of particles
void PIC_VP(vector<ParPos> *parpos, PICParams pp, OutersVP *outs); //Same as above, but particles already initialised

void PIC_VP_parahelper(const int tid, PICParams pp,
                       vector<ParPos> *parpos, FieldsVP *fields, LocsVP *locz,
                       ParaDataVP *pada,
                       OutersVP *outs); //Helper Function for parallelised part of VP PIC



#endif // PIC_VP_H_INCLUDED
