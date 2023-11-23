#define _USE_MATH_DEFINES //To use PI
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "steadystate.h"
#include "PIC_LVP.h"
#include "PIC_VP.h"

using namespace std; //Not good practice, but makes life easier

void datastst(); //This function gives an example for how to compute the isotropic polytropic steady state k=1=Rmax

void dataperiod(); //This function gives an example for how to compute & analyse the period function of the isotropic polytropic steady state k=1=Rmax

int main()
{
    ///Call Example of Steady State Function
    //datastst();

    ///Call Example of Period Function Function
    //dataperiod();

    ///Below are the function calls for the Particle-In-Cell Parts ; the parameters are read in from the various input .txt files

    int pictype;
    bool success = true;
    ifstream parameterin("in_basic.txt");
    if(!(parameterin >> pictype))
    {
        cout << endl << "Can't read in_basic.txt file" << endl;
        success = false;
    }
    parameterin.close();

    if(pictype == 0 && success) // pictype==0 means that we simulate non-linearised system
        PIC_VP_fromFiles();
    else if((pictype == 1 || pictype == 2) && success) // pictype==1 or 2 means that we simulate the linearised system or the pure transport equation
        PIC_LVP_fromFiles();
    else
        cout << endl << "Unexpected input from in_basic.txt file..." << endl;


    return 0;
}











void datastst()
{
    cout << endl << endl << "Create Isotropic Polytropic Steady State k=1=Rmax and write functions associated to it in files." << endl << endl;

    double drdef = 1e-6;
    int Nr = 1000;

    SteadyState ststk1 = polykappascaled(1., 0., drdef);
    ststk1.plotrad("k=1", Nr, 0., 2.);

    cout << endl << endl << "Finished Steady State k=1 Part." << endl << endl;

    return;
}


void dataperiod()
{
    cout << endl << endl << "Compute & analyse the Radial Period Function of the Isotropic Polytropic Steady State k=1=Rmax and write values of it in files." << endl << endl;

    double drdef = 1e-7;
    double dtdef = 1e-6;

    //use the following values for creating (E,L)-triangles
    int Esteps = 128;
    int Lsteps = 128;
    int minELs = 300000;

    SteadyState ststk1 = polykappascaled(1., 0., drdef);

    ststk1.ploteffpot("k=1", 7, 1000, 0., 2.);

    ststk1.plotGLHL("k=1", 6, 500);
    bool GLrefHLrefposk1 = ststk1.GLrefHLrefPos(50, 1000.*drdef, 1e-2, 1e-4);

    ststk1.plotTs_fixedL("k=1", 6, 300, dtdef);
    ststk1.plotTs_fixedE("k=1", 7, 300, dtdef);

    vector<EL> DELk1 = ststk1.createELtriangle(Esteps, Lsteps, minELs); //(E,L)-triangle for k=1
    vector<double> Tsk1 = ststk1.computeT(DELk1, dtdef);

    ststk1.plotTs("k=1", DELk1, Tsk1, dtdef);

    pair<double,double> Tminmaxk1 = ststk1.computeTminTmax(dtdef, Tsk1);
    cout << "For k=1:" << endl << "Tmin = " << Tminmaxk1.first << ", Tmax = " << Tminmaxk1.second << endl;

    vector<bool> Tmaxmonsk1 = ststk1.Tmaxmons(dtdef, DELk1, Tsk1, 5.*dtdef);
    cout << "T-Maximum attained at (E_0,L_0)=(" << ststk1.gete0() << ", " << ststk1.getL0() << "): " << Tmaxmonsk1.at(0) << endl;
    cout << "T_Max attained at boundary: " << Tmaxmonsk1.at(1) << endl;
    cout << "T increasing in E: " << Tmaxmonsk1.at(2) << endl;
    cout << "T decreasing in L: " << Tmaxmonsk1.at(3) << endl;
    cout << "GLref and HLref positive: " << GLrefHLrefposk1 << endl;

    cout << endl << "Scanning through near circular limits of partial_E T and partial_L T :" << endl;
    bool dETlimpos = true;
    bool dLTlimneg = true;
    for(int i = 0; i < Lsteps; i++)
    {
        double Lloc = ststk1.getLmax()*((double)i)/((double)Lsteps);
        if(ststk1.dETnearcirc(Lloc) < 0.)
            dETlimpos = false;
        if(ststk1.dLTnearcirc(Lloc) > 0.)
            dLTlimneg = false;
    }
    cout << "partial_E T limit always positive? " << dETlimpos << endl;
    cout << "partial_L T limit always negative? " << dLTlimneg << endl;

    cout << endl << endl << "Finished Period Function k=1 Part." << endl << endl;

    return;
}

