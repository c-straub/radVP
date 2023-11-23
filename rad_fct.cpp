/*
Source file for defining methods and constructor in RadialFct class.
*/

#include "rad_fct.h"


#include <cassert>
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


double RadFct::atr(double rloc)
{
    double r = rloc - rmin;
    double quot = r/dr;
    int indr = ceil(quot);
    double wtl = ((double)indr) - quot;
    double wtr = 1.-wtl;
    if(indr > 0 && indr < getnr())
        return wtl*atind(indr-1) + wtr*atind(indr);
    else if(indr == getnr())
    {
        cout << "Evaluation warning";
        return atind(indr-1);
    }
    else
        exit(7);
}

double RadFct::fctoverrsquareatr(double rloc)
{
    double r = rloc - rmin;
    double quot = r/dr;
    int indr = ceil(quot);
    double wtl = ((double)indr) - quot;
    double wtr = 1.-wtl;
    if(indr > 1)
        return wtl*((atind(indr-1)/(((double)(indr-1))*dr))/(((double)(indr-1))*dr)) + wtr*((atind(indr)/(((double)(indr))*dr))/(((double)(indr))*dr));
    else if(indr == 1)
        return wtr*(atind(1)/dr)/dr; //assuming that function over r^2 vanishes at r=0
    else
    {
 //       lock_guard<mutex> lock(muter);
        cout << endl << endl << "Evaluation error for fct over r*r; tried evaluation at r=" << r << " with index " << indr << endl << endl;
        exit(8);
    }
}

double RadFct::lpnorm(RadFct *other, double p, double rpower)
{
    if(p < numeric_limits<double>::epsilon())
        return 0.;

    double returner = 0.;

    for(int i = 0; i < getnr(); i++)
    {
        double rloc = getrmin() + getdr()*((double)i);

        returner += pow(rloc, rpower)*pow(fabs(atind(i)-other->atind(i)), p)*getdr();
    }

    return pow(returner, 1./p);
}
