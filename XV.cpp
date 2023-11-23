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

#include "XV.h"



XV::XV(double ri, double wi, double Li)
{
    if(ri < numeric_limits<double>::epsilon())
    {
        if(Li < numeric_limits<double>::epsilon())
        {
            x1=0.; x2=0.; x3=0.;
            v1=wi; v2=0.; v3=0.;
        }
        else
        {
            cout << endl << endl
                 << "Error; trying to create (x,v) for (r,w,L) with r=0 < L=" << Li << "..."
                 << endl << endl;
            exit(11);
        }
    }
    else
    {
        x1 = ri; x2 = 0.; x3 = 0.;
        v1 = wi; v2 = sqrt(Li)/ri; v3 = 0.;
    }
}


double XV::comp_w() //Computes w=(x/r) \cdot v
{
    return comp_w(comp_r());
}

double XV::comp_w(double r)
{
    if(r < numeric_limits<double>::epsilon())
    {
           cout << endl << endl
                << "Error; trying to compute w=(x dot v)/r for x=0 ?!"
                << endl << endl;
            exit(12);
    }
    return (x1/r)*v1 + (x2/r)*v2 + (x3/r)*v3;
}

void XV::rk4_update(XV k1, XV k2, XV k3, XV k4, double dt) //Updates (x,v) vector according to RK4
{
    x1 += dt*(k1.x1 + 2.*k2.x1 + 2.*k3.x1 + k4.x1)/6.;
    x2 += dt*(k1.x2 + 2.*k2.x2 + 2.*k3.x2 + k4.x2)/6.;
    x3 += dt*(k1.x3 + 2.*k2.x3 + 2.*k3.x3 + k4.x3)/6.;

    v1 += dt*(k1.v1 + 2.*k2.v1 + 2.*k3.v1 + k4.v1)/6.;
    v2 += dt*(k1.v2 + 2.*k2.v2 + 2.*k3.v2 + k4.v2)/6.;
    v3 += dt*(k1.v3 + 2.*k2.v3 + 2.*k3.v3 + k4.v3)/6.;
}


XV addWithFactor(const XV xv, double factor, const XV adder)
{
    return XV(xv.x1 + factor*adder.x1, xv.x2 + factor*adder.x2, xv.x3 + factor*adder.x3,
              xv.v1 + factor*adder.v1, xv.v2 + factor*adder.v2, xv.v3 + factor*adder.v3);
}

