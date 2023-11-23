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

#include "per_analyzer.h"

#include "steadystate.h" //For file writing things

vector<pair<double, double>> FileToVectorPair(string fname, int col1, int col2, int coltot,
                                              double shift2) //Saves values of columns col1 and col2 from input file as vector<pair>
{
    vector<pair<double, double>> returner;
    bool success = true;
    vector<double> inputter(coltot);

    ifstream infile(fname);

    double helper;

    //First skip lines starting with '#'
    string line;
    while(getline(infile, line))
    {
        if(line.empty() || line[0] == '#')
            continue; //Skip lines starting with #

        istringstream lineStream(line);

        for(int j = 0; j < coltot; j++)
        {
            if(!(lineStream >> inputter.at(j)))
                success = false;
        }

        pair<double, double> newelement(inputter.at(col1-1), inputter.at(col2-1)+shift2);
        returner.push_back(newelement);
    }
    infile.close();

    if(!success)
        cout << endl << "Warning: Failed converting file " << fname << "to vector<pair> :(" << endl;

    return returner;
}

double compPeriodSimple(vector<pair<double, double>> *dat, int returntype, string fname, bool withfile)
{
    vector<int> maxinds, mininds, zeroinds;

    double valmin2, valmin1, valnow;

    if(dat->size() <= 1)
        return -1.;

    //do first two steps explictly
    valmin1 = dat->at(0).second;
    valnow = dat->at(1).second;

    if((valmin1 < 0. && valnow > 0.) || (valnow < 0. && valmin1 > 0.))
        zeroinds.push_back(1);

    for(int i = 2; i < dat->size(); i++)
    {
        valmin2 = valmin1;
        valmin1 = valnow;
        valnow = dat->at(i).second;

        //Check for maxima and minima with three-point-method
        if(valmin2 >= valmin1 && valmin1 <= valnow)
            mininds.push_back(i-1); //found min
        if(valmin2 <= valmin1 && valmin1 >= valnow)
            maxinds.push_back(i-1); //found max

        //Search for zero with two-point-method
        if((valmin1 <= 0. && valnow >= 0.) || (valnow <= 0. && valmin1 >= 0.))
            zeroinds.push_back(i); //found zero
    }

    double l3p = 0.; //L3+, Period according to (at most) last 3 Maxima Distances (so consider 4 Maxima!)
    double l3m = 0.; //L3-, Period according to (at most) last 3 Minima Distances
    double l30 = 0.; //L30, Period according to (at most) last 3 Zeros Distances

    double l7p = 0.; //L7+, Period according to (at most) last 7 Maxima Distances
    double l7m = 0.; //L7-, Period according to (at most) last 7 Minima Distances
    double l70 = 0.; //L70, Period according to (at most) last 7 Zeros Distances

    double allp = 0.; //All+, Period according to all Maxima
    double allm = 0.; //All-, Period according to all Minima
    double all0 = 0.; //All0, Period according to all Zeros

    //Compute periods from Maxima data
    if(maxinds.size() <= 1)
    {
        l3p = -1.;
        l7p = -1.;
        allp = -1.;
    }
    else
    {
        int nrdifferences = maxinds.size()-1;

        int n3 = min(nrdifferences, 3);
        int n7 = min(nrdifferences, 7);

        double per3, per7, perall;
        for(int i = nrdifferences; i > nrdifferences-n3; i--)
            per3 +=  (dat->at(maxinds.at(i)).first-dat->at(maxinds.at(i-1)).first);
        per3 /= ((double)n3);

        for(int i = nrdifferences; i > nrdifferences-n7; i--)
            per7 +=  (dat->at(maxinds.at(i)).first-dat->at(maxinds.at(i-1)).first);
        per7 /= ((double)n7);

        for(int i = nrdifferences; i >= 1; i--)
            perall +=  (dat->at(maxinds.at(i)).first-dat->at(maxinds.at(i-1)).first);
        perall /= ((double)nrdifferences);

        l3p = per3;
        l7p = per7;
        allp = perall;
    }

    //Compute periods from Minima data
    if(mininds.size() <= 1)
    {
        l3m = -1.;
        l7m = -1.;
        allm = -1.;
    }
    else
    {
        int nrdifferences = mininds.size()-1;

        int n3 = min(nrdifferences, 3);
        int n7 = min(nrdifferences, 7);

        double per3, per7, perall;
        for(int i = nrdifferences; i > nrdifferences-n3; i--)
            per3 +=  (dat->at(mininds.at(i)).first-dat->at(mininds.at(i-1)).first);
        per3 /= ((double)n3);

        for(int i = nrdifferences; i > nrdifferences-n7; i--)
            per7 +=  (dat->at(mininds.at(i)).first-dat->at(mininds.at(i-1)).first);
        per7 /= ((double)n7);

        for(int i = nrdifferences; i >= 1; i--)
            perall +=  (dat->at(mininds.at(i)).first-dat->at(mininds.at(i-1)).first);
        perall /= ((double)nrdifferences);

        l3m = per3;
        l7m = per7;
        allm = perall;
    }

    //Compute periods from zeros data
    if(zeroinds.size() <= 1)
    {
        l30 = -1.;
        l70 = -1.;
        all0 = -1.;
    }
    else
    {
        int nrdifferences = zeroinds.size()-1;

        int n3 = min(nrdifferences, 3);
        int n7 = min(nrdifferences, 7);

        double per3, per7, perall;
        for(int i = nrdifferences; i > nrdifferences-n3; i--)
            per3 +=  (dat->at(zeroinds.at(i)).first-dat->at(zeroinds.at(i-1)).first);
        per3 /= ((double)n3);

        for(int i = nrdifferences; i > nrdifferences-n7; i--)
            per7 +=  (dat->at(zeroinds.at(i)).first-dat->at(zeroinds.at(i-1)).first);
        per7 /= ((double)n7);

        for(int i = nrdifferences; i >= 1; i--)
            perall +=  (dat->at(zeroinds.at(i)).first-dat->at(zeroinds.at(i-1)).first);
        perall /= ((double)nrdifferences);

        l30 = 2.*per3;
        l70 = 2.*per7;
        all0 = 2.*perall;
    }

    double returner;
    if(returntype == 1)
        returner = l3p;
    else if(returntype == 2)
        returner = l3m;
    else if(returntype == 3)
        returner = l30;
    else if(returntype == 4)
        returner = l7p;
    else if(returntype == 5)
        returner = l7m;
    else if(returntype == 6)
        returner = l70;
    else if(returntype == 7)
        returner = allp;
    else if(returntype == 8)
        returner = allm;
    else if(returntype == 9)
        returner = all0;
    else
        returner = -.1;

    if(withfile)
    {
        stringstream ss; //Used for output

        ss << returner << "\n";

        ss << "\n" << "Overview over evaluated periods:" << "\n\n";
        ss << "L3+ = " << l3p << ",\n";
        ss << "L3- = " << l3m << ",\n";
        ss << "L30 = " << l30 << ";\n\n";
        ss << "L7+ = " << l7p << ",\n";
        ss << "L7- = " << l7m << ",\n";
        ss << "L70 = " << l70 << ";\n\n";
        ss << "All+ = " << allp << ",\n";
        ss << "All- = " << allm << ",\n";
        ss << "All0 = " << all0 << ".\n\n";

        ss << "\n\n Raw Data: " << "\n\n\n\n";

        ss << "Maximas:\n\n";
        for(int i = 0; i < maxinds.size(); i++)
        {
            ss << "At " << dat->at(maxinds.at(i)).first << ", value " << dat->at(maxinds.at(i)).second;
            if(i > 0)
                ss << ", distance to previous max " << dat->at(maxinds.at(i)).first - dat->at(maxinds.at(i-1)).first;
            ss << ".\n";
        }

        ss << "\n\n\n Minimas:\n\n";
        for(int i = 0; i < mininds.size(); i++)
        {
            ss << "At " << dat->at(mininds.at(i)).first << ", value " << dat->at(mininds.at(i)).second;
            if(i > 0)
                ss << ", distance to previous min " << dat->at(mininds.at(i)).first - dat->at(mininds.at(i-1)).first;
            ss << ".\n";
        }

        ss << "\n\n\n Zeros:\n\n";
        for(int i = 0; i < zeroinds.size(); i++)
        {
            ss << "At " << dat->at(zeroinds.at(i)).first;
            if(i > 0)
                ss << ", distance to previous zero " << dat->at(zeroinds.at(i)).first - dat->at(zeroinds.at(i-1)).first << " = 0.5 * " << 2.*(dat->at(zeroinds.at(i)).first - dat->at(zeroinds.at(i-1)).first);
            ss << ".\n";
        }

        ss << "\n\n";

        clearOutputFile(fname, "out_per_simple");
        appendSSToOutputFile(fname, "out_per_simple", &ss);
    }

    return returner;
}

double compPeriodSimple(vector<pair<double, double>> *dat, int returntype)
{
    return compPeriodSimple(dat, returntype, "", false);
}

double compPeriodSimple(vector<pair<double, double>> *dat, string fname, int returntype)
{
    return compPeriodSimple(dat, returntype, fname, true);
}

double periodSimple_FromFile(string fnamein, int col1, int col2, int coltot, string fnameout, double shift2, int returntype)
{
    vector<pair<double, double>> perdat = FileToVectorPair(fnamein, col1, col2, coltot, shift2);
    return compPeriodSimple(&perdat, fnameout, returntype);
}



bool LowerPairFirst(const pair<double,double> &a,const pair<double,double> &b) //Helper function for sorting
{
       return a.first<b.first;
}

bool LargerPairSecond(const pair<double,double> &a,const pair<double,double> &b) //Helper function for sorting
{
       return a.second>b.second;
}

bool LargerDouble(const double &a, const double &b) //Helper function for sorting
{
       return a>b;
}

vector<double> dftNanalyse(vector<pair<double, double>> dat, string fname, double tstart)
{
    return dftNanalyse(dat, 10, fname, tstart);
}

vector<double> dftNanalyse(vector<pair<double, double>> dat, string fname)
{
    return dftNanalyse(dat, 10, fname);
}

vector<double> dftNanalyse(vector<pair<double, double>> dat, int initperiodgrid, string fname, double tstart)
{
    vector<pair<double,double>> newdat;
    for(int i=0; i < dat.size(); i++)
        if(dat.at(i).first >= tstart)
            newdat.push_back(dat.at(i));
    return dftNanalyse(newdat, initperiodgrid, fname);
}

vector<double> dftNanalyse(vector<pair<double, double>> dat, int initperiodgrid, string fname)
{
    int N = dat.size();
    double dt = dat.at(1).first - dat.at(0).first;

    //First go through data with initialperiodgrid size and compute amplitude of modes using DFT
    vector<pair<double, double>> initamps; //stores frequency and amplitude
    int initbound = ceil(((double)N)/((double)initperiodgrid));


    initamps.push_back(pair<double, double>(0, 0.));
    for(int i = 1; i < initbound; i++)
    {
        double realpart = 0.;
        double impart = 0.;

        for(int j = 0; j < N; j++)
        {
            realpart += dat.at(j).second*cos(2.*M_PI*((double)(i))*((double)j)/((double)N));
            impart -= dat.at(j).second*sin(2.*M_PI*((double)(i))*((double)j)/((double)N));
        }

        double amp = sqrt(realpart*realpart + impart*impart)/((double)N);

        pair<double, double> helper(i, amp);
        initamps.push_back(helper);
    }

    //Compute average amp
    int initsize = initamps.size();
    double initave = 0.;
    for(int i = 0; i < initsize; i++)
        initave += initamps.at(i).second;
    initave /= ((double)initsize);

    //Detect regions of interest, i.e., regions, where the amp is higher than 50 times the average value
    double interestfactor = 50.;
    double interestthreshold = interestfactor*initave;
    vector<pair<double, double>> interestregions;
    bool ininterestregion = false;
    double intereststarter = 0.;
    for(int i = 0; i < initsize; i++)
    {
        if((!ininterestregion) && initamps.at(i).second > interestthreshold)
        { // Start new interest region
            ininterestregion = true;
            intereststarter = initamps.at(i).first;
        }
        else if(ininterestregion && initamps.at(i).second < interestthreshold)
        { //end current interest region
            ininterestregion = false;
            pair<double,double> helper(intereststarter, initamps.at(i-1).first);
            interestregions.push_back(helper);
        }
    }
    if(ininterestregion)
    { //manually end current interest region
        pair<double,double> helper(intereststarter, initamps.at(N-1).first);
        if(initamps.at(N-1).second < interestthreshold && N >= 2)
            helper.second = initamps.at(N-2).first;
        interestregions.push_back(helper);
    }


/*    cout << "Found following interest regions; threshold is " << interestthreshold << endl;
    for(int i = 0; i < interestregions.size(); i++)
        cout << "From i=" << interestregions.at(i).first << " to " << interestregions.at(i).second << "." << endl;
    cout << endl << endl;*/

    if(interestregions.size() <= 0)
    { //If no interest regions have been found -> take maximizer
        int maxind = initamps.at(0).first;
        double maxamp = initamps.at(0).second;
        for(int i = 1; i < initamps.size(); i++)
        {
            if(initamps.at(i).second > maxamp)
            {
                maxind = initamps.at(i).first;
                maxamp = initamps.at(i).second;
            }
        }
        pair<double, double> helper(maxind, maxind);
        interestregions.push_back(helper);
    }

    double moreaccuratefactor = 500.;

    vector<double> returner;
    vector<double> returnervals;

    for(int i = 0; i < interestregions.size(); i++)
    {
        vector<pair<double, double>> intregamps; //amps in current interest region


        if(interestregions.at(i).first >= interestregions.at(i).second - numeric_limits<double>::epsilon())
        { //enlarge interestregion if it is only a singleton
            ;//interestregions.at(i).first = max(interestregions.at(i).first-1.,0.);
            //interestregions.at(i).second = min(interestregions.at(i).second+1.,N-1.);
            double periodhelper = dt*((double)N)/interestregions.at(i).first;
            returner.push_back(periodhelper);
            returnervals.push_back(-1);
            continue;
        }

        double frequhelper = interestregions.at(i).first; //+ 1./moreaccuratefactor;
        if(frequhelper < numeric_limits<double>::epsilon())
            frequhelper += 1./moreaccuratefactor;

        while(frequhelper < interestregions.at(i).second)
        {
            double realpart = 0.;
            double impart = 0.;

            for(int j = 0; j < N; j++)
            {
                realpart += dat.at(j).second*cos(2.*M_PI*frequhelper*((double)j)/((double)N));
                impart -= dat.at(j).second*sin(2.*M_PI*frequhelper*((double)j)/((double)N));
            }

            double amp = sqrt(realpart*realpart + impart*impart)/((double)N);

            pair<double, double> helper(frequhelper, amp);
            intregamps.push_back(helper);

            frequhelper += 1./moreaccuratefactor;
        }

        //Find largest amp in current interest region
        int maxind = 0;
        double maxval = intregamps.at(0).second;
        for(int j = 1; j < intregamps.size(); j++)
            if(intregamps.at(j).second > maxval)
            {
                maxind = j;
                maxval = intregamps.at(j).second;
            }

        //cout << "Largest amp in " << i << "th interest region is " << intregamps.at(maxind).second << endl;


        double periodhelper = dt*((double)N)/intregamps.at(maxind).first;

        returner.push_back(periodhelper);
        returnervals.push_back(intregamps.at(maxind).second);

        for(int j = 0; j < intregamps.size(); j++)
            initamps.push_back(intregamps.at(j));
    }

    sort(initamps.begin(),initamps.end(),LowerPairFirst);

    sort(returner.begin(), returner.end(), LargerDouble);

    stringstream ssana;

    ssana << "# Average: " << initave << "\n#\n";
    ssana << "# Found following periods: \n";
    for(int i = 0; i < returner.size(); i++)
        ssana << "# " << returner.at(i) << ", Fourier Amplitude is " << returnervals.at(i) <<  "\n";
    ssana << "# \n";
    ssana << "# Now raw data:\n#\n";

    for(int i = 0; i < initamps.size(); i++)
        ssana << initamps.at(i).first << " " << dt*((double)N)/initamps.at(i).first << " " << initamps.at(i).second << "\n";

    clearOutputFile(fname, "out_per_dft");
    appendSSToOutputFile(fname, "out_per_dft", &ssana);


    return returner;

}

vector<pair<double, double>> dft(vector<pair<double, double>> dat)
{
    int N = dat.size();

    vector<pair<double, double>> returner(N);

    for(int i = 0; i < N; i++)
    {
        double realpart = 0.;
        double impart = 0.;

        for(int j = 0; j < N; j++)
        {
            realpart += dat.at(j).second*cos(2.*M_PI*((double)i)*((double)j)/((double)N));
            impart -= dat.at(j).second*sin(2.*M_PI*((double)i)*((double)j)/((double)N));
        }

        pair<double, double> helper(realpart, impart);
        returner.at(i) = helper;
    }

    return returner;
}

void analyseDft(vector<pair<double, double>> dftdat, string fname, double dt)
{
    stringstream ss; //stringstream for raw data

    ss << "# Raw Data of DFT\n";
    ss << "#\n#\n# FORMAT\n# 1) Index of Period |2) Period |3) Real Part of Fourier Coeff. |4) Imaginary Part |5) Amplitude, i.e., normalized modulus of Fourier Coeff. |6) Phase\n";
    ss << "#\n#\n";

    int N = dftdat.size();

    for(int i = 0; i < N; i++)
    {
        double amp = sqrt(dftdat.at(i).first*dftdat.at(i).first + dftdat.at(i).second*dftdat.at(i).second)/((double)N);
        double phase = atan2(dftdat.at(i).second,dftdat.at(i).first);

        ss << i << " " << ((double)i)*dt << " " << dftdat.at(i).first << " " << dftdat.at(i).second << " " << amp << " " << phase << "\n";

    }

    clearOutputFile(fname, "out_dft_rawdat");
    appendSSToOutputFile(fname, "out_dft_rawdat", &ss);

    return;
}

vector<double> dftNanalyse_FromFile(string fnamein, int col1, int col2, int coltot,
                                    string fnameout, double tstart, double shift2)
{
    vector<pair<double, double>> perdat = FileToVectorPair(fnamein, col1, col2, coltot, shift2);
    vector<double> dftdat = dftNanalyse(perdat, fnameout, tstart);
    return dftdat;
}


bool detect_damping(vector<pair<double, double>> *dat, string fnameout, int gfactor)
{
    //First compute g_+, g_-,g_abs from given data on larger grid given by gfactor

    int indhelper = 0;

    vector<pair<double, double>> gplusdat;
    vector<pair<double, double>> gminusdat;
    vector<pair<double, double>> gabsdat;

    while(indhelper < dat->size())
    {
        double tloc = dat->at(indhelper).first;

        gplusdat.push_back(pair<double, double>(tloc, gplus(dat, tloc)));
        gminusdat.push_back(pair<double, double>(tloc, gminus(dat, tloc)));
        gabsdat.push_back(pair<double, double>(tloc, gabs(dat, tloc)));

        indhelper += gfactor;
    }

    ///TODO !


    stringstream ssout;

    //Write g-data into file

    ssout << "#\n# Now g-data. FORMAT: \n# 1) t |2) g_+ |3) g_- |4) g_abs\n#\n#\n";

    for(int i = 0; i < gabsdat.size(); i++)
        ssout << gabsdat.at(i).first << " " << gplusdat.at(i).second << " "
              << gminusdat.at(i).second << " " << gabsdat.at(i).second << "\n";


    clearOutputFile(fnameout, "out_per_damp");
    appendSSToOutputFile(fnameout, "out_per_damp", &ssout);


    ///TODO
    return true;


}

bool detect_damping_FromFile(string fnamein, int col1, int col2, int coltot,
                             string fnameout, int gfactor)
{
    vector<pair<double, double>> datain = FileToVectorPair(fnamein, col1, col2, coltot);
    return detect_damping(&datain, fnameout, gfactor);
}

double gplus(vector<pair<double, double>> *dat, double t)
{
    return max_ininterval(dat, t, numeric_limits<double>::max()).second;
}

double gminus(vector<pair<double, double>> *dat, double t)
{
    return min_ininterval(dat, t, numeric_limits<double>::max()).second;
}

double gabs(vector<pair<double, double>> *dat, double t)
{
    return (gplus(dat, t) - gminus(dat,t));
}

pair<double,double> max_ininterval(vector<pair<double, double>> *dat, double left, double right)
{
    double maxloc = numeric_limits<double>::lowest();
    int maxind = -1;
    for(int i = 0; i < dat->size(); i++)
        if(dat->at(i).first >= left && dat->at(i).first <= right && dat->at(i).second > maxloc)
        {
            maxloc = dat->at(i).second;
            maxind = i;
        }

    if(maxind < 0)
        return pair<double,double>(0.,0.);


    return dat->at(maxind);
}

pair<double,double> min_ininterval(vector<pair<double, double>> *dat, double left, double right)
{
    double minloc = numeric_limits<double>::max();
    int minind = -1;
    for(int i = 0; i < dat->size(); i++)
        if(dat->at(i).first >= left && dat->at(i).first <= right && dat->at(i).second < minloc)
        {
            minloc = dat->at(i).second;
            minind = i;
        }

    if(minind < 0)
        return pair<double,double>(0.,0.);


    return dat->at(minind);
}

double maxfirst(vector<pair<double, double>> *dat)
{
    double maxloc = numeric_limits<double>::lowest();
    for(int i = 0; i < dat->size(); i++)
        if(dat->at(i).first > maxloc)
            maxloc = dat->at(i).first;
    return maxloc;
}

double minfirst(vector<pair<double, double>> *dat)
{
    double minloc = numeric_limits<double>::max();
    for(int i = 0; i < dat->size(); i++)
        if(dat->at(i).first < minloc)
            minloc = dat->at(i).first;
    return minloc;
}

double meansecond(vector<pair<double, double>> *dat, double tstart)
{
    double meanloc = 0.;
    int counterhelp = 0;
    for(int i = 0; i < dat->size(); i++)
        if(dat->at(i).first > tstart)
        {
            meanloc += dat->at(i).second;
            counterhelp++;
        }
    if(counterhelp == 0)
        return 0.;
    else
        return meanloc/((double)counterhelp);
}

double meansecond_FromFile(string fnamein, int col1, int col2, int coltot, double tstart)
{
    vector<pair<double, double>> perdat = FileToVectorPair(fnamein, col1, col2, coltot);
    return meansecond(&perdat, tstart);
}
