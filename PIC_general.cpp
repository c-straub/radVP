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
#include <random>

#include "PIC_general.h"
#include "steadystate.h" //To use clearOutputFile and appendSSToOutputFile functions

void Outers::reset()
{
    clearAllStreams();

    clearOutputFile(fname, prefix_timecontroller);
    clearOutputFile(fname, prefix_globout);
    clearOutputFile(fname, prefix_rhomout);
    clearOutputFile(fname, prefix_trajout);
    clearOutputFile(fname, prefix_chessout);
}

void Outers::clearAllStreams()
{
    timecontroller.str(string());
    globout.str(string());
    rhomout.str(string());
    trajout.str(string());

    return;
}

void Outers::setrajinds(int am_trajparts, int numpart)
{
    trajinds.clear();

    if(am_trajparts > 0)
    {
        int trajfactor = numpart/am_trajparts;
        trajfactor++;

        while(!isPrime(trajfactor)) //We can sleep more peacefully at night if the indices are separated by prime numbers. Whether this actually has any impact, we do not know.
            trajfactor++;

        for(int i = 0; i < numpart; i += trajfactor)
            trajinds.push_back(i);

        int i = trajfactor/2;
        while(trajinds.size() < am_trajparts && i < numpart)
        {
            trajinds.push_back(i);
            i += trajfactor;
        }
    }
}

void Outers::writeConsole(string consoletext)
{
    cout << consoletext;

    timecontroller << consoletext;
    appendSSToOutputFile(fname, prefix_timecontroller, &timecontroller);
}

void Outers::writeConsoleAndClean(stringstream *consoletext)
{
    writeConsole(consoletext->str());
    consoletext->str(string());
}

void Outers::writeAllFiles()
{
    appendSSToOutputFile(fname, prefix_trajout, &trajout);
    appendSSToOutputFile(fname, prefix_rhomout, &rhomout);
    appendSSToOutputFile(fname, prefix_globout, &globout);
    stringstream chessss = chessout.outputstring();
    clearOutputFile(fname, prefix_chessout);
    appendSSToOutputFile(fname, prefix_chessout, &chessss);
}

void Outers::initwrite_rhomout(PICParams pp, Fields *fields, double rpmin, double rpmax)
{
    rhomout << "# Number of plotted radii: " << pp.Nrplot << "\n";
    rhomout << "# Minimal plotted radius: " << rpmin << "\n";
    rhomout << "# Maximal plotted radius: " << rpmax << "\n";
    rhomout << "# drplot: " << ((rpmax-rpmin)/((double)pp.Nrplot)) << "\n";
    rhomout << "# Minimal radius of fields: " << fields->getrmin() << "\n";
    rhomout << "# Maximal radius of fields: " << fields->getrmax() << "\n";
    rhomout << "# dr of fields: " << fields->getdr() << "\n";
    rhomout << "# Nr. of radial steps in fields: " << fields->getnr() << "\n";
    rhomout << "# Real t-step size: " << pp.dtoverdr*pp.dr << "\n";
    rhomout << "# Plot t-step size: " << pp.T/((double)pp.am_tplots) << "\n";
    rhomout << "# Nr. of t plot steps: " << pp.am_tplots << "\n";
    rhomout << "#\n#\n";
    rhomout << "# FORMAT: 1) r |2) rho(t,r) |3) m(t,r) |4) pot(t,r)\n";
    rhomout << "#\n#\n";
}

void Outers::write_rhomout(PICParams pp, double t, Fields *fields, double rpmin, double rpmax)
{
    if(pp.Nrplot <= 1)
        return;
    double drplot = ((rpmax-rpmin)/((double)(pp.Nrplot-1)));
    rhomout << "t=" << t << "\n";
    for(int i = 0; i < pp.Nrplot; i++)
    {
        double rloc = rpmin + ((double)i)*drplot;
        if(rloc-fields->getrmin() < numeric_limits<double>::epsilon())
            rhomout << rloc << " " << fields->rho.atind(0) << " " << fields->mass.atind(0) << " " << fields->pot.atind(0) << "\n";
        else
            rhomout << rloc << " " << fields->rho.atr(rloc) << " " << fields->mass.atr(rloc) << " " << fields->pot.atr(rloc) << "\n";
    }
    rhomout << "\n\n";
}

void Outers::initwrite_trajout(PICParams pp, int am_part)
{
    trajout << "# Tracked Particles: " << trajinds.size() << "\n";
    trajout << "# Total Particles: " << am_part << "\n";
    trajout << "# Portion of tracked particles: " << 100.*((double)trajinds.size())/((double)am_part) << " %\n";
    rhomout << "# Real t-step size: " << pp.dtoverdr*pp.dr << "\n";
    rhomout << "# Plot t-step size: " << pp.T/((double)pp.am_tplots) << "\n";
    trajout << "# Nr. of t plot steps: " << pp.am_tplots << "\n";
    trajout << "#\n#\n";
    trajout << "# FORMAT: 1) r |2) w |3) L |4) E |5) f(r,w,L) |6) volume of cell |7) E(r,w,L)\n";
    trajout << "#\n#\n";
}

void Outers::write_trajout(double t, vector<ParPos> *parpos)
{
    trajout << "t=" << t << "\n";
    for(int i = 0; i < trajinds.size(); i++)
    {
        ParPos parloc = parpos->at(trajinds.at(i));
        trajout << parloc.r << " " << parloc.w << " " << parloc.L << " "
                << parloc.E << " " << parloc.f << " " << parloc.vol << "\n";
    }
    trajout << "\n\n";
}

void ParPos::serialise(ofstream& stream) const
{
    stream.write(reinterpret_cast<const char*>(&r), sizeof(r));
    stream.write(reinterpret_cast<const char*>(&w), sizeof(w));
    stream.write(reinterpret_cast<const char*>(&L), sizeof(L));
    stream.write(reinterpret_cast<const char*>(&E), sizeof(E));
    stream.write(reinterpret_cast<const char*>(&f), sizeof(f));
    stream.write(reinterpret_cast<const char*>(&vol), sizeof(vol));
}

void ParPos::deserialise(istream& stream)
{
    stream.read(reinterpret_cast<char*>(&r), sizeof(r));
    stream.read(reinterpret_cast<char*>(&w), sizeof(w));
    stream.read(reinterpret_cast<char*>(&L), sizeof(L));
    stream.read(reinterpret_cast<char*>(&E), sizeof(E));
    stream.read(reinterpret_cast<char*>(&f), sizeof(f));
    stream.read(reinterpret_cast<char*>(&vol), sizeof(vol));
    deletepart = false;
}

void Outers::write_partoutFile(vector<ParPos> *parpos)
{
    string fnametxt = prefix_partout + "_" + fname + ".bin";
    ofstream partouts(fnametxt, ios::binary);
    for(int kp = 0; kp < parpos->size(); kp++)
        parpos->at(kp).serialise(partouts);
    partouts.close();
}

vector<ParPos> init_FromFile(string fnamesuff)
{
    string fnametxt = "outPART_" + fnamesuff + ".bin";
    vector<ParPos> returner;
    ifstream partins(fnametxt, ios::binary);
    while(partins)
    {
        ParPos parposloc;
        parposloc.deserialise(partins);
        parposloc.deletepart = false;
        returner.push_back(parposloc);
    }
    partins.close();
    return returner;
}

string conTabularRow(vector<string> row)
{
    stringstream returnerss;

    for(int i = 0; i < row.size(); i++)
        returnerss << " " << setw(11) << row.at(i) << " |";
    returnerss << "\n";

    return returnerss.str();
}

string conTabularSeparatorRow(int columns)
{
    stringstream returnerss;

    for(int i = 0; i < columns; i++)
        returnerss << setw(14) << setfill('_') << "_|";
    returnerss  << setfill(' ') << "\n";

    return returnerss.str();
}

string time_output(double secs)
{
    stringstream returner;

    if(secs < 60.)
        returner << setprecision(3) << secs << " s";
    else if(secs < 3600.)
        returner << setprecision(3) << secs/60. << " min";
    else if(secs < 86400.)
        returner << setprecision(3) << secs/3600. << " h";
    else
        returner << setprecision(3) << secs/86400. << " d";

    return returner.str();
}

string to_string_precision(double x, int prec)
{
    stringstream returner;
    returner << setprecision(prec) << x;
    return returner.str();
}

void Fields::addSingleParticle_NOVOL(ParPos *parpos)
{
    addSingleParticle_NOVOL_moms(parpos, &rho, 0);
}

void Fields::comprho_seq(vector<ParPos> *parpos)
{
    rho.setallzero();

    for(int kp = 0; kp < parpos->size(); kp++) //First compute values without volume element
        addSingleParticle_NOVOL(&parpos->at(kp));

    //Adjust volume element
    double dvolloc = 4.*M_PI*rho.getdr();
    rho.divall(dvolloc);
}

void Fields::compmass()
{
    mass.setallzero();
    double dr = getdr();
    int nrmax = mass.getnr()-1;

    mass.setind(0, 0.);

    double rloc1 = mass.getrmin() + dr;
    double rloc2 = mass.getrmin() + 2*dr;


    mass.setind(1, 2.*M_PI*dr*rloc1*rho.atind(1) - .5*M_PI*dr*rloc2*rho.atind(2)/3.);
    for(int kr = 2; kr < nrmax; kr++)
    {
        double rlocm2 = mass.getrmin() + ((double)(kr-2))*dr;
        double rlocm1 = mass.getrmin() + ((double)(kr-1))*dr;
        double rloc = mass.getrmin() + ((double)(kr))*dr;
        double rlocp1 = mass.getrmin() + ((double)(kr+1))*dr;
        mass.addtoprevind(kr, 4.*M_PI*dr*( - 1.*(rlocm2*rlocm2*rho.atind(kr-2) + rlocp1*rlocp1*rho.atind(kr+1))
                                       + 13.*(rlocm1*rlocm1*rho.atind(kr-1) + rloc*rloc*rho.atind(kr)) )/24.);
    }

    double rlocm2 = mass.getrmin() + ((double)(nrmax-2))*dr;
    double rlocm1 = mass.getrmin() + ((double)(nrmax-1))*dr;
    double rloc = mass.getrmin() + ((double)(nrmax))*dr;
    //Make sure that every rho appears in the overall mass with (added) weight 1.
    mass.addtoprevind(nrmax, 4.*M_PI*dr*(-rlocm2*rlocm2*rho.atind(nrmax-2)
                                         + 12.*rlocm1*rlocm1*rho.atind(nrmax-1)
                                         + 25.*rloc*rloc*rho.atind(nrmax))/24.);

    //Now set potprime:
    potprime.setallzero();
    potprime.setind(0, 0.);
    for(int kr = 1; kr < mass.getnr(); kr++)
    {
        double rloc = mass.getrmin() + ((double)(kr))*dr;
        potprime.setind(kr, (mass.atind(kr)/rloc)/rloc);
    }
}

double Fields::getmtot()
{
    return mass.get_last_entry();
}

double Fields::getEpot()
{
    double dr = mass.getdr();
    int nrmax = mass.getnr()-1;

    double rmax = dr*((double)nrmax) + mass.getrmin();
    double mtot = mass.atind(nrmax);
    double epot = -.5*(mtot*mtot)/rmax; //potential energy outside of radius rmax

    for(int kr = nrmax-1; kr >= 1; kr--)
    {
        double rloc0 = ((double)kr)*dr + mass.getrmin();
        epot -= .5*((mass.atind(kr)/rloc0)*(mass.atind(kr)/rloc0))*dr;
    }

    return epot;
}

double Fields::getEkin_seq(vector<ParPos> *parpos)
{
    double returner = 0.;
    for(int kp = 0; kp < parpos->size(); kp++)
        addSingleParticle_num(&parpos->at(kp), &returner, 2);

    return .5*returner;
}



void Fields::comppot()
{
    pot.setallzero();
    double dr = mass.getdr();
    int nrmax = mass.getnr()-1;

    double rmax = dr*((double)nrmax) + mass.getrmin();
    double mtot = mass.atind(nrmax);

    pot.setind(nrmax, -mtot/rmax);
    for(int kr = nrmax-1; kr >= 1; kr--)
    {
        double rloc = ((double)kr)*dr + mass.getrmin();
        double rlocp1 = ((double)(kr+1))*dr + mass.getrmin();

        pot.addtonextind(kr, -.5*dr*(mass.atind(kr)/(rloc*rloc) +  mass.atind(kr+1)/(rlocp1*rlocp1)));
    }
    double rloc = dr + mass.getrmin();
    pot.addtonextind(0, -.5*dr*mass.atind(1)/(rloc*rloc));
}

double Fields::getpot0()
{
    return pot.atind(0);
}

void Locs::setpartpartition(int am_parts)
{
    partpartition.clear();

    if(((double)am_parts)/((double)number_of_threads) <= 1.)
    {
        cout << endl << endl << "Error! There are too few particles or too much threads..." << endl;
        exit(563103);
    }
    for(int kthr = 0; kthr < number_of_threads; kthr++)
        partpartition.push_back(floor((((double)am_parts)/((double)number_of_threads))*((double)kthr)));
    partpartition.push_back(am_parts); //Last thread takes the rest
}

void Locs::setrpartition(int fieldslen)
{
    rpartition.clear();
    for(int kthr = 0; kthr < number_of_threads; kthr++)
        rpartition.push_back(floor(((double)(kthr*fieldslen))/((double)number_of_threads)));
    rpartition.at(0) = 0; //Make sure that nothing nasty happens at 0
    rpartition.push_back(fieldslen); //Last thread takes the rest
}

void Locs::comprho_parahelper(int tid, const vector<ParPos> *parpos)
{
    rholoc.at(tid).setallzero();
    //Compute rho using the batch of particles assigned to this thread
    for(int kp = partpartition.at(tid); kp < partpartition.at(tid+1); kp++)
        addSingleParticle_NOVOL_moms(&parpos->at(kp), &rholoc.at(tid), 0);
}

void Locs::mergerholocs_parahelper(int tid, RadFct *rho)
{
    for(int kr = rpartition.at(tid); kr < rpartition.at(tid+1); kr++) //Sets rho in radial part of the thread to zero
        rho->setind(kr, 0.);

    double overdvolloc = 4.*M_PI*rho->getdr();

    for(int kr = rpartition.at(tid); kr < rpartition.at(tid+1); kr++)
    {
        for(int kthr = 0; kthr < number_of_threads; kthr++)
            rho->addind(kr, rholoc.at(kthr).atind(kr));

        rho->divind(kr, overdvolloc); //divide radial volume element, as it is already contained in whole volume element
    }

    if(rpartition.at(tid) == 0)
        rho->multind(0, 2.);
}

void Locs::compEkin_parahelper(int tid, const vector<ParPos> *parpos)
{
    ekinloc.at(tid) = 0.;
    //Compute Ekin using the batch of particles assigned to this thread
    for(int kp = partpartition.at(tid); kp < partpartition.at(tid+1); kp++)
        addSingleParticle_num(&parpos->at(kp), &ekinloc.at(tid), 2);
    ekinloc.at(tid) *= .5;
}

double Locs::get_Ekin_all()
{
    double returner = 0.;
    for(int i = 0; i < ekinloc.size(); i++)
        returner += ekinloc.at(i);
    return returner;
}

void Locs::compEpot_parahelper(int tid, const vector<ParPos> *parpos, Fields *fields)
{
    epotloc.at(tid) = 0.;
    //Compute Epot using the batch of particles assigned to this thread
    for(int kp = partpartition.at(tid); kp < partpartition.at(tid+1); kp++)
        addSingleParticle_fac(&parpos->at(kp), &epotloc.at(tid), &fields->pot);
    epotloc.at(tid) *= .5;
}

double Locs::get_Epot_all()
{
    double returner = 0.;
    for(int i = 0; i < epotloc.size(); i++)
        returner += epotloc.at(i);
    return returner;
}

void Locs::rminmaxparahelper(int tid, const vector<ParPos> *parpos)
{
    double rminloc = parpos->at(partpartition.at(tid)).r;
    double rmaxloc = rminloc;

    for(int kp = partpartition.at(tid)+1; kp < partpartition.at(tid+1); kp++)
    {
        double rloc = parpos->at(kp).r;
        rminloc = min(rminloc, rloc);
        rmaxloc = max(rmaxloc, rloc);
    }

    rminmaxloc.at(tid).first = rminloc;
    rminmaxloc.at(tid).second = rmaxloc;
}

pair<double,double> Locs::getrminmax_all()
{
    double rminloc = rminmaxloc.at(0).first;
    double rmaxloc = rminmaxloc.at(0).second;

    for(int i = 1; i < rminmaxloc.size(); i++)
    {
        rminloc = min(rminloc, rminmaxloc.at(i).first);
        rmaxloc = max(rmaxloc, rminmaxloc.at(i).second);
    }

    return pair<double,double>(rminloc,rmaxloc);
}

void Locs::wminmaxparahelper(int tid, const vector<ParPos> *parpos)
{
    double wminloc = parpos->at(partpartition.at(tid)).w;
    double wmaxloc = wminloc;

    for(int kp = partpartition.at(tid)+1; kp < partpartition.at(tid+1); kp++)
    {
        double wloc = parpos->at(kp).w;
        wminloc = min(wminloc, wloc);
        wmaxloc = max(wmaxloc, wloc);
    }

    wminmaxloc.at(tid).first = wminloc;
    wminmaxloc.at(tid).second = wmaxloc;
}

void Locs::deleteparts(int tid, vector<ParPos> *parpos)
{
    nrdelpartsloc.at(tid) = 0;
    for(int kp = partpartition.at(tid)+1; kp < partpartition.at(tid+1); kp++)
    {
        if(parpos->at(kp).deletepart)
        {
            parpos->at(kp).f = 0.;

            int ko = rand()%parpos->size(); //Index of some other particle
            while(parpos->at(ko).deletepart) //Make sure that other particle is not among deleted ones
            {
                ko++;
                ko = ko%parpos->size();
            }

            parpos->at(kp).r = parpos->at(ko).r;
            parpos->at(kp).w = parpos->at(ko).w;
            parpos->at(kp).L = parpos->at(ko).L;
            parpos->at(kp).E = parpos->at(ko).E;
            parpos->at(kp).vol = parpos->at(ko).vol;

            parpos->at(kp).deletepart = false;

            nrdelpartsloc.at(tid)++;
        }
    }
}

void Locs::nrdelpartstot_update()
{
    for(int kthr = 0; kthr < number_of_threads; kthr++)
    {
        nrdelpartstot += nrdelpartsloc.at(kthr);
        nrdelpartsloc.at(kthr) = 0;
    }
    return;
}

pair<double,double> Locs::getwminmax_all()
{
    double wminloc = wminmaxloc.at(0).first;
    double wmaxloc = wminmaxloc.at(0).second;

    for(int i = 1; i < wminmaxloc.size(); i++)
    {
        wminloc = min(wminloc, wminmaxloc.at(i).first);
        wmaxloc = max(wmaxloc, wminmaxloc.at(i).second);
    }

    return pair<double,double>(wminloc,wmaxloc);
}


void ParaData::update_time(int stepscur)
{
    pasteps = stepscur;
    t = ((double)(stepscur+1))*dt;
}

bool ParaData::console_thisround(PICParams pp)
{
    return ((pasteps < 5) || rhomtrajout_thisround(pp));
}

bool ParaData::globout_thisround(PICParams pp)
{
    int frach = tsteps/(pp.am_tplots*100);
    if(frach <= 0)
        frach = 1;
    return ( (pasteps+1)%frach == 0 || pasteps == tsteps-1);
}

bool ParaData::rhomtrajout_thisround(PICParams pp)
{
    int frach = tsteps/pp.am_tplots;
    if(frach <= 0)
        frach = 1;
    return ( (pasteps+1)%frach == 0 || pasteps == tsteps-1);
}

bool ParaData::writeFiles(PICParams pp)
{
    int frach = (20*tsteps)/pp.am_tplots;
    if(frach <= 0)
        frach = 1;
    return ( (pasteps+1)%frach == 0 || pasteps == tsteps-1);
}

bool ParaData::writeParticleFile(PICParams pp)
{
    if(!pp.partsout)
        return false;
    return writeFiles(pp); //Always do both

    int frach = (100*tsteps)/pp.am_tplots;
    if(frach <= 0)
        frach = 1;
    return ( (pasteps+1)%frach == 0 || pasteps == tsteps-1);
}

double ParaData::get_mtoterrorrel()
{
    return ((mtot-mtoti)/mtoti);
}

void ParaData::startclock()
{
    starttime = chrono::steady_clock::now();
}

void ParaData::setloutime()
{
    loutime = chrono::steady_clock::now();
}

double ParaData::get_perTloop()
{
    return 100.*((double)(pasteps+1))/((double)(tsteps));
}

string ParaData::get_tslo()
{
    chrono::steady_clock::time_point curtime = chrono::steady_clock::now();
    double tslosec = (std::chrono::duration_cast<std::chrono::milliseconds>(curtime-loutime).count())*.001;
    return time_output(tslosec);
}

string ParaData::get_etl()
{
    chrono::steady_clock::time_point curtime = chrono::steady_clock::now();
    double etlsec = (std::chrono::duration_cast<std::chrono::milliseconds>(curtime-starttime).count())*.001*(((double)(tsteps))/((double)(pasteps+1)) - 1.);
    return time_output(etlsec);
}

void ParaData::setrminparts_seq(vector<ParPos> *parpos)
{
    double rminloc = parpos->at(0).r;
    for(int i = 1; i < parpos->size(); i++)
        rminloc = min(rminloc,parpos->at(i).r);
    rminparts = rminloc;
}

void ParaData::setrmaxparts_seq(vector<ParPos> *parpos)
{
    double rmaxloc = parpos->at(0).r;
    for(int i = 1; i < parpos->size(); i++)
        rmaxloc = max(rmaxloc,parpos->at(i).r);
    rmaxparts = rmaxloc;
}

void ParaData::setwminparts_seq(vector<ParPos> *parpos)
{
    double wminloc = parpos->at(0).r;
    for(int i = 1; i < parpos->size(); i++)
        wminloc = min(wminloc,parpos->at(i).w);
    wminparts = wminloc;
}

void ParaData::setwmaxparts_seq(vector<ParPos> *parpos)
{
    double wmaxloc = parpos->at(0).r;
    for(int i = 1; i < parpos->size(); i++)
        wmaxloc = max(wmaxloc,parpos->at(i).w);
    wmaxparts = wmaxloc;
}

void addSingleParticle_NOVOL_moms(const ParPos *parpos, RadFct *rho, int wmom)
{
    double rloc = parpos->r;
    double wloc = parpos->w;
    double Lloc = parpos->L;
    double fvloc = parpos->f*parpos->vol;

    double dr = rho->getdr();

    double rhelper = rloc-rho->getrmin();
    int indr = ceil(rhelper/dr);
    if(indr == 0)
    {
        cout << "r=" << rloc << "\n";
    }
    double wtl = ((double)indr) - rhelper/dr;
    double wtr = 1. - wtl;

    double wfactor = 1.;
    if(wmom != 0)
    {
        if(wmom < 0)
            exit(665); //Negative Moments are bad due to integrability...

        wfactor = wloc;
        for(int j = 1; j < wmom; j++)
            wfactor *= wloc;
    }

    rho->addind(indr, (wtr*wfactor*(fvloc/rloc))/rloc);
    rho->addind(indr-1, (wtl*wfactor*(fvloc/rloc))/rloc);
}

void addSingleParticle_num(const ParPos *parpos, double *res, int vmom)
{
    double rloc = parpos->r;
    double wloc = parpos->w;
    double Lloc = parpos->L;
    double fvloc = parpos->f*parpos->vol;

    double vfactor = 1.;
    if(vmom != 0)
    {
        if(vmom < 0)
            exit(563101); //Negative Moments are bad due to integrability...

        double vabs = wloc*wloc;
        if(rloc > numeric_limits<double>::epsilon() && Lloc > numeric_limits<double>::epsilon())
            vabs += Lloc/(rloc*rloc);
        vabs = sqrt(vabs);
        vfactor = vabs;
        for(int j = 1; j < vmom; j++)
            vfactor *= vabs;
    }

    *res += vfactor*fvloc;
}

void addSingleParticle_fac(const ParPos *parpos, double *res, RadFct *radfac)
{
    double rloc = parpos->r;
    double fvloc = parpos->f*parpos->vol;

    *res += fvloc*radfac->atr(rloc);
}

/*
void comprho_parahelper(const vector<ParPos> *parpos, vector<RadFct> *rholoc, const vector<int> *partpartition, int tid)
{
    rholoc->at(tid).setallzero();
    //Compute rho using the batch of particles assigned to this thread
    for(int kp = partpartition->at(tid); kp < partpartition->at(tid+1); kp++)
        addSingleParticle_NOVOL_moms(&parpos->at(kp), &rholoc->at(tid), 0);
}


void merge_locrhos(vector<RadFct> *rholoc, RadFct *rho, const vector<int> *rpartition, int tid, int num_thr)
{
    for(int kr = rpartition->at(tid); kr < rpartition->at(tid+1); kr++) //Sets rho in radial part of the thread to zero
        rho->setind(kr, 0.);

    double overdvolloc = 4.*M_PI*rho->getdr();

    for(int kr = rpartition->at(tid); kr < rpartition->at(tid+1); kr++)
    {
        for(int kthr = 0; kthr < num_thr; kthr++)
            rho->addind(kr, rholoc->at(kthr).atind(kr));

        rho->divind(kr, overdvolloc); //divide radial volume element, as it is already contained in whole volume element
    }
}

void compEkin_parahelper(const vector<ParPos> *parpos, vector<double> *ekinloc, const vector<int> *partpartition, int tid)
{
    ekinloc->at(tid) = 0.;
    //Compute Ekin using the batch of particles assigned to this thread
    for(int kp = partpartition->at(tid); kp < partpartition->at(tid+1); kp++)
        addSingleParticle_num(&parpos->at(kp), &ekinloc->at(tid), 2);
    ekinloc->at(tid) *= .5;
}
*/

PICParams::PICParams()
{ //default numerical parameters; optimized for 3D steady state with radius R=1
    T = 200.;
    dtoverdr = 2.;
    am_tplots = 200;
    Nrplot = 200;
    am_trajparts = 1;
    with_chess = false;
    partsout = false;
    fname = "0";
}

void PICParams::readInParams(string parameterfilename) //reads in parameters from file; if no file is found, dont change anything
{
    ifstream parameterin(parameterfilename);
    double Tn, dtoverdrn;
    int am_tplotsn, Nrplotn, am_trajpartsn, with_chessn, partsoutn;
    string fnamen;

    bool success = true;

    if(! (parameterin >> Tn))
        success = false;
    if(! (parameterin >> dtoverdrn))
        success = false;
    if(! (parameterin >> am_tplotsn >> Nrplotn >> am_trajpartsn >> with_chessn))
        success = false;
    if(! (parameterin >> partsoutn))
        success = false;
    if(! (parameterin >> fnamen))
        success = false;
    parameterin.close();

    if(with_chessn != 1 && with_chessn != 0)
        success = false;
    if(partsoutn != 1 && partsoutn != 0)
        success = false;

    if(success)
    {
        T = Tn;
        am_tplots = am_tplotsn;
        dtoverdr = dtoverdrn;
        Nrplot = Nrplotn;
        am_trajparts = am_trajpartsn;
        fname = fnamen;
        if(with_chessn == 1)
            with_chess = true;
        else
            with_chess = false;
        if(partsoutn == 1)
            partsout = true;
        else
            partsout = false;

        cout << endl << "Using numerical parameters from file." << endl;
    }
    else
        cout << endl << "Can't read numerical parameters from file; using default parameters." << endl;
}


InitParams::InitParams()
{
    dr = .002;
    usteps = 100;
    asteps = 50;
    min_parts_mio = 5;
    fromFile = false;
    parallel = false;
}

void InitParams::readInParams(string parameterfilename)
{
    ifstream parameterin(parameterfilename);
    double drn;
    int ustepsn, astepsn, min_parts_mion, fromFilen;

    bool success = true;

    if(! (parameterin >> drn >> ustepsn >> astepsn))
        success = false;
    if(! (parameterin >> min_parts_mion))
        min_parts_mion = -1;
    if(! (parameterin >> fromFilen))
        fromFilen = 0;
    parameterin.close();

    if(success)
    {
        dr = drn;
        usteps = ustepsn;
        asteps = astepsn;
        min_parts_mio = min_parts_mion;
        if(fromFilen == 1)
            fromFile = true;
        else
            fromFile = false;
        parallel = false;

        cout << endl << "Using initial value parameters from file." << endl;
    }
    else
        cout << endl << "Can't read initial value parameters from file; using default parameters." << endl;
}

ChessPlotData::ChessPlotData()
{
    vector<Fields> fields();
    vector<double> ts();

    vector<vector<double>> rho1();
    vector<vector<double>> rho2();
    vector<vector<double>> m1();
    vector<vector<double>> Up2();
    vector<vector<double>> U1();
    vector<vector<double>> U2();
}

void ChessPlotData::newData(Fields *newfield, double t)
{
    Fields fieldloc = *newfield;

    fields.push_back(fieldloc);
    ts.push_back(t);

    vector<double> vechelper(ts.size(), 0.);

    rho1.push_back(vechelper);
    rho2.push_back(vechelper);
    m1.push_back(vechelper);
    Up2.push_back(vechelper);
    U1.push_back(vechelper);
    U2.push_back(vechelper);

    for(int i = 0; i < ts.size(); i++)
    {
        double helper = 4.*M_PI*fields.at(i).rho.lpnorm(&fields.at(ts.size()-1).rho, 1., 2.);
        rho1.at(ts.size()-1).at(i) = helper;
        if(i < ts.size()-1)
            rho1.at(i).push_back(helper);

        helper = 4.*M_PI*fields.at(i).rho.lpnorm(&fields.at(ts.size()-1).rho, 2., 2.);
        rho2.at(ts.size()-1).at(i) = helper;
        if(i < ts.size()-1)
            rho2.at(i).push_back(helper);

        helper = 4.*M_PI*fields.at(i).mass.lpnorm(&fields.at(ts.size()-1).mass, 1., 2.);
        m1.at(ts.size()-1).at(i) = helper;
        if(i < ts.size()-1)
            m1.at(i).push_back(helper);

        helper = 4.*M_PI*fields.at(i).potprime.lpnorm(&fields.at(ts.size()-1).potprime, 2., 2.);
        Up2.at(ts.size()-1).at(i) = helper;
        if(i < ts.size()-1)
            Up2.at(i).push_back(helper);

        helper = 4.*M_PI*fields.at(i).pot.lpnorm(&fields.at(ts.size()-1).pot, 1., 2.);
        U1.at(ts.size()-1).at(i) = helper;
        if(i < ts.size()-1)
            U1.at(i).push_back(helper);

        helper = 4.*M_PI*fields.at(i).pot.lpnorm(&fields.at(ts.size()-1).pot, 2., 2.);
        U2.at(ts.size()-1).at(i) = helper;
        if(i < ts.size()-1)
            U2.at(i).push_back(helper);
    }

    return;
}


stringstream ChessPlotData::outputstring()
{
    stringstream ss;

    ss << "# Data for L^p Differences of Several Field Quantities at different times\n";
    ss << "# Using same time step size as for rhomout & trajout\n";
    ss << "#\n#\n";
    ss << "# FORMAT: 1) t_1 |2) t_2 |3) rho L^1 |4) rho L^2 |5) m L^1 |6) U' L^2 |7) U L^1 |8) U L^2\n";
    ss << "#\n#\n";
    for(int i = 0; i < ts.size(); i++)
        for(int j = 0; j < ts.size(); j++)
            ss << ts.at(i) << " " << ts.at(j) << " " << rho1.at(i).at(j) << " " << rho2.at(i).at(j) << " "
               << m1.at(i).at(j) << " " << Up2.at(i).at(j) << " " << U1.at(i).at(j) << " " << U2.at(i).at(j) << "\n";

    return ss;
}


bool isPrime(int n)
{
    if (n <= 1)
        return false;

    for (int k = 2; k < n; k++)
        if (n % k == 0)
            return false;

    return true;
}
