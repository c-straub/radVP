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

#include "PIC_VP.h"
#include "vp_barrier.h"


void PIC_VP_fromFiles()
{
    bool success = true;

    //Basic Input:
    int pictype;
    ifstream parameterin("in_basic.txt");
    if(!(parameterin >> pictype))
    {
        cout << endl << "Can't read in_basic.txt file" << endl;
        return;
    }
    parameterin.close();

    bool with_response = true;
    if(pictype != 0)
    {
        cout << endl << "Unexpected input from in_basic.txt file..." << endl;
        return;
    }

    //Steady State Input:
    SteadyState stst = ststfromfile("in_equili.txt");

    //Initialisation Input:
    InitParams ip;
    ip.readInParams("in_initParams.txt");

    //Perturbation Input:
    double pertamp = 1.;
    ifstream pertin("in_pert.txt");
    if(!(pertin >> pertamp))
    {
        cout << endl << "Can't read in_pert.txt file" << endl;
        return;
    }
    pertin.close();

    //PIC Parameter Input:
    PICParams pp;
    pp.readInParams("in_picParams.txt");

    //Call Function:
    PIC_VP(ip, pertamp, &stst, pp);
}



vector<ParPos> init_VP(InitParams params, double pertamp, SteadyState *stst, Outers *outs)
{
    double rmax = stst->getrmax();
    double rmin = stst->getrmin();
    double e0 = stst->gete0();
    double dr = params.dr;

    int rstepsst = (int)((rmax-rmin)/dr);
    if(rmax > ((double)rstepsst)*dr + dr*.5 + rmin)
        rstepsst++;
    int usteps = params.usteps;
    int asteps = params.asteps;

    stringstream ssloc;
    ssloc << "\n" << "Start initialising particles." << "\n";
    ssloc << "We have " << rstepsst << " r-steps of dr=" << dr << " for the steady state to reach R=" << rmax << "; matter starts at r_min=" << rmin << ".\n";
    outs->writeConsoleAndClean(&ssloc);

    vector<ParPos> parpos; //array for numerical particles
    int oldsize = 0; //old size of parpos array; used to detect whether particle adding stagnates

    do
    { //Initialise particles until we have reached enough
        parpos.clear();

        double dug;
        if(e0 - stst->potati(0) > 0.)
            dug = sqrt(2.*e0 - 2.*stst->potati(0))/((double)usteps);
        else
        {
            cout << endl << endl << "Fatal error! Illegal and/or bad steady state..." << endl;
            exit(563003);
        }
        double dag = M_PI/((double)asteps);

        for(int kr = 0; kr < rstepsst; kr++)
        {
            double rg = ((double)kr + .5)*dr + rmin;

            double potr = stst->pot_atr(rg);

            double umaxr = 0.;

            if(e0 - potr > 0.)
                umaxr = sqrt(2.*e0 - 2.*potr);
            else
            {
                cout << endl << endl << "Fatal error! Illegal u-bound at r=" << rg << ", because " << e0 - potr << "<0 ..." << endl;
                exit(563004);
            }

            int ustepsr = ceil(umaxr/dug);

            double dugr = umaxr/((double)ustepsr);

            for(int ku = 0; ku < ustepsr; ku++)
            {
                double ug = ((double)ku + .5)*dugr;

                double energy = .5*ug*ug + potr; //E(rg, ug, alpha), independent of alpha

                //Now create new particle
                for(int kalpha = 0; kalpha < asteps; kalpha++)
                {
                    double ag = ((double)kalpha + .5)*dag;

                    ParPos parpos_loc;
                    parpos_loc.r = rg;
                    parpos_loc.w = ug*cos(ag);
                    parpos_loc.L = rg*rg*ug*ug*sin(ag)*sin(ag);
                    parpos_loc.E = energy;
                    parpos_loc.f = pertamp*stst->atrwl(rg, ug*cos(ag), rg*rg*ug*ug*sin(ag)*sin(ag));
                    parpos_loc.vol = 8.*M_PI*M_PI*dr*dugr*dag*rg*rg*ug*ug*sin(ag);
                    parpos_loc.deletepart = false;

                    if(fabs(parpos_loc.f) > numeric_limits<double>::epsilon()) //makes sure that we are only on support of initial condition
                        parpos.push_back(parpos_loc);
                }
            }
        } //end of r-init loop

        ssloc << "usteps=" << usteps << " and alphasteps=" << asteps << " lead to " << parpos.size() << " particles.\n";
        outs->writeConsoleAndClean(&ssloc);

        if(parpos.size() <= oldsize)
        {
            ssloc << "\n\n\nParticle adding process isnt making any progress...\n\n\n";
            outs->writeConsoleAndClean(&ssloc);

            exit(572801);
        }
        oldsize = parpos.size();

        //Increase asteps and usteps for case where we need more particles
        int astepsold = asteps;
        asteps = astepsold + 10;
        usteps += (usteps*asteps-usteps*astepsold)/astepsold;
    }
    while(parpos.size()+1 < params.min_parts_mio*1000000 && params.min_parts_mio > 0);

    return parpos;
}

void PIC_VP(InitParams ip, double pertamp, SteadyState* stst, PICParams pp)
{
    OutersVP outs(pp.fname);

    vector<ParPos> parpos;
    if(ip.fromFile)
        parpos = init_FromFile(pp.fname);
    else
        parpos = init_VP(ip, pertamp, stst, &outs);

    pp.dr = ip.dr;

    PIC_VP(&parpos, pp, &outs);
}

void PIC_VP(vector<ParPos> *parpos, PICParams pp, OutersVP *outs)
{
    stringstream ssloc;
    ssloc << "\n\nStart PIC for (non-linearised) VP!\n\n";
    outs->writeConsoleAndClean(&ssloc);

    outs->setrajinds(pp.am_trajparts, parpos->size()); //Set up indices of particles we track
    ssloc << "Tracking " << outs->trajinds.size() << " particles.\n";
    outs->writeConsoleAndClean(&ssloc);

    //Some abbreviations:
    double dr = pp.dr;

    /////////////////////////////////////////////////
    /// Initialise fields at t=0 (non-parallelised)
    /////////////////////////////////////////////////

    ParaDataVP padaloc; //Local instance of pada, used to compute rmin and rmax
    padaloc.setrminparts_seq(parpos);
    padaloc.setrmaxparts_seq(parpos);
    double rminpart = padaloc.rminparts;
    double rmaxpart = padaloc.rmaxparts;
    int rstepsinit = (int)floor((rminpart-rminpart)/dr);

    double rmin = max(rminpart - 100.*pp.buffer*dr, 0.);
    FieldsVP fields(dr, rstepsinit+200*pp.buffer, rmin);

    fields.comprho_seq(parpos);

    fields.compmass();
    double mtot = fields.getmtot();

    fields.comppot();
    double pot0 = fields.getpot0();

    //Compute energies initially (fully sequential)
    double ekin = fields.getEkin_seq(parpos);
    double epot = fields.getEpot();

    /////////////////////////////////////////////////
    /// Set Up & Start Parallelised Part
    /////////////////////////////////////////////////

    LocsVP locz(parpos->size(), &fields.rho);

    double dt = pp.dr*pp.dtoverdr;
    int tsteps = (int)ceil(pp.T/dt);
    double plotradialfurther = .1*(rmaxpart-rminpart);
    double rminplot = max({0., rminpart-plotradialfurther, fields.getrmin()});
    double rmaxplot = min(rmaxpart+plotradialfurther, fields.getrmax()-2.*fields.getdr());
    ParaDataVP pada(mtot, ekin, epot, pot0, dt, tsteps, rminplot, rmaxplot);

    //Set minimal and maximal radius and w of particles initially (fully sequential)
    pada.setrminparts_seq(parpos);
    pada.setrmaxparts_seq(parpos);
    pada.setwminparts_seq(parpos);
    pada.setwmaxparts_seq(parpos);

    //Write Streams & Files initially
    outs->initwrite_globout(pp);
    outs->initwrite_rhomout(pp, &fields, pada.rminplot, pada.rmaxplot);
    outs->initwrite_trajout(pp, parpos->size());

    outs->write_globout(&pada);
    outs->write_rhomout(pp, pada.t, &fields, pada.rminplot, pada.rmaxplot);
    outs->write_trajout(pada.t, parpos);

    outs->writeConsole(conTabularVPHead());
    outs->data2console(&pada);
    outs->writeConsole(conTabularVPSeparator());

    outs->chessout.newData(&fields, pada.t);

    outs->writeAllFiles();

    //Create an array of threads of size number_of_threads and start these threads
    vector<thread> thr(locz.number_of_threads);
    for(int kthr = 0; kthr < locz.number_of_threads; kthr++) //Start the threads
        thr.at(kthr) = thread(PIC_VP_parahelper, kthr, pp,
                              parpos, &fields, &locz,
                              &pada,
                              outs);

    for(int kthr = 0; kthr < locz.number_of_threads; kthr++)
        thr.at(kthr).join();

    outs->writeConsole(string("\n\n\nDONE with non-linearised PIC :)\n\n"));

    return;
}

void PIC_VP_parahelper(const int tid, PICParams pp,
                        vector<ParPos> *parpos, FieldsVP *fields, LocsVP *locz,
                        ParaDataVP *pada,
                        OutersVP *outs)
{
    /////////////////////////////////////////////////
    /// The time loop (only one)
    /////////////////////////////////////////////////
    for(int kt = 0; kt < pada->tsteps; kt++)
    {
        if(tid == 0)
            pada->update_time(kt);

        wait_for_all_threads_VP();


        /////////////////////////////////////////////////
        /// First RK4 Step -- Half step from orig pos a la explicit Euler; Save (x,v)^dot as k1
        /////////////////////////////////////////////////

        locz->rholoc.at(tid).setallzero();

        for(int kp = locz->partpartition.at(tid); kp < locz->partpartition.at(tid+1); kp++)
        {
            ParPos parposloc = parpos->at(kp);

            double rloc = parposloc.r;
            double wloc = parposloc.w;
            double Lloc = parposloc.L;

            XV xvdot = fields->RHScharsys(&parposloc);
            locz->par_k1.at(kp) = xvdot;

            XV xvnew = addWithFactor(XV(rloc, wloc, Lloc), .5*pada->dt, xvdot); //half Euler Step

            parposloc.r = xvnew.comp_r();
            parposloc.w = xvnew.comp_w(parposloc.r);

            addSingleParticle_NOVOL_moms(&parposloc, &locz->rholoc.at(tid), 0);
        }

        wait_for_all_threads_VP();

        locz->mergerholocs_parahelper(tid, &fields->rho);

        wait_for_all_threads_VP();

        if(tid == 0)
            fields->compmass(); //also sets potprime

        wait_for_all_threads_VP();

        /////////////////////////////////////////////////
        /// Second RK4 Step -- Half step from orig pos a la implicit Euler using fields & positions from first step
        /////////////////////////////////////////////////

        locz->rholoc.at(tid).setallzero();

        for(int kp = locz->partpartition.at(tid); kp < locz->partpartition.at(tid+1); kp++)
        {
            ParPos parposloc = parpos->at(kp);

            double rloc = parposloc.r;
            double wloc = parposloc.w;
            double Lloc = parposloc.L;

            //Repeat k1-step, so that we can evaluate RHS at position which fits the updated time
            XV xvk1 = addWithFactor(XV(rloc, wloc, Lloc), .5*pada->dt, locz->par_k1.at(kp));

            XV xvdot = fields->RHScharsys(&xvk1);
            locz->par_k2.at(kp) = xvdot;

            XV xvnew = addWithFactor(XV(rloc, wloc, Lloc), .5*pada->dt, xvdot); //half Euler Step

            parposloc.r = xvnew.comp_r();
            parposloc.w = xvnew.comp_w(parposloc.r);

            addSingleParticle_NOVOL_moms(&parposloc, &locz->rholoc.at(tid), 0);
        }

        wait_for_all_threads_VP();

        locz->mergerholocs_parahelper(tid, &fields->rho);

        wait_for_all_threads_VP();

        if(tid == 0)
            fields->compmass();

        wait_for_all_threads_VP();

        /////////////////////////////////////////////////
        /// Third RK4 Step -- Full step from orig pos a la midpoint using fields & positions from second step
        /////////////////////////////////////////////////

        locz->rholoc.at(tid).setallzero();

        for(int kp = locz->partpartition.at(tid); kp < locz->partpartition.at(tid+1); kp++)
        {
            ParPos parposloc = parpos->at(kp);

            double rloc = parposloc.r;
            double wloc = parposloc.w;
            double Lloc = parposloc.L;

            //Repeat k2-step, so that we can evaluate RHS at position which fits the updated time
            XV xvk2 = addWithFactor(XV(rloc, wloc, Lloc), .5*pada->dt, locz->par_k2.at(kp));

            XV xvdot = fields->RHScharsys(&xvk2);
            locz->par_k3.at(kp) = xvdot;

            XV xvnew = addWithFactor(XV(rloc, wloc, Lloc), pada->dt, xvdot); //full midpoint Step

            parposloc.r = xvnew.comp_r();
            parposloc.w = xvnew.comp_w(parposloc.r);

            addSingleParticle_NOVOL_moms(&parposloc, &locz->rholoc.at(tid), 0);
        }

        wait_for_all_threads_VP();

        locz->mergerholocs_parahelper(tid, &fields->rho);

        wait_for_all_threads_VP();

        if(tid == 0)
            fields->compmass();

        wait_for_all_threads_VP();

        /////////////////////////////////////////////////
        /// Fourth & Final RK4 Step -- Full step from orig pos a la implicit Euler using fields & positions from third step,
        ///              and make final step using k1, k2, and k3
        /////////////////////////////////////////////////

        locz->rholoc.at(tid).setallzero();

        for(int kp = locz->partpartition.at(tid); kp < locz->partpartition.at(tid+1); kp++)
        {
            ParPos parposloc = parpos->at(kp);

            double rloc = parposloc.r;
            double wloc = parposloc.w;
            double Lloc = parposloc.L;

            //Repeat k3-step, so that we can evaluate RHS at position which fits the updated time
            XV xvk3 = addWithFactor(XV(rloc, wloc, Lloc), pada->dt, locz->par_k3.at(kp));

            XV xvdot = fields->RHScharsys(&xvk3); //This is k4

            //Make final RK4-Step:
            XV xvnew = XV(rloc, wloc, Lloc);
            xvnew.rk4_update(locz->par_k1.at(kp), locz->par_k2.at(kp), locz->par_k3.at(kp), xvdot, pada->dt);

            ///XV xvnew = addWithFactor(XV(rloc, wloc, Lloc), pada->dt, fields->RHScharsys(&parposloc)); ///TEST: full Euler Step

            parposloc.r = xvnew.comp_r();
            parposloc.w = xvnew.comp_w(parposloc.r);

            //Check whether particle is to be deleted
            if(parposloc.r > .9*(fields->getrmax()-fields->getrmin()))
                parposloc.deletepart = true;
            else
                parposloc.deletepart = false;

            parpos->at(kp) = parposloc; //Save new position
        }

        wait_for_all_threads_VP();

        locz->deleteparts(tid, parpos); //Delete particles if necessary
        locz->nrdelpartstot_update();

        wait_for_all_threads_VP();

        /////////////////////////////////////////////////
        /// Compute rho, mass, U, and Energies
        /////////////////////////////////////////////////

        bool globoutthisround = pada->globout_thisround(pp);
        bool conoutthisround = pada->console_thisround(pp);
        bool rhomtrajoutthisround = pada->rhomtrajout_thisround(pp);

        locz->comprho_parahelper(tid, parpos); //Compute local rhos from particles

        if(globoutthisround || conoutthisround)
        {
            locz->compEkin_parahelper(tid, parpos);
            locz->rminmaxparahelper(tid, parpos);
            locz->wminmaxparahelper(tid, parpos);
        }

        wait_for_all_threads_VP();

        locz->mergerholocs_parahelper(tid, &fields->rho); //Compute rho from local rhos; parallelised in r

        double mtotcur, pot0cur, ekincur, epotcur, rminparts, rmaxparts, wminparts, wmaxparts;

        if(tid == 0 && (globoutthisround || conoutthisround))
        {
            ekincur = locz->get_Ekin_all();

            pair<double,double> rminmaxparts = locz->getrminmax_all();
            rminparts = rminmaxparts.first;
            rmaxparts = rminmaxparts.second;

            pair<double,double> wminmaxparts = locz->getwminmax_all();
            wminparts = wminmaxparts.first;
            wmaxparts = wminmaxparts.second;
        }

        wait_for_all_threads_VP();

        if(tid == 0)
        {
            fields->compmass(); //Compute mass from rho, non-parallelised :/
            mtotcur = fields->getmtot();
            if(globoutthisround || conoutthisround)
            {
                fields->comppot(); //Compute Potential U from Mass
                pot0cur = fields->getpot0();
                epotcur = fields->getEpot(); //Computes potential energy from mass

                pada->update_vars(kt, mtotcur, pot0cur, ekincur, epotcur,
                                  rminparts, rmaxparts, wminparts, wmaxparts, locz->nrdelpartstot); //Update pada
            }
        }

        wait_for_all_threads_VP();

        /////////////////////////////////////////////////
        /// Output
        /////////////////////////////////////////////////

        if(tid == 0)
        { //globout
            if(globoutthisround)
                outs->write_globout(pada);
        }
        if(tid == 1 || (tid == 0 && locz->number_of_threads <= 1))
        { //rhomout
            if(rhomtrajoutthisround)
                outs->write_rhomout(pp, pada->t, fields, pada->rminplot, pada->rmaxplot);
        }
        if(tid == 2 || (tid == 0 && locz->number_of_threads <= 2))
        { //trajout
            if(rhomtrajoutthisround)
                outs->write_trajout(pada->t, parpos);
        }
        if(tid == 3 || (tid == 0 && locz->number_of_threads <= 3))
        { //console
            if(conoutthisround)
            {
                outs->data2console(pada);
                pada->setloutime();
            }

        }
        if(tid == 4 || (tid == 0 && locz->number_of_threads <= 4))
        { //console
            if(rhomtrajoutthisround && pp.with_chess)
                outs->chessout.newData(fields, pada->t);
        }

        wait_for_all_threads_VP();

        if(tid == 0 && pada->writeFiles(pp))
        {
            outs->writeAllFiles(); //Write stringstreams into file
            outs->writeConsole(conTabularVPMidHead());

            //Read globout file and analyse periods
            string fnamein = outs->prefix_globout + "_" + outs->fname + ".txt";

            double ekinmean = meansecond_FromFile(fnamein, 1, 5, 13, pada->t*.5); //Mean over ekin
            //Kinetic Energy
            periodSimple_FromFile(fnamein, 1, 5, 13, "ekin", -ekinmean);
            dftNanalyse_FromFile(fnamein, 1, 5, 13, "ekin", 0., -ekinmean);
            dftNanalyse_FromFile(fnamein, 1, 5, 13, "ekin_2ndHalf", pada->t*.5, -ekinmean);
            //U(t,0)
            periodSimple_FromFile(fnamein, 1, 9, 13, "pot0");
            dftNanalyse_FromFile(fnamein, 1, 9, 13, "pot0");
            dftNanalyse_FromFile(fnamein, 1, 9, 13, "pot0_2ndHalf", pada->t*.5, 0.);
        }

        if(tid == 1 || (tid == 0 && locz->number_of_threads <= 1)) //Write particles to File if desired
        {
            if(pada->writeParticleFile(pp))
            {
                outs->write_partoutFile(parpos);
                outs->writeConsole(conTabularVPSeparator());
                outs->writeConsole(conTabularVPSeparator());
            }
        }
    }//End of Time-Loop (only one)


    return;

}


XV FieldsVP::RHScharsys(ParPos *papo)
{
    double rloc = papo->r;
    double wloc = papo->w;
    double Lloc = papo->L;

    return RHScharsys(rloc, wloc, Lloc);
}

XV FieldsVP::RHScharsys(XV *xv)
{
    double rloc = xv->comp_r();
    double Uprimeloc = potprime.atr(rloc);

    if(rloc < numeric_limits<double>::epsilon())
        return XV(xv->v1, xv->v2, xv->v3, 0., 0., 0.);
    else
        return XV(xv->v1, xv->v2, xv->v3, (-Uprimeloc*xv->x1)/rloc, (-Uprimeloc*xv->x2)/rloc, (-Uprimeloc*xv->x3)/rloc);
}


XV FieldsVP::RHScharsys(double rloc, double wloc, double Lloc)
{
    double Uprimeloc = potprime.atr(rloc);

    if(rloc < numeric_limits<double>::epsilon())
    {
        if(Lloc < numeric_limits<double>::epsilon())
            return XV(wloc, 0., 0., 0., 0., 0.);
        else
        {
            cout << endl << endl
                 << "Error; trying to compute RHS of char sys for " << rloc << "=r=0 < L=" << Lloc << "..."
                 << endl << endl;
            exit(572802);
        }
    }
    else
        return XV(wloc, sqrt(Lloc)/rloc, 0., -Uprimeloc, 0., 0.);
}

LocsVP::LocsVP(int am_parts, RadFct *rho)
{
    number_of_threads = number_of_threads_constVP;
    setpartpartition(am_parts);
    setrpartition(rho->getnr());

    rholoc.clear();
    for(int i = 0; i < number_of_threads; i++)
        rholoc.push_back(RadFct(rho->getdr(), rho->getnr(), rho->getrmin()));

    ekinloc.clear();
    ekinloc.resize(number_of_threads, 0.);

    epotloc.clear();
    epotloc.resize(number_of_threads, 0.);

    par_k1 = vector<XV>(am_parts);
    par_k2 = vector<XV>(am_parts);
    par_k3 = vector<XV>(am_parts);

    rminmaxloc = vector<pair<double,double>>(number_of_threads);
    wminmaxloc = vector<pair<double,double>>(number_of_threads);

    nrdelpartsloc = vector<int>(number_of_threads, 0);
    nrdelpartstot = 0;
}



ParaDataVP::ParaDataVP(double mtotin, double ekinin, double epotin, double pot0in,
                       double dti, int tstepsi,
                       double rminplotin, double rmaxplotin)
{
    t = 0.;
    pasteps = 0;

    mtot = mtotin;
    mtoti = mtotin;

    pot0 = pot0in;

    etot = ekinin + epotin;
    etoti = etot;
    ekin = ekinin;
    epot = epotin;
    ekini = ekinin;

    dt = dti;
    tsteps = tstepsi;

    starttime = chrono::steady_clock::now();
    loutime = chrono::steady_clock::now();

    rminplot = rminplotin;
    rmaxplot = rmaxplotin;

    nrdelparts = 0;
}

double ParaDataVP::get_etoterrorrel()
{
    return ((etot-etoti)/etoti);
}

void ParaDataVP::update_vars(int stepscur, double mtotcur, double pot0cur,
                             double ekincur, double epotcur,
                             double rmincur, double rmaxcur, double wmincur, double wmaxcur,
                             int nrdelpartscur)
{
    update_time(stepscur);

    mtot = mtotcur;
    pot0 = pot0cur;

    etot = ekincur + epotcur;
    ekin = ekincur;
    epot = epotcur;

    rminparts = rmincur;
    rmaxparts = rmaxcur;

    wminparts = wmincur;
    wmaxparts = wmaxcur;

    nrdelparts = nrdelpartscur;
}

OutersVP::OutersVP(string fnamei)
{
    fname = fnamei;
    reset();
}

void OutersVP::initwrite_globout(PICParams pp)
{
    globout << "#\n#\n";
    globout << "# FORMAT: 1) t |2) t/M(0) |3) M(t) |4) (M(t)-M(0))/M(0)"
            << "|5) Ekin |6) Epot |7) Etot |8) (Etot(t)-Etot(0))/Etot(0)"
            << "|9) U(t,0) |10) Rmin |11) Rmax |12) wmin |13) wmax\n";
    globout << "#\n#\n";
}

void OutersVP::write_globout(ParaDataVP *pada)
{
    globout << setprecision(9) << pada->t << " " << pada->t/pada->mtoti << " " << pada->mtot << " " << pada->get_mtoterrorrel() << " "
            << pada->ekin << " " << pada->epot << " " << pada->etot << " " << pada->get_etoterrorrel() << " "
            << pada->pot0 << " " << pada->rminparts << " " << pada->rmaxparts << " " << pada->wminparts << " " << pada->wmaxparts << setprecision(6) << "\n";
}

void OutersVP::data2console(ParaDataVP *pada)
{
    vector<string> cols;
    cols.push_back(to_string_precision(pada->t,4));
    cols.push_back(to_string_precision(pada->get_perTloop(),3));
    cols.push_back(pada->get_tslo());
    cols.push_back(pada->get_etl());
    cols.push_back(to_string(pada->mtot));
    cols.push_back(to_string(pada->etot));
    cols.push_back(to_string(pada->ekin));
    cols.push_back(to_string(pada->pot0));
    cols.push_back(to_string(pada->rmaxparts));
    cols.push_back(to_string(pada->nrdelparts));

    writeConsole(conTabularRow(cols));
}

string conTabularVPHead()
{
    vector<string> colhe;
    colhe.push_back(string("t"));
    colhe.push_back(string("% t-loop"));
    colhe.push_back(string("TSLO"));
    colhe.push_back(string("ETL"));
    colhe.push_back(string("Total Mass"));
    colhe.push_back(string("Tot. En."));
    colhe.push_back(string("Kin. En."));
    colhe.push_back(string("U_f(t,0)"));
    colhe.push_back(string("Rmax"));
    colhe.push_back(string("# Del.P."));

    return (conTabularRow(colhe) + conTabularSeparatorRow(colhe.size()));
}

string conTabularVPSeparator()
{
    return conTabularSeparatorRow(10);
}

string conTabularVPMidHead()
{
    return (conTabularVPSeparator() + conTabularVPHead());
}



