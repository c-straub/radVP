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

#include "PIC_LVP.h"
#include "lvp_barrier.h"

void PIC_LVP_fromFiles()
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
    if(pictype == 1)
        with_response = true;
    else if(pictype == 2)
        with_response = false;
    else
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
    LinPert *pert = pertfromfile("in_pert.txt");
    if(pertfromfile_nr("in_pert.txt") == 3)
        ip.parallel = true; //Only use parallelised initialisation if we have perturbation type 3

    //PIC Parameter Input:
    PICParams pp;
    pp.readInParams("in_picParams.txt");

    //Call Function:
    PIC_LVP(ip, pert, &stst, pp, with_response);
}

vector<ParPos> init_LVP(InitParams params, LinPert *pert, SteadyState *stst, Outers *outs)
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

                //Because evaluating the perturbation may take very long, we parallelise this step if necessary
                vector<double> fvals(asteps, 0.);

                if(params.parallel)
                {
                    vector<thread> thr(asteps);
                    for(int kalpha = 0; kalpha < asteps; kalpha++)
                    {
                        double ag = ((double)kalpha + .5)*dag;
                        thr.at(kalpha) = thread(atrua_parahelper, rg, ug, ag, pert, stst, &fvals.at(kalpha));
                    }
                    for(int kalpha = 0; kalpha < asteps; kalpha++)
                        thr.at(kalpha).join();
                }

                //Now create new particle
                for(int kalpha = 0; kalpha < asteps; kalpha++)
                {
                    double ag = ((double)kalpha + .5)*dag;

                    ParPos parpos_loc;
                    parpos_loc.r = rg;
                    parpos_loc.w = ug*cos(ag);
                    parpos_loc.L = rg*rg*ug*ug*sin(ag)*sin(ag);
                    parpos_loc.E = energy;
                    if(params.parallel)
                        parpos_loc.f = fvals.at(kalpha);
                    else
                        parpos_loc.f = pert->atrwL(rg, ug*cos(ag), rg*rg*ug*ug*sin(ag)*sin(ag), stst);
                    parpos_loc.vol = 8.*M_PI*M_PI*dr*dugr*dag*rg*rg*ug*ug*sin(ag);
                    parpos_loc.deletepart = false;

                    if(fabs(stst->atrwl(parpos_loc.r, parpos_loc.w, parpos_loc.L)) > numeric_limits<double>::epsilon()) //makes sure that we are only on steady state support
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

            exit(563005);
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

void PIC_LVP(InitParams ip, LinPert *pert, SteadyState *stst, PICParams pp, bool with_response)
{
    OutersLVP outs(pp.fname);

    vector<ParPos> parpos;
    if(ip.fromFile)
        parpos = init_FromFile(pp.fname);
    else
        parpos = init_LVP(ip, pert, stst, &outs);

    pp.dr = ip.dr;

    PIC_LVP(&parpos, stst, pp, &outs, with_response);
}

void PIC_LVP(vector<ParPos> *parpos, SteadyState *stst, PICParams pp, OutersLVP *outs,  bool with_response)
{
    stringstream ssloc;
    ssloc << "\n\nStart PIC for LVP!\nSimulating ";
    if(with_response)
        ssloc << "(full) linearised ";
    else
        ssloc << "pure transport ";
    ssloc << "system.\n\n";
    outs->writeConsoleAndClean(&ssloc);

    outs->setrajinds(pp.am_trajparts, parpos->size()); //Set up indices of particles we track
    ssloc << "Tracking " << outs->trajinds.size() << " particles.\n";
    outs->writeConsoleAndClean(&ssloc);

    //Some abbreviations:
    double rmaxst = stst->getrmax();
    double rminst = stst->getrmin();
    double dr = pp.dr;
    int rstepsst = (int)floor((rmaxst-rminst)/dr);

    /////////////////////////////////////////////////
    /// Initialise fields at t=0 (non-parallelised)
    /////////////////////////////////////////////////

    double rmin = max(rminst - pp.buffer*dr, 0.);
    FieldsLVP fields(dr, rstepsst+2*pp.buffer, rmin);

    fields.comprho_seq(parpos);

    fields.compmass();
    double mtot = fields.getmtot();

    fields.comppot();
    double pot0 = fields.getpot0();

    //Compute energies initially (fully sequential)
    double ekinfree = fields.getEkinfree_seq(parpos, stst);
    double epotfree = fields.getEpot();
    double ekinlin = fields.getEkin_seq(parpos);
    double epotlin = fields.getEpotlin_seq(parpos, stst);

    /////////////////////////////////////////////////
    /// Set Up & Start Parallelised Part
    /////////////////////////////////////////////////

    LocsLVP locz(parpos->size(), &fields.rho);

    double dt = pp.dr*pp.dtoverdr;
    int tsteps = (int)ceil(pp.T/dt);
    ParaDataLVP pada(mtot, ekinfree, epotfree, ekinlin, epotlin, stst->getmtot(), pot0, dt, tsteps);

    //Set minimal and maximal radius and w of particles initially (fully sequential)
    pada.setrminparts_seq(parpos);
    pada.setrmaxparts_seq(parpos);
    pada.setwminparts_seq(parpos);
    pada.setwmaxparts_seq(parpos);


    //Write Streams & Files initially
    outs->initwrite_globout(pp);
    outs->initwrite_rhomout(pp, &fields, fields.getrmin(), fields.getrmax());
    outs->initwrite_trajout(pp, parpos->size());

    outs->write_globout(&pada);
    outs->write_rhomout(pp, pada.t, &fields, fields.getrmin(), fields.getrmax());
    outs->write_trajout(pada.t, parpos);

    outs->writeConsole(conTabularLVPHead());
    outs->data2console(&pada);
    outs->writeConsole(conTabularLVPSeparator());

    outs->chessout.newData(&fields, pada.t);

    outs->writeAllFiles();

    //Create an array of threads of size number_of_threads and start these threads
    vector<thread> thr(locz.number_of_threads);
    for(int kthr = 0; kthr < locz.number_of_threads; kthr++) //Start the threads
        thr.at(kthr) = thread(PIC_LVP_parahelper, kthr, pp, with_response,
                              stst,
                              parpos, &fields, &locz,
                              &pada,
                              outs);

    for(int kthr = 0; kthr < locz.number_of_threads; kthr++)
        thr.at(kthr).join();

    outs->writeConsole(string("\n\n\nDONE with linearised PIC :)\n\n"));

    return;
}

void PIC_LVP_parahelper(const int tid, PICParams pp, bool withresponse,
                        SteadyState *stst,
                        vector<ParPos> *parpos, FieldsLVP *fields, LocsLVP *locz,
                        ParaDataLVP *pada,
                        OutersLVP *outs)
{

    double dtstst = .5*pada->dt; //Use this fixed dt when solving the steady state characteristic system

    /////////////////////////////////////////////////
    /// The time loop (only one)
    /////////////////////////////////////////////////
    for(int kt = 0; kt < pada->tsteps; kt++)
    {
        if(tid == 0)
            pada->update_time(kt);

        if(withresponse)
        {
            wait_for_all_threads();

            /////////////////////////////////////////////////
            /// First RK4 Step -- Half step from orig pos a la explicit Euler; save "integrand" as k1
            /////////////////////////////////////////////////

            locz->rholoc.at(tid).setallzero();

            for(int kp = locz->partpartition.at(tid); kp < locz->partpartition.at(tid+1); kp++)
            {
                ParPos parposloc = parpos->at(kp);

                double rloc = parposloc.r;
                double wloc = parposloc.w;
                double Lloc = parposloc.L;
                double gloc = parposloc.f;

                double integral_helper = 0.;
                pair<double, double> rwfull = stst->trajectoriesolver(rloc, wloc, Lloc, dtstst, .5*pada->dt, &integral_helper, &fields->potprime);
                integral_helper *= stst->phiprimeatrwl(rloc, wloc, Lloc);

                locz->int_k1.at(kp) = 2.*integral_helper; //Factor 2 since integral contains "factor" .5*dt
                gloc += integral_helper;

                parposloc.r = rwfull.first;
                parposloc.w = rwfull.second;
                parposloc.f = gloc;

                addSingleParticle_NOVOL_moms(&parposloc, &locz->rholoc.at(tid), 0);
            }

            wait_for_all_threads();

            locz->mergerholocs_parahelper(tid, &fields->rho);

            wait_for_all_threads();

            if(tid == 0)
                fields->compmass(); //also sets potprime

            wait_for_all_threads();

            /////////////////////////////////////////////////
            /// Second RK4 Step -- Half step from orig pos a la implicit Euler using k1-data; save "integrand" as k2
            /////////////////////////////////////////////////

            locz->rholoc.at(tid).setallzero();

            for(int kp = locz->partpartition.at(tid); kp < locz->partpartition.at(tid+1); kp++)
            {
                ParPos parposloc = parpos->at(kp);

                double rloc = parposloc.r;
                double wloc = parposloc.w;
                double Lloc = parposloc.L;
                double gloc = parposloc.f;

                double integral_helper = 0.;
                pair<double, double> rwfull = stst->trajectoriesolver(rloc, wloc, Lloc, dtstst, .5*pada->dt, &integral_helper, &fields->potprime);
                integral_helper *= stst->phiprimeatrwl(rloc, wloc, Lloc);

                locz->int_k2.at(kp) = 2.*integral_helper; //Factor 2 since integral contains "factor" .5*dt
                gloc += integral_helper;

                parposloc.r = rwfull.first;
                parposloc.w = rwfull.second;
                parposloc.f = gloc;

                addSingleParticle_NOVOL_moms(&parposloc, &locz->rholoc.at(tid), 0);
            }

            wait_for_all_threads();

            locz->mergerholocs_parahelper(tid, &fields->rho);

            wait_for_all_threads();

            if(tid == 0)
                fields->compmass();

            wait_for_all_threads();

            /////////////////////////////////////////////////
            /// Third RK4 Step -- Full step from orig pos a la midpoint method using k2-data; save "integrand" as k3
            /////////////////////////////////////////////////

            locz->rholoc.at(tid).setallzero();

            for(int kp = locz->partpartition.at(tid); kp < locz->partpartition.at(tid+1); kp++)
            {
                ParPos parposloc = parpos->at(kp);

                double rloc = parposloc.r;
                double wloc = parposloc.w;
                double Lloc = parposloc.L;
                double gloc = parposloc.f;

                double integral_helper = 0.;
                pair<double, double> rwfull = stst->trajectoriesolver(rloc, wloc, Lloc, dtstst, pada->dt, &integral_helper, &fields->potprime);
                integral_helper *= stst->phiprimeatrwl(rloc, wloc, Lloc);

                locz->int_k3.at(kp) = integral_helper;
                gloc += integral_helper;

                parposloc.r = rwfull.first;
                parposloc.w = rwfull.second;
                parposloc.f = gloc;

                addSingleParticle_NOVOL_moms(&parposloc, &locz->rholoc.at(tid), 0);
            }

            wait_for_all_threads();

            locz->mergerholocs_parahelper(tid, &fields->rho);

            wait_for_all_threads();

            if(tid == 0)
                fields->compmass();

        } //end of if(withresponse)

        wait_for_all_threads();

        /////////////////////////////////////////////////
        /// Particle Propagation (Full-Step)
        /// In RK4-Case: Forth RK4 Step -- Full step from orig pos a la implicit Euler using k3-data,
        ///              and make final step using k1, k2, and k3
        /////////////////////////////////////////////////

        for(int kp = locz->partpartition.at(tid); kp < locz->partpartition.at(tid+1); kp++)
        {
            ParPos parposloc = parpos->at(kp);

            double rloc = parposloc.r;
            double wloc = parposloc.w;
            double Lloc = parposloc.L;
            double gloc = parposloc.f;

            //Full Step
            pair<double, double> rwfull;
            if(withresponse)
            {
                double integral_helper = 0.;
                rwfull = stst->trajectoriesolver(rloc, wloc, Lloc, dtstst, pada->dt, &integral_helper, &fields->potprime);
                integral_helper *= stst->phiprimeatrwl(rloc, wloc, Lloc);

                gloc += (locz->int_k1.at(kp) + 2.*locz->int_k2.at(kp)
                          + 2.*locz->int_k3.at(kp) + integral_helper)/6.; //Make final RK4 Step
            }
            else
                rwfull = stst->trajectoriesolver(rloc, wloc, Lloc, dtstst, pada->dt);

            parposloc.r = rwfull.first;
            parposloc.w = rwfull.second;
            parposloc.f = gloc;

            parpos->at(kp) = parposloc;
        } //end of particle propagation loop

        wait_for_all_threads();

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
            locz->compEkinfree_parahelper(tid, parpos, stst);
            locz->compEpotlin_parahelper(tid, parpos, stst);
            locz->rminmaxparahelper(tid, parpos);
        }

        wait_for_all_threads();

        locz->mergerholocs_parahelper(tid, &fields->rho); //Compute rho from local rhos; parallelised in r

        double mtotcur, pot0cur, ekinfreecur, epotfreecur, ekinlincur, epotlincur, rminparts, rmaxparts;

        if(tid == 0 && (globoutthisround || conoutthisround))
        {
            ekinlincur = locz->get_Ekin_all();
            ekinfreecur = locz->get_Ekinfree_all();
            epotlincur = locz->get_Epotlin_all();

            pair<double,double> rminmaxparts = locz->getrminmax_all();
            rminparts = rminmaxparts.first;
            rmaxparts = rminmaxparts.second;
        }

        wait_for_all_threads();

        if(tid == 0)
        {
            fields->compmass(); //Compute mass from rho, non-parallelised :/
            mtotcur = fields->getmtot();
            if(globoutthisround || conoutthisround)
            {
                fields->comppot(); //Compute Potential U from Mass
                pot0cur = fields->getpot0();
                epotfreecur = fields->getEpot(); //Computes potential energy from mass

                pada->update_vars(kt, mtotcur, pot0cur, ekinfreecur, epotfreecur, ekinlincur, epotlincur,
                                  rminparts, rmaxparts); //Update pada
            }
        }

        wait_for_all_threads();

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
                outs->write_rhomout(pp, pada->t, fields, fields->getrmin(), fields->getrmax());
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

        wait_for_all_threads();

        if(tid == 0 && pada->writeFiles(pp))
        {
            outs->writeAllFiles(); //Write stringstreams into file
            outs->writeConsole(conTabularLVPMidHead());

            //Read globout file and analyse periods
            string fnamein = outs->prefix_globout + "_" + outs->fname + ".txt";

            //Kinetic Free Energy
            periodSimple_FromFile(fnamein, 1, 5, 15, "freekin");
            dftNanalyse_FromFile(fnamein, 1, 5, 15, "freekin");
            dftNanalyse_FromFile(fnamein, 1, 5, 15, "freekin_2ndHalf", pada->t*.5);
            //Kinetic Lin Energy
            periodSimple_FromFile(fnamein, 1, 9, 15, "linekin");
            dftNanalyse_FromFile(fnamein, 1, 9, 15, "linekin");
            dftNanalyse_FromFile(fnamein, 1, 9, 15, "linekin_2ndHalf", pada->t*.5);
            //U(t,0)
            periodSimple_FromFile(fnamein, 1, 13, 15, "pot0");
            dftNanalyse_FromFile(fnamein, 1, 13, 15, "pot0");
            dftNanalyse_FromFile(fnamein, 1, 13, 15, "pot0_2ndHalf", pada->t*.5);
        }

        if(tid == 1 || (tid == 0 && locz->number_of_threads <= 1)) //Write particles to File if desired
        {
            if(pada->writeParticleFile(pp))
            {
                outs->write_partoutFile(parpos);
                outs->writeConsole(conTabularLVPSeparator());
                outs->writeConsole(conTabularLVPSeparator());
            }
        }


    }//End of Time-Loop (only one)


    return;
}


double FieldsLVP::getEkinfree_seq(const vector<ParPos> *parpos, SteadyState *stst)
{
    double returner = 0.;
    for(int kp = 0; kp < parpos->size(); kp++)
        addSingleParticle_num(&parpos->at(kp), &returner, 2, true, false, stst);

    return -.5*returner;
}

double FieldsLVP::getEpotlin_seq(const vector<ParPos> *parpos, SteadyState *stst)
{
    double returner = 0.;
    for(int kp = 0; kp < parpos->size(); kp++)
        addSingleParticle_num(&parpos->at(kp), &returner, 1, false, true, stst);

    return returner;
}

LocsLVP::LocsLVP(int am_parts, RadFct *rho)
{
    number_of_threads = number_of_threads_const;
    setpartpartition(am_parts);
    setrpartition(rho->getnr());

    rholoc.clear();
    for(int i = 0; i < number_of_threads; i++)
        rholoc.push_back(RadFct(rho->getdr(), rho->getnr(), rho->getrmin()));

    ekinloc.clear();
    ekinloc.resize(number_of_threads, 0.);

    epotloc.clear();
    epotloc.resize(number_of_threads, 0.);

    ekinfreeloc.clear();
    ekinfreeloc.resize(number_of_threads, 0.);

    epotlinloc.clear();
    epotlinloc.resize(number_of_threads, 0.);

    int_k1 = vector<double>(am_parts, 0.);
    int_k2 = vector<double>(am_parts, 0.);
    int_k3 = vector<double>(am_parts, 0.);

    rminmaxloc = vector<pair<double,double>>(number_of_threads);
    wminmaxloc = vector<pair<double,double>>(number_of_threads);

    nrdelpartsloc = vector<int>(number_of_threads, 0);
    nrdelpartstot = 0;
}

void LocsLVP::compEkinfree_parahelper(int tid, const vector<ParPos> *parpos, SteadyState *stst)
{
    ekinfreeloc.at(tid) = 0.;
    //Compute Ekin using the batch of particles assigned to this thread
    for(int kp = partpartition.at(tid); kp < partpartition.at(tid+1); kp++)
        addSingleParticle_num(&parpos->at(kp), &ekinfreeloc.at(tid), 2, true, false, stst);
    ekinfreeloc.at(tid) *= -.5;
}

double LocsLVP::get_Ekinfree_all()
{
    double returner = 0.;
    for(int i = 0; i < ekinfreeloc.size(); i++)
        returner += ekinfreeloc.at(i);
    return returner;
}

void LocsLVP::compEpotlin_parahelper(int tid, const vector<ParPos> *parpos, SteadyState *stst)
{
    epotlinloc.at(tid) = 0.;
    //Compute Ekin using the batch of particles assigned to this thread
    for(int kp = partpartition.at(tid); kp < partpartition.at(tid+1); kp++)
        addSingleParticle_num(&parpos->at(kp), &epotlinloc.at(tid), 1, false, true, stst);
}

double LocsLVP::get_Epotlin_all()
{
    double returner = 0.;
    for(int i = 0; i < epotlinloc.size(); i++)
        returner += epotlinloc.at(i);
    return returner;
}

ParaDataLVP::ParaDataLVP(double mtotin, double ekinfreein, double epotfreein, double ekinlinin, double epotlinin,
                         double mststin, double pot0in,
                         double dti, int tstepsi)
{
    t = 0.;
    pasteps = 0;

    mtot = mtotin;
    mtoti = mtotin;

    pot0 = pot0in;

    efree = ekinfreein + epotfreein;
    efreei = efree;
    ekinfree = ekinfreein;
    epotfree = epotfreein;

    elin = ekinlinin + epotlinin;
    elini = elin;
    ekinlin = ekinlinin;
    epotlin = epotlinin;

    mstst = mststin;

    dt = dti;
    tsteps = tstepsi;

    starttime = chrono::steady_clock::now();
    loutime = chrono::steady_clock::now();
}

double ParaDataLVP::get_efreeerrorrel()
{
    return ((efree-efreei)/efreei);
}

double ParaDataLVP::get_elinerrorrel()
{
    return ((elin-elini)/elini);
}

void ParaDataLVP::update_vars(int stepscur, double mtotcur, double pot0cur, double ekinfreecur, double epotfreecur, double ekinlincur, double epotlincur,
                              double rmincur, double rmaxcur)
{
    update_time(stepscur);

    mtot = mtotcur;
    pot0 = pot0cur;

    efree = ekinfreecur + epotfreecur;
    ekinfree = ekinfreecur;
    epotfree = epotfreecur;

    elin = ekinlincur + epotlincur;
    ekinlin = ekinlincur;
    epotlin = epotlincur;

    rminparts = rmincur;
    rmaxparts = rmaxcur;
}

void addSingleParticle_num(const ParPos *parpos, double *res, int fpow, bool overphiprime, bool factorU0, SteadyState *stst)
{
    double rloc = parpos->r;
    double wloc = parpos->w;
    double Lloc = parpos->L;

    double floc = 1.;
    if(fpow < 0)
        exit(563102); //Negative powers are bad due to integrability
    if(fpow > 0)
    {
        for(int i = 0; i < fpow; i++)
            floc *= parpos->f;
    }
    double fvloc = floc*parpos->vol;

    if(overphiprime)
    {
        double phiprimedenom = stst->phiprimeatrwl(rloc, wloc, Lloc);
        if(fabs(phiprimedenom) > numeric_limits<double>::epsilon())
            fvloc /= phiprimedenom;
        else
            fvloc = 0.;
    }

    if(factorU0)
        fvloc *= stst->pot_atr(rloc);

    *res += fvloc;
}

OutersLVP::OutersLVP(string fnamei)
{
    fname = fnamei;
    reset();
}

void OutersLVP::initwrite_globout(PICParams pp)
{
    globout << "#\n#\n";
    globout << "# FORMAT: 1) t |2) t/M_0 |3) M(t) |4) (M(t)-M(0))/M(0)"
            << "|5) Free Ekin |6) Free Epot |7) Free Etot |8) (FreEtot(t)-FreEtot(0))/FreEtot(0)"
            << "|9) Lin. Ekin |10) Lin. Epot |11) Lin. Etot |12) (LinEtot(t)-LinEtot(0))/LinEtot(0)"
            << "|13) U(t,0) |14) Rmin |15) Rmax\n";
    globout << "#\n#\n";
}

void OutersLVP::write_globout(ParaDataLVP *pada)
{
    globout << setprecision(9) << pada->t << " " << pada->t/pada->mstst << " " << pada->mtot << " " << pada->get_mtoterrorrel() << " "
            << pada->ekinfree << " " << pada->epotfree << " " << pada->efree << " " << pada->get_efreeerrorrel() << " "
            << pada->ekinlin << " " << pada->epotlin << " " << pada->elin << " " << pada->get_elinerrorrel() << " "
            << pada->pot0 << " " << pada->rminparts << " " << pada->rmaxparts << setprecision(6) << "\n";
}

void OutersLVP::data2console(ParaDataLVP *pada)
{
    vector<string> cols;
    cols.push_back(to_string_precision(pada->t,4));
    cols.push_back(to_string_precision(pada->get_perTloop(),3));
    cols.push_back(pada->get_tslo());
    cols.push_back(pada->get_etl());
    cols.push_back(to_string(pada->mtot));
    cols.push_back(to_string(pada->efree));
    cols.push_back(to_string(pada->ekinfree));
    cols.push_back(to_string(pada->pot0));
    cols.push_back(to_string(pada->rmaxparts));

    writeConsole(conTabularRow(cols));
}

string conTabularLVPHead()
{
    vector<string> colhe;
    colhe.push_back(string("t"));
    colhe.push_back(string("% t-loop"));
    colhe.push_back(string("TSLO"));
    colhe.push_back(string("ETL"));
    colhe.push_back(string("Total Mass"));
    colhe.push_back(string("Free En."));
    colhe.push_back(string("Kin. FreEn."));
    colhe.push_back(string("U_f(t,0)"));
    colhe.push_back(string("Rmax"));

    return (conTabularRow(colhe) + conTabularSeparatorRow(colhe.size()));
}

string conTabularLVPSeparator()
{
    return conTabularSeparatorRow(9);
}

string conTabularLVPMidHead()
{
    return (conTabularLVPSeparator() + conTabularLVPHead());
}


double LPwpowerphi::atrwL(double r, double w, double L, SteadyState *stst)
{
    double wpart = 0.;
    double C = 1.;
    if(j < 0)
        C = -1.;
    else if(j == 0)
        wpart = 1.;
    else
    {
        wpart = w;
        for(int i = 1; i < std::abs(j); i++)
            wpart *= w;
    }

    return C*wpart*stst->atrwl(r,w,L);
}

double LPwpowerphiprime::atrwL(double r, double w, double L, SteadyState *stst)
{
    double wpart = 0.;
    double C = 1.;
    if(j < 0)
        C = -1.;
    else if(j == 0)
        wpart = 1.;
    else
    {
        wpart = w;
        for(int i = 1; i < std::abs(j); i++)
            wpart *= w;
    }

    return C*wpart*stst->phiprimeatrwl(r,w,L);
}

double LPsqrtrwphiprime::atrwL(double r, double w, double L, SteadyState *stst)
{
    double wpart = 0.;
    double rpart = 0.;
    double C = 1.;
    if(j < 0)
        C = -1.;
    else if(j == 0)
    {
        rpart = 1.;
        wpart = 1.;
    }
    else
    {
        rpart = r;
        wpart = w;
        for(int i = 1; i < std::abs(j); i++)
        {
            rpart *= r;
            wpart *= w;
        }
    }

    return C*sqrt(fabs(rpart+wpart))*stst->phiprimeatrwl(r,w,L);
}

double smoothcutoffer(double x)
{
    if(x < std::numeric_limits<double>::epsilon())
        return 0.;
    else
        return exp(-1./x);
}

double LPsqrtrwphiprimeawaywell::atrwL(double r, double w, double L, SteadyState *stst)
{
    double wpart = 0.;
    double rpart = 0.;
    double C = 1.;
    if(j < 0)
        C = -1.;
    else if(j == 0)
    {
        rpart = 1.;
        wpart = 1.;
    }
    else
    {
        rpart = r;
        wpart = w;
        for(int i = 1; i < std::abs(j); i++)
        {
            rpart *= r;
            wpart *= w;
        }
    }

    double cutofferarg = (stst->energy(r,w,L) - (stst->Emin(L) + convanish*(stst->gete0() - stst->potati(0)))) / (stst->gete0() - stst->potati(0));

    return C*sqrt(fabs(rpart+wpart))*stst->phiprimeatrwl(r,w,L)*smoothcutoffer(cutofferarg);
}

LinPert* pertfromfile(string parameterfilename)
{
    LinPert* returner;
    ifstream parameterin(parameterfilename);

    bool success = true;

    int perttype;

    if(! (parameterin >> perttype))
        success = false;

    if(perttype == 0)
    {
        int j = 0;
        if(! (parameterin >> j))
            success = false;
        returner = new LPwpowerphi(j);
    }
    else if(perttype == 1)
    {
        int j = 0;
        if(! (parameterin >> j))
            success = false;
        returner = new LPwpowerphiprime(j);
    }
    else if(perttype == 2)
    {
        int j = 0;
        if(! (parameterin >> j))
            success = false;
        returner = new LPsqrtrwphiprime(j);
    }
    else if(perttype == 3)
    {
        int j = 0; double convanish = .5;
        if(! (parameterin >> j >> convanish))
            success = false;
        returner = new LPsqrtrwphiprimeawaywell(j, convanish);
    }
    else
    {
        cout << endl << "Illegal Perturbation Type for LVP Perturbation." << endl << endl;
        exit(563001);
    }

    if(!success)
    {
        cout << endl << "Error when reading LVP Perturbation File." << endl << endl;
        exit(563002);
    }

    parameterin.close();

    return returner;
}

int pertfromfile_nr(string parameterfilename)
{
    ifstream parameterin(parameterfilename);

    bool success = true;

    int perttype;

    if(! (parameterin >> perttype))
        success = false;

    parameterin.close();

    if(success)
        return perttype;
    return -1;
}

void atrua_parahelper(double r, double u, double a, LinPert *pert, SteadyState *stst, double *resval)
{
    double wloc = u*cos(a);
    double Lloc = r*r*u*u*sin(a)*sin(a);

    *resval = pert->atrwL(r, wloc, Lloc, stst);
    return;
}
