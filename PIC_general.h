#ifndef PIC_GENERAL_H_INCLUDED
#define PIC_GENERAL_H_INCLUDED

/////////////////////////////////////////////////
/// File for PIC Stuff used for both linearised and non-linearised part
/////////////////////////////////////////////////

#include <string>
#include <vector>
#include <sstream>
#include <chrono> //For time measurement
#include <fstream>
#include <iostream>

#include "rad_fct.h"
#include "per_analyzer.h"

using namespace std;


class ParPos { //Particle Position in (r,w,L), including values of energy, f, and cell size
    public:
        double r; //Radius
        double w; //Radial Momentum
        double L; //Squared Modolus of Angular Momentum
        double E; //Particle Energy
        double f; //f(r,w,L)
        double vol; //Volume of cell

        bool deletepart; //Store whether this particle is to be "deleted"

        void serialise(ofstream& stream) const; //Bring ParPos Element into format for binary output file
        void deserialise(istream& stream); //Translate ParPos Element from binary file format back
};

class Fields //Wrapper Class for Potential, Mass, and Rho; Parent Class for more specialised children
{
    public:
        RadFct rho;
        RadFct mass;
        RadFct potprime;
        RadFct pot;

        Fields(double dri, int nri, double rmini)
        {
            rho = RadFct(dri, nri, rmini);
            mass = RadFct(dri, nri, rmini);
            potprime = RadFct(dri, nri, rmini);
            pot = RadFct(dri, nri, rmini);
        }

        double getdr() {return rho.getdr();}
        int getnr() {return rho.getnr();}
        double getrmax() {return rho.getrmax();}
        double getrmin() {return rho.getrmin();}

        void addSingleParticle_NOVOL(ParPos *parpos); //Adds influence of single particle to rho; but without adjusting volume element (keep values low for summation...)
        void comprho_seq(vector<ParPos> *parpos); //Computes rho from vector of particles; sequentially

        void compmass(); //Computes mass from rho (non-parallelised) and sets potprime
        double getmtot(); //Returns last entry of mass function; assumes that mass-array is already set

        double getEpot(); //Returns potential energy, i.e., -1/2 \int_0^\infty (m(r)/r)^2 dr = -1/(8Pi) \int_{\R^3} |\partial_x U(x)|^2 dx

        double getEkin_seq(vector<ParPos> *parpos); //Returns kinetic energy computed sequentially, i.e., 1/2 \int_{\R^3x\R^3} |v|^2 f(x,v) d(x,v)

        void comppot(); //Computes potential from mass (non-parallelised)
        double getpot0(); //Returns U(0)
};

class Locs { //Class for local quantities appearing in PIC codes, usually nedded throughout threads. Parent Class for more specialised children
    public:
        int number_of_threads;
        vector<int> partpartition; //Partition of numerical particles
        vector<int> rpartition; //Partition of radial range (save indices)
        vector<RadFct> rholoc; //Local Rho Vectors, used to compute rho parallelised
        vector<double> ekinloc; //Local ekin values, used to compute ekin parallelised
        vector<double> epotloc; //Local epot values, used to computed epot parallelised (using particles)
        vector<pair<double,double>> rminmaxloc; //Local Rmin and Rmax, used to compute Rmin & Rmax parallelised
        vector<pair<double,double>> wminmaxloc; //Local wmin and wmax, used to compute wmin & wmax parallelised

        vector<int> nrdelpartsloc; //Number of deleted particles per thread
        int nrdelpartstot;

        void setpartpartition(int am_parts);
        void setrpartition(int fieldslen);

        void comprho_parahelper(int tid, const vector<ParPos> *parpos); //Computes rho from particles WITHOUT VOLUME ELEMENT ADJUST
        void mergerholocs_parahelper(int tid, RadFct *rho); //Merges local rhos into global one  & adjust volume element; parallelised over r-range

        void compEkin_parahelper(int tid, const vector<ParPos> *parpos); //Computes ekin=ekinlin from particles
        double get_Ekin_all(); //Returns sum of all ekinlocs

        void compEpot_parahelper(int tid, const vector<ParPos> *parpos, Fields *fields); //Computes epot=epotfree from particles
        double get_Epot_all(); //Returns sum of all epotlocs

        void rminmaxparahelper(int tid, const vector<ParPos> *parpos); //Computes local RMin & Rmax
        pair<double,double> getrminmax_all(); //Returns pair (Rmin,Rmax), computed used rminmaxloc vector

        void wminmaxparahelper(int tid, const vector<ParPos> *parpos); //Computes local wMin & wmax
        pair<double,double> getwminmax_all(); //Returns pair (wmin,wmax), computed used wminmaxloc vector

        void deleteparts(int tid, vector<ParPos> *parpos); //"Deletes" particles if necessary, i.e., set f to 0 and position to the one of another particle
        void nrdelpartstot_update(); //Updates total number of deleted parts
};

void addSingleParticle_NOVOL_moms(const ParPos *parpos, RadFct *rho, int wmom); //Adds influence of single particle to rho, which is now a parameter, WITHOUT VOLUME ELEMENT ADJUST
    //with optional w-Moment, i.e., factor w^j
void addSingleParticle_num(const ParPos *parpos, double *res, int vmom); //Adds influence of single particle to number res
    //with optional v-Moment, i.e., factor |v|^j
void addSingleParticle_fac(const ParPos *parpos, double *res, RadFct *radfac); //Adds influence of single particle to number res
    //with radial function as factor, i.e., integrate f(x,v)*radfac(x)

class PICParams { //general parameters for particle-in-cell algorithm
    public:
        double T; //final time T
        double dtoverdr; //dt/dr, where dt is the time step size and dr radial step size (used for initialisation)
        int am_tplots; //amount of time plots
        int Nrplot; //number of dr steps taken for radial outputs
        int am_trajparts; //amount of particles tracked and saved in output file
        bool with_chess; //Determines whether we print the chess-plot
        bool partsout; //Determines whether we create output file with particle positions
        string fname; //Filename (suffix) for all output files

        int buffer = 200; //Buffer for fields; make them a bit larger to cope with particles flying a bit outside the original support
        double dr = 0.; //dr from initialisation; needed to set dt

        double getdt() {return dtoverdr*dr;}

        PICParams(); //default constructor, containing some default parameters for PIC

        void readInParams(string parameterfilename); //read in parameters from given file
};

class InitParams { //Parameters for initialisation
    public:
        double dr; //radial step size
        double usteps; //amount of usteps
        double asteps; //amount of asteps
        int min_parts_mio; //Minimal amount of overall particles for initialisation, in millions; reached by increasing usteps and asteps iteratively
        bool fromFile; //Determines whether we read in particles from a particle-file

        bool parallel; //Variable which decides whether initialisation is performed parallelised (slower for "trivial" perturbations!)

        InitParams(); //default constructor, containing some default parameters for PIC

        void readInParams(string parameterfilename); //read in parameters from given file
};

class ParaData { //Class for data used in parallelised part of PIC; parent for more specialised children
    public:
        double t; //Current time
        int pasteps; //Passed t-steps, not counting the current one (i.e., after first iteration this number is 0)

        double mtot; //Total Mass
        double mtoti; //Initial Total Mass

        double pot0; //Value of Potential at r=0=x

        double rminparts; //Minimal Radius accross all particles
        double rmaxparts; //Maximal Radius accross all particles

        double wminparts; //Minimal w-value of all particles
        double wmaxparts; //Maximal w-value of all particles

        double dt; //Time step size
        int tsteps; //Amount of t-steps needed to reach final time

        chrono::steady_clock::time_point starttime; //Start Time of Algorithm
        chrono::steady_clock::time_point loutime; //Time of last output

    public:
        void update_time(int stepscur);

        bool console_thisround(PICParams pp); //Returns whether we print data in console in this round
        bool globout_thisround(PICParams pp); //Returns whether we write globout in this round
        bool rhomtrajout_thisround(PICParams pp); //Returns whether we write rhomout & trajout in this round
        bool writeFiles(PICParams pp); //Returns whether we write sstreams to files in this round
        bool writeParticleFile(PICParams pp); //Returns whether we write particle file in this round

        double get_mtoterrorrel(); //Returns relative Mass Error

        void startclock(); //Sets starttime to current time
        void setloutime(); //Sets time of last out

        double get_perTloop(); //Returns % of done t-loop
        string get_tslo(); //Returns time since last output as string
        string get_etl(); //Returns etl (estimated time left) as string

        void setrminparts_seq(vector<ParPos> *parpos); //Sets rmin by scanning through particles sequentially
        void setrmaxparts_seq(vector<ParPos> *parpos); //Sets rmax by scanning through particles sequentially

        void setwminparts_seq(vector<ParPos> *parpos); //Sets wmin by scanning through particles sequentially
        void setwmaxparts_seq(vector<ParPos> *parpos); //Sets wmax by scanning through particles sequentially
};

class ChessPlotData { //Class for Data used to create "Chess Board Plot", i.e., ||U(t_1)-U(t_2)||_p etc
    public:
        vector<Fields> fields; //Saves field data used to create new plots
        vector<double> ts; //Saves t-values corresponding to above data

    private: //Save computed data
        vector<vector<double>> rho1; //L^1-norm of \rho differences
        vector<vector<double>> rho2; //L^2-norm of \rho differences
        vector<vector<double>> m1; //L^1-Norm of mass differences
        vector<vector<double>> Up2; //L^2-Norm of U' differences
        vector<vector<double>> U1; //L^1-Norm of U differences
        vector<vector<double>> U2; //L^2-Norm of U differences

    public:
        ChessPlotData(); //Default constructor

        void newData(Fields *newfield, double t); //Adds new data and adjust all of above attributes

        stringstream outputstring(); //Puts current data into string which can be used for plotting
};

class Outers //Wrapper Class for Several Output-Streams; Parent Class for more specialised children
{
    public:
        string fname; //Filename used for all files

        stringstream timecontroller; //Stringstream for same output as in console
        stringstream globout; //Stringstream for global data
        stringstream rhomout; //Stringstream for radial data
        stringstream trajout; //Stringstream for particle/trajectory data

        ChessPlotData chessout; //Data for ChessPlot

        string prefix_timecontroller = "out_tcontrol"; //Prefix for filename of timecontroller file
        string prefix_globout = "out_glob"; //Prefix for filename of globout file
        string prefix_rhomout = "out_rhom"; //Prefix for filename of rhomout file
        string prefix_trajout = "out_traj"; //Prefix for filename of trajout file
        string prefix_chessout = "out_chess"; //Prefix for filename of chessout file
        string prefix_partout = "outPART"; //Prefix for filename of particle file

        void reset(); //Clears all stringstreams and clears files
        void clearAllStreams(); //Clears all stringstreams

        vector<int> trajinds; //Vector of indices of particles we track
        void setrajinds(int am_trajparts, int numpart); //Sets indices of particles to track

        void writeConsole(string consoletext); //Write text in Console as well as timecontroller File
        void writeConsoleAndClean(stringstream *consoletext); //Same as above, but also clean stringstream

        void writeAllFiles(); //Write stringstreams to output files for globout, rhomout, trajout, chessout (timecontroller gets written in writeConsole)

        void initwrite_rhomout(PICParams pp, Fields *fields, double rpmin, double rpmax); //Writes initial information into rhomout stream
        void write_rhomout(PICParams pp, double t, Fields *fields, double rpmin, double rpmax); //Writes functions rho, m, and potential in rhomout stream

        void initwrite_trajout(PICParams pp, int am_part); //Writes initial information into trajout stream
        void write_trajout(double t, vector<ParPos> *parpos); //Writes positions of selected particles in stream

        void write_partoutFile(vector<ParPos> *parpos); //Writes positions of all particles into file
};

vector<ParPos> init_FromFile(string fnamesuff); //Reads in Particle Position Vector from binary input file

string conTabularRow(vector<string> row); //Returns row of tabular
string conTabularSeparatorRow(int columns); //Returns separator row of the neat little console tabular

string time_output(double secs); //Given certain number of seconds, return string of time in suitable unit (s, min, h, or d)

string to_string_precision(double x, int prec); //Converts given number x to string with given precision

bool isPrime(int n);  //Returns whether given integer is prime number; useful when choosing batch of particles

#endif // PIC_GENERAL_H_INCLUDED
