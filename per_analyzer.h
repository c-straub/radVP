#ifndef PER_ANALYZER_H_INCLUDED
#define PER_ANALYZER_H_INCLUDED

/////////////////////////////////////////////////
/// File for Period & Damping Analyzer
/////////////////////////////////////////////////

using namespace std;

vector<pair<double, double>> FileToVectorPair(string fname, int col1, int col2, int coltot,
                                              double shift2=0.); //Saves values of columns col1 and col2 from input file as vector<pair>
                                                                 //Second column gets shift given by shift2

double compPeriodSimple(vector<pair<double, double>> *dat, int returntype=7); //Return period of given data in simple way
    //Returntype determines precise way of computing period; good default values are 7 (period based on all local maxima) or 8 (period based on all zeros); see implementation for other values
double compPeriodSimple(vector<pair<double, double>> *dat, string fname, int returntype=7); //Same as above; writes detailed analysis in file

double periodSimple_FromFile(string fnamein, int col1, int col2, int coltot, string fnameout, double shift2=0., int returntype=7); //Executes above function using data from file, and writes data in file

vector<pair<double, double>> dft(vector<pair<double, double>> dat); //Computes (Non-Fast) Discrete Fourier Transformation; input are t-equidistant (t,x) pairs; output are complex Fourier Coefficients
void analyseDft(vector<pair<double, double>> dftdat, string fname, double dt); //analyzes dft coefficients by creating suitable plot files

vector<double> dftNanalyse(vector<pair<double, double>> dat, int initperiodgrid, string fname); //Analyzes data in vector using adaptive DFT scheme and put detailed mode analysis in file; returns periods found
vector<double> dftNanalyse(vector<pair<double, double>> dat, string fname); //use default value of 10 for initial grid size
vector<double> dftNanalyse(vector<pair<double, double>> dat, int initperiodgrid, string fname, double tstart); //Only consider data starting from given t
vector<double> dftNanalyse(vector<pair<double, double>> dat, string fname, double tstart);

vector<double> dftNanalyse_FromFile(string fnamein, int col1, int col2, int coltot,
                                    string fnameout, double tstart=0.,double shift2=0.); //Performes DFT analysis with data from file

bool detect_damping(vector<pair<double, double>> *dat, string fnameout, int gfactor=10);
bool detect_damping_FromFile(string fnamein, int col1, int col2, int coltot,
                             string fnameout, int gfactor=10);


double gplus(vector<pair<double, double>> *dat, double t); //Computes g_+(t)=max( )_{\infty,[t,\infty[} of given data
double gminus(vector<pair<double, double>> *dat, double t); //Computes g_-(t)=min( )_{\infty,[t,\infty[} of given data
double gabs(vector<pair<double, double>> *dat, double t); //Computes g_{abs}(t)=g_+(t)-g_-(t) of given data

pair<double,double> max_ininterval(vector<pair<double, double>> *dat, double left, double right); //Returns pair with maximal value of second value for first entries in given interval
pair<double,double> min_ininterval(vector<pair<double, double>> *dat, double left, double right); //Same as above with minimum
double maxfirst(vector<pair<double, double>> *dat); //Returns Maximum in First Variable
double minfirst(vector<pair<double, double>> *dat); //Returns Minimum in First Variable
double meansecond(vector<pair<double, double>> *dat, double tstart); //Returns mean in second variable for first variable >= tstart
double meansecond_FromFile(string fnamein, int col1, int col2, int coltot, double tstart);

#endif // PER_ANALYZER_H_INCLUDED
