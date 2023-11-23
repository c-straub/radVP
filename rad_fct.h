#ifndef RAD_FCT_H_INCLUDED
#define RAD_FCT_H_INCLUDED

#include <vector>
#include <iostream>

using namespace std;

class RadFct {  //class for radial functions, e.g., rho, U, mass, etc.; defined on a radial grid of uniform step length starting at some radius
   protected:
        double dr; //step length
        int nr; //amount of steps
        vector<double> val; //given values at k*dr for 0 <= k < nr
        double rmin;
    public:
        double getdr() {return dr;}
        int getnr() {return nr;}
        vector<double> getvalues() {return val;}
        int getvaluessize() {return val.size();}
        double getrmin() {return rmin;}

        RadFct() //default constructor
        {
            dr = 0.;
            nr = 0;
            val = vector<double>();
            rmin = 0.;
        }

        RadFct(double drinit, int nrinit, double rmini)
        {
            dr = drinit;
            nr = nrinit;
            val = vector<double>(nr, 0.);
            rmin = rmini;
        }

        RadFct(double drinit, int nrinit, double rmini, vector<double> valinit)
        {
            dr = drinit;
            nr = nrinit;
            rmin = rmini;
            val = valinit;
        }

        void setdr(double drnew) {dr = drnew; return;}
        void setnr(int nrnew) {nr = nrnew; return;}
        void setvals(vector<double> newvals) {val = newvals; return;}

        double atind(int ind)
        {
            try {return val.at(ind);}
            catch(...)
            {
                cout << endl << endl << "Error when accessing radial function values..." << endl << "Index is " << ind << ", array length " << val.size() << "..." << endl;
                exit(4);
            }
        }

        void addind(int ind, double x)
        {
            try {val.at(ind) += x;}
            catch(...)
            {
                cout << endl << endl << "Error when accessing radial function values for summation..." << endl << "Index is " << ind << ", array length " << val.size() << "..." << endl;
                exit(6);
            }
        }

        void multind(int ind, double x)
        {
            try {val.at(ind) *= x;}
            catch(...)
            {
                cout << endl << endl << "Error when accessing radial function values for multiplication..." << endl << "Index is " << ind << ", array length " << val.size() << "..." << endl;
                exit(6);
            }
        }

        void divall(double x)
        {
            for(int i = 0; i < val.size(); i++)
                val.at(i) /= x;
        }

        void divind(int ind, double x)
        {
            try {val.at(ind) /= x;}
            catch(...)
            {
                cout << endl << endl << "Error when accessing radial function values for division..." << endl << "Index is " << ind << " divisor is " << x << ", array length " << val.size() << "..." << endl;
                exit(6);
            }
        }

        void setallzero()
        {
            for(int i = 0; i < val.size(); i++)
                val.at(i) = 0.;
        }

        void setind(int ind, double x)
        {
            try {val.at(ind) = x;}
            catch(...)
            {
                cout << endl << endl << "Error when accessing radial function values for setting value..." << endl << "Index is " << ind << ", array length " << val.size() << "..." << endl;
                exit(6);
            }
        }

        void addtoprevind(int ind, double x)
        {
            try {val.at(ind) = val.at(ind-1)+x;}
            catch(...)
            {
                cout << endl << endl << "Error when accessing radial function values for summation to previous index..." << endl << "Indices are " << ind << " and " << ind-1 << ", array length " << val.size() << "..." << endl;
                exit(6);
            }
        }

        void addtonextind(int ind, double x)
        {
            try {val.at(ind) = val.at(ind+1)+x;}
            catch(...)
            {
                cout << endl << endl << "Error when accessing radial function values for summation to previous index..." << endl << "Indices are " << ind << " and " << ind-1 << ", array length " << val.size() << "..." << endl;
                exit(6);
            }
        }

        bool inrange(double rloc) //Return whether the radius is in the radial range of the radfct
        {
            double r = rloc - rmin;
            return (r>0 && r < (val.size()-1)*dr && r < (nr-1)*dr);
        }

        double get_last_entry()
        {
            double returner;

            try {returner = val.at(val.size()-1);}
            catch(...)
            {
                std::cout << std::endl << std::endl << "Error when accessing last radial function entry..." << std::endl << "Last index is " << val.size()-1 << ", array length " << val.size() << "..." << std::endl;
                exit(6);
            }

            return returner;
        }

        double getrmax()
        {
            return rmin + (nr-1)*dr;
        }

        double atr(double r); //evaluation at radius r
        double fctoverrsquareatr(double r); //evaluation of function over r^2 at radius r by linear interpolation

        double lpnorm(RadFct *other, double p, double rpower); //Computes l^p Norm between this function and other radial fct,
            //i.e., (int r^rpower |this(r) - other(r)|^p dr)^(1/p)
            //Use radial grid of this radial function
};


#endif // RADIAL_FCT_CLASS_H_INCLUDED
