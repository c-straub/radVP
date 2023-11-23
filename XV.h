#ifndef XV_H_INCLUDED
#define XV_H_INCLUDED

/////////////////////////////////////////////////
/// (Helper) class for (x,v)\in\R^3\times\R^3; useful for solving characteristic system
/////////////////////////////////////////////////

#include <cmath>

using namespace std;

class XV //helper class for (x,v)\in\R^3\times\R^3; useful for solving characteristic system
{
    public:
        double x1;
        double x2;
        double x3;
        double v1;
        double v2;
        double v3;

        XV()
        {
            x1=0.; x2=0.; x3=0.;
            v1=0.; v2=0.; v3=0.;
        }

        XV(double x1i, double x2i, double x3i, double v1i, double v2i, double v3i)
        {
            x1 = x1i; x2 = x2i; x3 = x3i;
            v1 = v1i; v2 = v2i; v3 = v3i;
        }

        XV(double ri, double wi, double Li); //Creates (non-unique) XV vector for given (r,w,L)

        double comp_r() //Computes r=|x|
        {
            return sqrt(x1*x1 + x2*x2 + x3*x3);
        }

        double comp_w(); //Computes w=(x/r) \cdot v
        double comp_w(double r); //Computes w, where r is already computed

        void mult_factor(double factor) //Multiply whole (x,v) vector with scalar factor
        {
            x1 *= factor; x2 *= factor; x3 *= factor;
            v1 *= factor; v2 *= factor; v3 *= factor;
        }

        void rk4_update(XV k1, XV k2, XV k3, XV k4, double dt); //Updates (x,v) vector according to RK4
};

XV addWithFactor(const XV xv, double factor, const XV adder);


#endif // XV_H_INCLUDED
