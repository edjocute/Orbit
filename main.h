//#include "allvar.h"
#include <vector>


#ifndef HEADER
#define HEADER

#define GRAVITY 6.6738e-20*1.989e40
#define NLIM    12
#define LLIM    6
#define NMAX    14
#define LMAX    10

#define SQRT4PI 3.5449077018110318

typedef std::vector<double> state_type;

struct Indata {
    double Knlm[LMAX][LMAX][NMAX][2];
    double scalerad, virialrad;
};

class calcAcc{

    struct Indata *var;

    void cartToSph(state_type &x);
    void sphToCart(state_type &x);
    void cartVec(const state_type x, state_type &vec);

    public:
        calcAcc(struct Indata *invar){var = invar;}
        void getSphAcc(const state_type sph, state_type &acc);
        void getCartAcc(const state_type x, state_type &dxdt);
        void operator() (const state_type x, state_type &dxdt, const double){
            getCartAcc(x,dxdt);
        }
};


int loadHdf5Input(char *filename, struct Indata *var);
//void getAcceleration(const state_type sph, const struct Indata var, state_type &acc);



#endif

