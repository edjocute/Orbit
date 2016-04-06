//#include "allvar.h"

#ifndef HEADER
#define HEADER



#define GRAVITY 6.6738e-20*1.989e40
#define NLIM    12
#define LLIM    6
#define SQRT4PI 3.5449077018110318


struct Indata {
    double Knlm[20][10][10][2];
    double scalerad, virialrad;
};

int load_hdf5_input(char *filename, struct Indata *var);
int get_acceleration(double sph[3], struct Indata var, double *acc);



#endif

