#include <math.h>
#include <stdio.h>
#include "main.h"
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

int main(){

    char filename[150];
    struct Indata infile; //to store coefficients
    state_type acc(6), sph(6);
    sph[0]=1.75151514e+17; sph[1]=-1.48955833e+16; sph[2]=-1.68645397e+17;
    sph[3]=3.83157600e+02; sph[4]=8.23211098e+01;  sph[5]=7.36204910e+01;

    snprintf(filename,150,"./var227.hdf5");
    loadHdf5Input(filename, &infile);

    fprintf(stdout,"G=%e\n",GRAVITY);
    fprintf(stdout,"a=%e\n",infile.scalerad);

    calcAcc ACC (&infile);
    ACC(sph, acc,sph[0]);
    fprintf(stdout,"acc=%e,%e,%e,%e,%e,%e\n",acc[0],acc[1],acc[2],acc[3],acc[4],acc[5]);


    double dt;
    std::cin >> dt;
    //typedef symplectic_rkn_sb3a_mclachlan< state_type > stepper_type;
    //typedef runge_kutta4< state_type > stepper_type;
    //typedef runge_kutta_cash_karp54< state_type > stepper_type;
    typedef runge_kutta_dopri5< state_type > stepper_type;
    //typedef runge_kutta_fehlberg78< state_type > stepper_type;
    integrate_const( stepper_type() , ACC , sph ,
                0.0 , 1e18 , dt , streaming_observer( std::cout ) );

    fprintf(stdout,"Done!!\n");
return 0;
}

