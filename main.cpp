#include <math.h>
#include <stdio.h>
#include "main.h"


int main(){

    char filename[150];
    struct Indata infile; //to store coefficients
    state_type acc(6), sph(6);
    sph[0]=0.1; sph[1]=-0.2; sph[2]=0.3;

    snprintf(filename,150,"./var227.hdf5");
    loadHdf5Input(filename, &infile);

    calcAcc ACC (&infile);
    ACC(sph, acc,sph[0]);

    fprintf(stdout,"G=%e\n",GRAVITY);
    fprintf(stdout,"a=%e\n",infile.scalerad);
    fprintf(stdout,"acc=%e,%e,%e,%e,%e,%e\n",acc[0],acc[1],acc[2],acc[3],acc[4],acc[5]);
    fprintf(stdout,"Done!!\n");

    //write_output();

return 0;
}

