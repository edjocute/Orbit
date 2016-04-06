#include <math.h>
#include <stdio.h>
#include "main.h"


int main(){

    char filename[150];
    struct Indata var;
    double acc[3];
    //double sph[3];

    snprintf(filename,150,"var227.hdf5");
    load_hdf5_input(filename, &var);

    double sph[3]={0.1,-0.2,0.3};
    get_acceleration(sph, var, acc);

    fprintf(stdout,"G=%e\n",GRAVITY);
    fprintf(stdout,"a=%e\n",var.scalerad);
    fprintf(stdout,"acc=%e,%e,%e\n",acc[0],acc[1],acc[2]);
    fprintf(stdout,"Done!!\n");



    /*
    for (int i=0; i<20;i++){
        for (int j=0; j<10; j++){
            for (int m=0; m<10; m++){
                fprintf(stdout,"%d,%d,%d,Knlm=%f,%f\n",
                        i,j,m,var.Knlm[i][j][m][0],var.Knlm[i][j][m][1]);
            }
        }
    }
    fprintf(stdout,"Scale radius=%f\n",var.scalerad);
    fprintf(stdout,"R_200c=%f\n",var.virialrad);
    */

    //integrate();
    //write_output();


return 0;
}

