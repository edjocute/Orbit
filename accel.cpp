#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_legendre.h>
#include "main.h"

/*Calculates unnormalized accelerations from positions in spherical coordinates
 * Requires normalization by G/a**2
 */
int get_acceleration(double sph[3], struct Indata var, double *acc){

    int ndim=(LLIM+2)*(LLIM+3)/2;//array dimension required for to Legendre polynomials
    double cosm[LLIM],sinm[LLIM],legen[ndim],legendiff[ndim];
    double gegen[NLIM],gegenm1[NLIM],Phi[NLIM],Phidiff[NLIM];
    double CDEF[4];
    int z;

    const double sintheta=sqrt(1-pow(sph[1],2));
    const double sinthetainv=1/sintheta;
    const double r=sph[0];
    const double xi=(r-1)/(r+1);

    /*
    fprintf(stdout,"%f,%f,%f\n",sph[0],sph[1],sph[2]);
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


    /* Initialize required arrays */
    acc[0]=0;acc[1]=0;acc[2]=0;

    for (int m=0;m<LLIM;m++){
        cosm[m]=cos(m*sph[2]);
        sinm[m]=sin(m*sph[2]);
    }

    /* Calculate generalized Legendre polynomials P_lm(x).
     * Conventions tested to be same as Python implementation.
     * This is an array containing non-zeros elements with 0<l<llim and 0<m<=l.
     * To access array elements, use index=l(l+1)/2+m */
    gsl_sf_legendre_deriv_array_e(GSL_SF_LEGENDRE_NONE,LLIM,sph[1],-1,legen,legendiff);

    /* Start loops over l */
    fprintf(stdout,"starting loops... \n");
    for (int ll=0;ll<LLIM;ll++){

        /* Calculate ultraspherical harmonics (Gegenbauer polynomials)
         * in array of length nlim for 0<n<nlim */
        gsl_sf_gegenpoly_array(NLIM,2*ll+1.5,xi,gegen);
        gsl_sf_gegenpoly_array(NLIM-1,2*ll+2.5,xi,gegenm1);
        double Phifac = -1*SQRT4PI*pow(r,ll)*pow(1+r,-2*ll-1);

        Phi[0]=Phifac;
        Phidiff[0]=Phifac*(ll/r-(2*ll+1)/(1+r));
        for (int n=1;n<NLIM;n++){
            Phi[n]=Phifac*gegen[n];
            Phidiff[n]=Phifac * ( (8*ll+6)*pow(1+r,-2)*gegenm1[n-1] + (ll/r-(2*ll+1)/(1+r))*gegen[n]);
        }
        
        for (int mm=0;mm<=ll;mm++){

            /* Reinitialize C_lm, D_lm etc. to 0 for each l,m 
             * and sum over n*/
            CDEF[0]=0; CDEF[1]=0; CDEF[2]=0; CDEF[3]=0;

            for (int n=0;n<NLIM;n++){
                CDEF[0]+=Phi[n]     *var.Knlm[n][ll][mm][0];
                CDEF[1]+=Phi[n]     *var.Knlm[n][ll][mm][1];
                CDEF[2]+=Phidiff[n] *var.Knlm[n][ll][mm][0];
                CDEF[3]+=Phidiff[n] *var.Knlm[n][ll][mm][1];
                }
            
            z = ll*(ll+1)/2 +  mm;
            //fprintf(stdout,"%d,%d,%f,%f\n",ll,mm,legen[z],legendiff[z]);
            acc[0]-=legen[z] * (CDEF[2]*cosm[mm] + CDEF[3]*sinm[mm]);
            acc[1]+=sintheta/r * legendiff[z]* (CDEF[0]*cosm[mm] + CDEF[1]*sinm[mm]);
            if (mm!=0){
                acc[2]-= mm*legen[z]/(sintheta*r)*(CDEF[1]*cosm[mm] - CDEF[0]*sinm[mm]);
            }
        }
    }

    /* Normalize the accelerations */
    for (int i =0;i<3;i++){
        acc[i]*=GRAVITY*pow(var.scalerad,-2);
    }

    fprintf(stdout,"Loops done\n");
    return 0;
}




                






