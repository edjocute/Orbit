#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include "main.h"

#define SQRT4PI (2*M_SQRTPI)

/* Converts spherical to cartesian coordinates */
void calcAcc::sphToCart(state_type &x){
    double sintheta=sqrt(1-pow(x[1],2));
    double xx = x[0]*sintheta*cos(x[2]);
    double yy = x[0]*sintheta*sin(x[2]);
    double zz = x[0]*x[1];
    x[0]=xx; x[1]=yy; x[2]=zz;
}

/* Converts cartesian coordinates to spherical coordinates */
inline void calcAcc::cartToSph(state_type &x){
            double r = gsl_hypot3(x[0],x[1],x[2]);
            if (r==0){
                x[0]=0; x[1]=1e-5; x[2]=1e-5;
            }  
            else{
                double costheta=x[2]/r;
                double phi=atan2(x[1],x[0]);
                x[0]=r;
                x[1]=costheta;
                x[2]=phi;
            }
}

/* Convert vector in spherical coordinates to cartesian coorinates */
void calcAcc::cartVec(const state_type x, state_type &vec){
    double costheta=x[1];
    double sintheta=sqrt(1-pow(costheta,2));
    double cosphi=cos(x[2]);
    double sinphi=sin(x[2]);

    state_type temp(3);
    temp[0] = sintheta*cosphi*vec[0] + costheta*cosphi*vec[1] - sinphi*vec[2];
    temp[1] = sintheta*sinphi*vec[0] + costheta*sinphi*vec[1] + cosphi*vec[2];
    temp[2] = costheta*vec[0] - sintheta*vec[1];
    vec[0]=temp[0]; vec[1]=temp[1]; vec[2]=temp[2];
}



/* Wrapper around getSphAcc to return dx/dt=(x',x'') from x=(x,v).
 * Reads and return accelerations in cartesian coordinates
 * Allows for both 1 and 2 component models
 */
void calcAcc::getCartAcc(const state_type x, state_type &dxdt){
    dxdt[0]=x[3]; dxdt[1]=x[4]; dxdt[2]=x[5]; //write velocities into dxdt
    state_type pos(3),acc(3);

    /* Calculate acceleration for 1st component (i.e. DM) */
    double ainv=1./var->scalerad;
    pos[0]=x[0]*ainv; pos[1]=x[1]*ainv; pos[2]=x[2]*ainv;
    cartToSph(pos); //convert positions to spherical coordinates
    getSphAcc(var,pos,acc); //calculate acceleration in spherical coordinates
    cartVec(pos,acc); //convert accelerations to cartesian
    dxdt[3]=acc[0]; dxdt[4]=acc[1]; dxdt[5]=acc[2];

    /* Do 2nd component if it exists */
    if (varfp){
        acc[0]=0;   acc[1]=0;   acc[2]=0; //reinitialize acc to zeros
        ainv=1./varfp->scalerad;
        pos[0]=x[0]*ainv; pos[1]=x[1]*ainv; pos[2]=x[2]*ainv;
        cartToSph(pos);
        getSphAcc(varfp,pos,acc);
        cartVec(pos,acc);
        dxdt[3]+=acc[0]; dxdt[4]+=acc[1]; dxdt[5]+=acc[2];
    }
}

/*Calculates unnormalized accelerations from positions in spherical coordinates
 * Requires normalization by G/a**2
 */
void calcAcc::getSphAcc(const struct Indata *Var, const state_type sph, state_type &acc){

    //int ndim=(LLIM+2)*(LLIM+3)/2;//array dimension required for Legendre polynomials
    //size_t ndim=gsl_sf_legendre_array_n(LLIM);
    int z; //index for Legendre polynomials
    double cosm[LLIM],sinm[LLIM],legen[NDIM],legendiff[NDIM];
    double gegen[NLIM] __attribute__((aligned(32))), gegenm1[NLIM]  __attribute__((aligned(32)));
    double Phi[NLIM]  __attribute__((aligned(32))),Phidiff[NLIM]  __attribute__((aligned(32)));
    double CDEF[4],Phifac,Phidiffac1,Phidiffac2;

    const double sintheta=sqrt(1.-gsl_pow_2(sph[1]));
    const double r=sph[0];
    const double xi=(r-1.)/(r+1.);

    /* Initialize required arrays */
    #pragma omp simd
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
    for (int ll=0;ll<LLIM;ll++){

        /* Calculate ultraspherical harmonics (Gegenbauer polynomials)
         * in array of length nlim for 0<n<nlim */
        gsl_sf_gegenpoly_array(NLIM,2*ll+1.5,xi,gegen);
        gsl_sf_gegenpoly_array(NLIM-1,2*ll+2.5,xi,gegenm1);
        Phifac = -1*SQRT4PI*gsl_pow_int(r,ll)/gsl_pow_int(1.0+r,2*ll+1);
        Phidiffac1 = (8*ll+6)/gsl_pow_2(1.0+r);
        Phidiffac2 = (ll/r-(2*ll+1)/(1.0+r));

        //Phi[0]=Phifac;
        Phidiff[0]=Phifac*Phidiffac2;
        /*for (int n=1; n<NLIM; n++){
            Phi[n]=Phifac*gegen[n];
            //Phidiff[n]=Phifac * ( (8*ll+6)/gsl_pow_2(1+r)*gegenm1[n-1] + (ll/r-(2*ll+1)/(1+r))*gegen[n]);
            Phidiff[n]=Phifac* (Phidiffac1*gegenm1[n-1] + Phidiffac2*gegen[n]);
        }*/
        Phi[:]=Phifac*gegen[:];
        Phidiff[1:NLIM]=Phifac*(Phidiffac1*gegenm1[:] + Phidiffac2*gegen[1:NLIM]);

        for (int mm=0;mm<=ll;mm++){
            /* Reinitialize C_lm, D_lm etc. to 0 for each l,m 
             * and sum over n*/
            //CDEF[0]=0; CDEF[1]=0; CDEF[2]=0; CDEF[3]=0;
            CDEF[:]=0;

            #pragma novector
            for (int n=0; n<NLIM; n++){
                CDEF[0]+=Phi[n]     *Var->Knlm[ll][mm][n][0];//C
                CDEF[1]+=Phi[n]     *Var->Knlm[ll][mm][n][1];//D
                CDEF[2]+=Phidiff[n] *Var->Knlm[ll][mm][n][0];//E
                CDEF[3]+=Phidiff[n] *Var->Knlm[ll][mm][n][1];//F
                }
            
            z = (ll*(ll+1))/2 +  mm;
            acc[0]-= legen[z] *     (CDEF[2]*cosm[mm] + CDEF[3]*sinm[mm]);
            acc[1]+= legendiff[z]*  (CDEF[0]*cosm[mm] + CDEF[1]*sinm[mm]);
            if (mm!=0){
            acc[2]-= mm*legen[z]*   (CDEF[1]*cosm[mm] - CDEF[0]*sinm[mm]);
            }
        }
    }

    /* Normalize the accelerations */
    double norm=GRAVITY/gsl_pow_2(Var->scalerad);
    acc[0]*= norm;
    acc[1]*= norm*sintheta/r;
    acc[2]*= norm/(sintheta*r);
}


void calcAcc::getPotential(const state_type x, double &pot){
    state_type pos(3);
    double temppot;

    /* Calculate potential for 1st component (i.e. DM) */
    double ainv=1./var->scalerad;
    pos[0]=x[0]*ainv; pos[1]=x[1]*ainv; pos[2]=x[2]*ainv;
    cartToSph(pos); //convert positions to spherical coordinates
    getSphPot(var,pos,temppot); //calculate potential
    pot=temppot;

    /* Do 2nd component if it exists */
    if (varfp){
        ainv=1./varfp->scalerad;
        pos[0]=x[0]*ainv; pos[1]=x[1]*ainv; pos[2]=x[2]*ainv;
        cartToSph(pos);
        getSphPot(varfp,pos,temppot);
        pot+=temppot;
    }
}

void calcAcc::getSphPot(const struct Indata *Var, const state_type sph, double &pot){
    //int ndim=(LLIM+2)*(LLIM+3)/2;//array dimension required for to Legendre polynomials
    double cosm[LLIM],sinm[LLIM],legen[NDIM],legendiff[NDIM];
    double gegen[NLIM],Phi[NLIM];
    double C,D,Phifac;
    int z;

    const double sintheta=sqrt(1.-gsl_pow_2(sph[1]));
    const double r=sph[0];
    const double xi=(r-1.)/(r+1.);

    /* Initialize required arrays */
    pot=0;
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
    for (int ll=0;ll<LLIM;ll++){

        /* Calculate ultraspherical harmonics (Gegenbauer polynomials)
         * in array of length nlim for 0<n<nlim */
        gsl_sf_gegenpoly_array(NLIM,2*ll+1.5,xi,gegen);
        Phifac = -1*SQRT4PI*gsl_pow_int(r,ll)/gsl_pow_int(1+r,2*ll+1);

        Phi[0]=Phifac;
        Phi[1:NLIM]=Phifac*gegen[1:NLIM];
        /*for (int n=1; n<NLIM; n++){
            Phi[n]=Phifac*gegen[n];
        }*/
        
        for (int mm=0;mm<=ll;mm++){
            /* Reinitialize C_lm, D_lm etc. to 0 for each l,m 
             * and sum over n*/
            C=0; D=0;

            #pragma novector
            for (int n=0; n<NLIM; n++){
                C +=Phi[n]     *Var->Knlm[ll][mm][n][0];
                D +=Phi[n]     *Var->Knlm[ll][mm][n][1];
            }
            
            z = (ll*(ll+1))/2 +  mm;
            pot += legen[z] * (C*cosm[mm] + D*sinm[mm]);
        }
    }
    /* Normalize the potentials */
    pot *= GRAVITY/Var->scalerad;
}
