#include <vector>
#include <array>
//#include <boost/array.hpp>
//#include <boost/multi_array.hpp>
#include <iostream>
#include <math.h>
#include <fftw3.h>
#include <complex>
#include "readfile.h"
#include <assert.h>

#ifndef HEADER
#define HEADER

#define GRAVITY (6.6738e-20*1.989e40)
#define NLIM    12
#define LLIM    6
#define NMAX    16
#define LMAX    10
#define NDIM    ((LLIM+1)*(LLIM+6)/2)

//#define SQRT4PI   3.5449077018110318
//#define SQRT4PI   3.544907701811032054596334966682290365595098912244774256427

typedef std::vector<double> state_type;
typedef std::vector<std::complex<double>>  complex_type;
typedef std::array<double,3> array3;
typedef std::array<double,6> array6;

extern Readparams allparams;

/* Structure for storing input params from Python init.py*/
struct Indata {
    double Knlm[LMAX][LMAX][NMAX][2];
    double scalerad, virialrad;
};

class calcAcc{
    protected:
    double temppot;
    struct Indata *var,*varfp;
    void cartToSph(state_type &x);
    void sphToCart(state_type &x);
    void cartVec(const state_type &x, state_type &vec);
    void getSphAcc(const struct Indata *Var, const state_type &sph, state_type &acc, const double t);
    void getSphPot(const struct Indata *Var, const state_type sph, double &pot);


    public:
        calcAcc() {};
        calcAcc(struct Indata *invar){var = invar; varfp=NULL;}//Nbody
        calcAcc(struct Indata *invar, struct Indata *invarfp){//Hydro
            assert(allparams.NumComponents == 2);
            var=invar;
            varfp=invarfp;
        }
        void getPotential(const array6 &x, double &pot);
        void getEnergy(const array6 &x, double &pot){
            temppot=0;
            getPotential(x, temppot);
            pot=temppot+0.5*(pow(x[3],2)+pow(x[4],2)+pow(x[5],2));
        }
        void getCartAcc(const array6 &x, array6 &dxdt, const double t);
        void operator() (const array6 &x, array6 &dxdt, const double t){
            getCartAcc(x,dxdt,t);
        }
};

class calcAccDecay: public calcAcc{
    protected:
    void getSphAcc(const struct Indata *Var, const state_type &sph, state_type &acc, const double t);

};

class calcOrb{
    private:
        calcAcc ACC;
        void getperiod(fftw_plan plan, const std::vector<array6>& Xpriv, state_type &v, 
                complex_type &fv, double &maxamp, int &maxi);

    public:
        calcOrb(struct Indata &invar, struct Indata &invarfp){
            ACC=(allparams.twocomp) ? calcAcc(&invar,&invarfp) : calcAcc(&invar);
        }
        void test(std::vector<array6> const &xinit);
        void integrate(state_type &XX, state_type &T, std::vector<array6> &xinit_run);
        void integrateDecay(state_type &XX, state_type &T, std::vector<array6> &xinit_run);
};

struct saveStates{
    std::vector<array6>& m_states;
    std::vector<double>& m_times;

    saveStates( std::vector<array6> &states ,std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const array6 &x , double t ){
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

//[ streaming_observer
struct streaming_observer
{
    std::ostream& m_out;
    streaming_observer( std::ostream &out ) : m_out( out ) { }

    template< class State >
    void operator()( const State &x , double t ) const
    {
        m_out << t;
        for( size_t i=0 ; i<x.size() ; ++i ) m_out << "\t" << x[i];
        m_out << "\n";
    }
};
//]

int loadHdf5Input(char *filename, struct Indata *var);
int loadHdf5Init(char *filename, std::vector<array6> &init);
int saveHdf5(state_type &OUT, state_type &TIME);

#endif

