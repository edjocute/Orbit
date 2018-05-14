#include <vector>
#include <array>
//#include <boost/array.hpp>
//#include <boost/multi_array.hpp>
#include <iostream>
#include <math.h>
#include <fftw3.h>
#include <complex>
#include "readfile.h"

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

extern Readparams allparams;

/* Structure for storing input params from Python init.py*/
struct Indata {
    double Knlm[LMAX][LMAX][NMAX][2];
    double scalerad, virialrad;
};

class calcAcc{
    double temppot;
    struct Indata *var,*varfp;
    void cartToSph(state_type &x);
    void sphToCart(state_type &x);
    void cartVec(const state_type &x, state_type &vec);
    void getSphAcc(const struct Indata *Var, const state_type &sph, state_type &acc);
    void getSphPot(const struct Indata *Var, const state_type sph, double &pot);

    public:
        calcAcc() {};
        calcAcc(struct Indata *invar){var = invar; varfp=NULL;}//Nbody
        calcAcc(struct Indata *invar, struct Indata *invarfp){//Hydro
            var=invar;
            varfp=invarfp;
        }
        void getPotential(const state_type x, double &pot);
        void getEnergy(const state_type x, double &pot){
            temppot=0;
            getPotential(x, temppot);
            pot=temppot+0.5*(pow(x[3],2)+pow(x[4],2)+pow(x[5],2));
        }
        void getCartAcc(const state_type x, state_type &dxdt);
        void operator() (const state_type x, state_type &dxdt, const double){
            getCartAcc(x,dxdt);
        }
};

class calcOrb{
    private:
        calcAcc ACC;
        //Readparams *allparams;
        void getperiod(fftw_plan plan, const std::vector<state_type>& Xpriv, state_type &v, 
                complex_type &fv, double &maxamp, int &maxi);

    public:
        //calcOrb(Readparams *params, struct Indata &invar, struct Indata &invarfp){
        calcOrb(struct Indata &invar, struct Indata &invarfp){
            //allparams=params;
            ACC=(allparams.twocomp) ? calcAcc(&invar,&invarfp) : calcAcc(&invar);
        }
        void test(std::vector<state_type> const &xinit);
        void integrate(state_type &XX, state_type &T, std::vector<state_type> &xinit_run);
};

struct saveStates{
    std::vector<state_type>& m_states;
    std::vector<double>& m_times;

    saveStates( std::vector<state_type> &states ,std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t ){
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
int loadHdf5Init(char *filename, std::vector<state_type> &init);
int saveHdf5(state_type &OUT, state_type &TIME);

#endif

