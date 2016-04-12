//#include "allvar.h"
#include <vector>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <iostream>


#ifndef HEADER
#define HEADER

#define GRAVITY 6.6738e-20*1.989e40
#define NLIM    12
#define LLIM    6
#define NMAX    14
#define LMAX    10

#define SQRT4PI 3.5449077018110318

typedef std::vector<double> state_type;
typedef boost::multi_array<double,3> store_array_type;

/* Structure for storing input params from Python init.py*/
struct Indata {
    double Knlm[LMAX][LMAX][NMAX][2];
    double scalerad, virialrad;
};

class calcAcc{
    struct Indata *var,*varfp;
    void cartToSph(state_type &x);
    void sphToCart(state_type &x);
    void cartVec(const state_type x, state_type &vec);
    void getSphAcc(const struct Indata *Var, const state_type sph, state_type &acc);

    public:
        calcAcc(struct Indata *invar){var = invar; varfp=NULL;}//Nbody
        calcAcc(struct Indata *invar, struct Indata *invarfp){//Hydro
            var=invar;
            varfp=invarfp;
        }
        void getCartAcc(const state_type x, state_type &dxdt);
        void operator() (const state_type x, state_type &dxdt, const double){
            getCartAcc(x,dxdt);
        }
};

struct saveStates
{
    std::vector<state_type>& m_states;
    std::vector<double>& m_times;

    saveStates( std::vector<state_type> &states ,std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
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
int saveHdf5(char *filename, state_type &OUT, state_type &TIME);
//void getAcceleration(const state_type sph, const struct Indata var, state_type &acc);

#endif

