#ifndef NEURON_H
#define NEURON_H

/**
    @author Pinar Oz <poz.neuro AT gmail.com>
    @author Michael Kreissl <mig80 AT gmx.net>
    */
#include <iostream>
#include <string>
#include <cmath>
#include "WB.h"

using namespace std;

const double twopi = 2.*M_PI;
/*! \class neuron_part 
\brief The class where the construction of the single part of the neuron model is done.*/
/*! \class neuron 
\brief  The class where the actual calculations (i.e. membrane potential and difusion) are done.*/

class neuron_part
{

  friend class neuron;
  public:
    neuron_part();
    ~neuron_part();
    void init(string, int, double, double, double, double, double, double[3], double[2], double); 
    /*! \fn void init(string name, int comp, double l, double a_in, double a_out, double c_m, double r_l, double gi[3], double Vi[2],double pr)
    \brief Initiates all the parameters for the part. Only new thing here is the parameter pr.

    	\param name file name and place of storage
    	\param comp number of comaprtments
	\param l length of each compartment
	\param a_in initial radius of the part
	\param a_out final radius of the part
	\param c_m membrane capacitance
	\param c_m Intracellular resistivity
	\param gi[3] the conductance of leak, Na and K channels (see WB.h)
	\param Vi[2] the voltage shift (only when a voltage shift in the channel gating is used)
	\param pr fraction (percentage) of the coupling channels among the whole population in a given area. Always between 0 and 1.
*/

    void print_paras();
    void print_params(ostream& str);
    static void print_param_names(ostream& str);
    
/*!\brief  The function that takes the values for the injected current and the compartment of injection

\param I injected current (nA)
\param c site of current injection
*/
void current (double, int); 
    double pr;                
  // double gj;
private:
    string name;
    int comp;                   
    double l;                  
    double c_m; /**< membrane capacitance, nF/mum^2*/ 
    double r_l; /**<Intracellular resistivity , KOhm.mm ( MOhm.mum );*/
    double *gi, *Vi, Ei[3]; 

    double I_inj;              	
    int comp_inj;		
    double *R_L1,*R_L2,*area,*radius_in; ///Resistance, area, radii of the compartment
};

class neuron
{
  public:
    neuron(int, neuron_part*, double, double);
    ~neuron();
    void timestep();

    int compartments();
    bool check_latency(double, double);
    double *pot, *dudt, *latency;
    double freq;
    double ap[4];
    double *lengths;
    double *i_Na;		/**<Sodium current storage*/
    double *i_lateral;		/**<Lateral current storage*/
  //  int comp_soma;		
    double *mcstore,*mstore,*tau_mstore,*tau_mcstore;

  private:
    int parts, comps,comp_cut,comp_cutend;
    neuron_part *neuronparts;
    inline void diffusion_implicit();  /**< Function that calculates the diffusion among compartments*/
    double V_init, dt;
    double *g_plus, *g_minus; /**<Flux to the neighboring compartments*/
    double *c_m, *i_e;
    double *alpha, *beta, *gamma, *delta; /**< Abbrevations used in the calculations for messy sums and products*/
    double *dudt_last, dudt_max,*pot_last;
    int max, comp_lat;
    WB *membrane;
    double g[3], sum_gi, sum_giEi,sum_L,sum_K,sum_Na; 
};

#endif
