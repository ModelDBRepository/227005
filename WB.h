#ifndef WB_H
#define WB_H
/**
    @author Pinar Oz <poz.neuro AT gmail.com>
    @author Michael Kreissl <mig80 AT gmx.net>

 */
#include <cmath>

using namespace std;

/*! leak, Na, K */

/*! /muS/mum^2 */
static const double g_bar_WB[3] = {1e-6, 0.3e-3, 0.15e-3};  

static const double E_WB[3] = {-65., 55., -90.};  /*!original values (mV ) from paper. 
/**!This part will be used for the implementation of the cooperativity */
static const double KJ=900.; /*! coupling strength in mV*/
static const double k_A= 4; /*! used only in the V-clamp mode*/
static const double tau_0=0.1; /*! a constant, used in \tau_m*/


class WB
{
public:
  
WB();
  ~WB(); 
  void init(double);
  double calculate_Na(double, double);
  double calculate_K(double, double); 

/*! these two are new!*/  
double calculate_Naj(double, double); /// Cooperative !
  double mcs, hcs,ms,hs,ns,tau_mcs,tau_ms;

private:

  inline void pde_time(double*); 
  inline void store_h(double*);
  inline void store_m(double*,double*);
  inline void store_n(double*);
  double V, dt,t;
  double m,h,n; 

  double alpha, beta, tau, limit;   //temporal use
  inline void m_tau_limit();
  inline void h_tau_limit();
  inline void n_tau_limit();

  inline void calculate_m();
  inline void calculate_h();
  inline void calculate_n();


  inline void store_mj(double*,double*);
  inline void store_hj(double*);
  inline void mj_tau_limit();
  inline void hj_tau_limit();
  inline void calculate_mj();
  inline void calculate_hj();
  double VNa;
  double mj,hj;
 };

#endif
