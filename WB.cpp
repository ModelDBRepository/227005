/**
	@author Pinar Oz <poz.neuro AT gmail.com>
  @author Michael Kreissl <mig80 AT gmx.net>
	
 */
#include "WB.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

WB::WB(){}

WB::~WB(){}

void WB::init(double V_init)
{
  V = V_init;
  m_tau_limit();
  m = limit;  
  h_tau_limit();
  h = limit;  
  n_tau_limit();
  n = limit; 
  mj_tau_limit();
  mj = limit;
  hj_tau_limit();
  hj = limit;
}

/// Calculating the gating probability for non-cooperative Na^+ gating (m^3*h)
double WB::calculate_Na(double V, double dt)
{
  this -> V = V;
   this -> dt = dt;
  calculate_m();
  calculate_h();
  return  pow(m,3)*h;
 }

/// Calculating the gating probability for cooperative Na^+ gating (m^3_{coop}*h_{coop})
double WB::calculate_Naj(double V, double dt) 
{
   this -> V = V;
   this -> dt = dt;
  calculate_mj();
  calculate_hj();
  return pow(mj,3)*h;//hj;

}

/// Calculating the gating probability for K^+ gating (n^4)
double WB::calculate_K(double V, double dt)
{this -> V = V;
   this -> dt = dt;
  calculate_n();
  return n*n*n*n; 
}

/// Calculating the Wang-Buzsaki type kinetics for classical "m"
inline void WB::m_tau_limit()
  { 
	if(V==-35)
{
		alpha=1; //=0.1/0.1;
	} else{
  alpha = 0.1*(V + 35.)/(1.-exp(-0.1*(V+35.)));
	}
  beta = 4*exp(-0.0556*(V + 60));
  
tau = tau_0/(alpha + beta); /// tau_0 is 1 in original.
limit =alpha/(alpha + beta);
}

/// Calculating the Wang-Buzsaki type kinetics for COOPERATIVE "m"
inline void WB::mj_tau_limit()
  { 
  if(KJ>0){
      VNa = V + (KJ*hcs*pow(mcs,3)); 
      if(V==-35){
	  alpha=0.1/0.1;
      }else{
	  alpha = 0.1*(VNa + 35.)/(1.-exp(-0.1*(VNa+35.)));
	}
      beta = 4*exp(-0.0556*(VNa + 60));
  }else{
      if(V==-35)
      {
	alpha=0.1/0.1;
      }else{
	alpha = 0.1*(V + 35.)/(1.-exp(-0.1*(V+35.)));
      }
      beta = 4*exp(-0.0556*(V + 60));
  }

tau = tau_0/(alpha + beta);
limit = alpha/(alpha + beta);
 // double theta = 1 + exp(-(VNa+35)/k_A);
//  limit = 1 / theta;  

}

/// Calculating the Wang-Buzsaki type kinetics for classical "h"
inline void WB::h_tau_limit()
{ 
  alpha = 0.35*exp(-0.05*(V+58.));
  beta = 5/(1. + exp(-0.1*(V + 28.)));
 
 tau = 1./(alpha + beta);
 limit =alpha/(alpha + beta);
// limit = tau * alpha;
}

///Calculating the Wang-Buzsaki type kinetics for COOPERATIVE "h"
inline void WB::hj_tau_limit()
{ 
 if(KJ>0){
    VNa = V + (KJ*hcs*pow(mcs,3));
    alpha = 0.35*exp(-0.05*(VNa+58.));
    beta = 5/(1. + exp(-0.1*(VNa + 28.)));
  }else{
    alpha = 0.35*exp(-0.05*(V+58.));
    beta = 5/(1. + exp(-0.1*(V + 28.)));
  }
   tau = 1/(alpha + beta);
//  double theta = 1+exp((VNa+58.)/k_A);
 // limit= 1 /theta ;
   limit =alpha/(alpha + beta);
   }

/// Calculating the Wang-Buzsaki type kinetics for classical "n"
inline void WB::n_tau_limit()
{ 
	if(V==-34)
{
		alpha=0.05/0.1;
	}else
{ 
		alpha = 0.05*(V + 34.)/(1. - exp(-0.1*(V + 34.)));
 	}
  beta = 0.625*exp(-0.0125*(V + 44));
  tau = 1./(alpha + beta);
  limit = tau * alpha;
}

/** Calculating the PDEs */

inline void WB::calculate_m()
{
  m_tau_limit();
   m = limit;
// pde_time(&m);
store_m(&ms,&tau_ms);
  
}

inline void WB::calculate_mj()
{
  mj_tau_limit();
 mj = limit;
//  pde_time(&mj);
store_mj(&mcs,&tau_mcs);
}

inline void WB::calculate_h()
{
  h_tau_limit();
  pde_time(&h);
store_h(&hs);
}

inline void WB::calculate_hj()
{
  hj_tau_limit();
  pde_time(&hj);
store_hj(&hcs);
}

inline void WB::calculate_n()
{
  n_tau_limit();
  pde_time(&n);
}

inline void WB::pde_time(double *x)
{
  *x = limit + (*x - limit) * exp(-dt/tau);
}

inline void WB::store_m(double *x, double*tx)
{
 *x = m;
 *tx = tau;
}
inline void WB::store_mj(double *x,double*tx)
{
 *x = mj;
 *tx = tau;
}
inline void WB::store_h(double *x)
{
 *x = h;
}

inline void WB::store_hj(double *x)
{
 *x = hj;
}

inline void WB::store_n(double *x)
{
 *x = n;
}