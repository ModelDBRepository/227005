 /**
	
	@author Pinar Oz <poz.neuro AT gmail.com>
	@author Michael Kreissl <mig80 AT gmx.net>


 */ 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "neuron.h"
#include "WB.h"
#include <ctime>
#include <cstdlib>
#include <errno.h>
/*! \mainpage The Multi-compartmental Cooperative Axon Initial Segment Model
 * \section intro Introduction
 * \section Cond Conductance Calculations
 	* \subsection WB Wang-Buzsaki Model
 	*\subsection Coop Implementation of Cooperativity
 *\section Neuron Neuron Construction
	*\subsection Build Neuron Geometry
	*\subsection Pot Membrane Potential Estimations
 */

using namespace std;

double gaussrand();//!< This is a random number generator. A pseudo random number generator is seeded from the clock and set in a new function which further randomizes the returned value!

/*! @file multicomp.cpp The units are :
mV
mum
nA
nF
ms
muS, MOhm
*/


int main(int argc, char *argv[])
{

    /*! \brief First part in the main code is the Neuron Construction. Passive electrical parameters are as below
      Simulation time =  100 ms.
      Time step = 10 microseconds;
      r_L = 5.0  KOhm.mm ( MOhm.mum );
      Cm = 1.0e-5  nF/mum^2;
      Cm_myelin = 1.0e-7   nF/mum^2;
      rL_myelin = 5.0 KOhm.mm ( MOhm.mum );
      a_in = 0.5  micron ; a_out = 0.5 micron	*/

      // ==== FIRST PART OF FILE NAME ===
      string root = "";	

      //This part only if these loops are commented out
      int xa = 25;
      int gFactor =3;
      // -----------

      int run_total = 1;
      double frequencies[]={ 0.001, 0.002, 0.005, 0.010, 0.020, 0.040, 0.060, 0.080, 0.100, 0.200, 0.300, 0.400, 0.500, 0.800}; 

/*	for(int xa = 25; xa<26; xa+=5)*/
	for(int r = 0; r<run_total; r++)
	{
	//	for(int gFactor = 3; gFactor<4; gFactor+=10)
	    for(int f_ind = 0; f_ind<1; f_ind++)
	    {                                                                                                                                                                                            
   		int V12 = 0;
		// for(int V12=0; V12<11; V12+=2)
		//   {
  		 // ==== 2ND PART OF FILE NAME ===
		stringstream temp;
  		temp << "filepath"<<gFactor;//<<"fI"<<frequencies[f_ind]*1000<<"_"<<r;
		      //"x" << xa<<"runtime_rloop_ctrl_gF"<<gFactor; //<<"_V12_"<<V12;
  		string prefix;
  		prefix = temp.str();

		      cout << "data will be written to " << prefix << "..." << endl;

		//ofstream fMO((prefix+"_MO.dat").c_str());
		//ofstream fMP((prefix+"_MP.dat").c_str());
		//ofstream fsoma((prefix+"_soma.dat").c_str());
		//ofstream flat((prefix+"_lat.dat").c_str());
		//ofstream fIinj((prefix+"_Iinj.dat").c_str());
 	 	ofstream fspike((prefix+"_spike.dat").c_str());
		//if (!fsoma.is_open())
		/*if(!fIinj.is_open())
		{
			cerr << "Can't open file for output.\n";
			exit(1);
		}*/

	  	double time1, tstart;      // time measurment variables
	  	tstart = clock();              // start

		// constants//!< micron
	  	double time=1000; // ms
	  	double dt = 0.01; //ms

	  	double V_init = E_WB[0];;

	  	double r_L = 5.0;           //**<  KOhm.mm =  MOhm.mum*/
		//r_L = 3.0;

	  	double Cm = 1.0e-5;          //< 1.e-5 nF/mum^2*/
	  	double Cm_myelin = 1.0e-7;   // nF/mum^2*/
	  	double rL_myelin = 5.0;	// <  KOhm.mm =  MOhm.mum*/

	  	double a_in = 0.5;          // < micron
  		double a_out = 0.5;		//< micron

 	 	double dudt_ap = 10.;

		//conductances, reversal energies
 	 	double g[3], dV0[2],dV1[2],dV2[2],g_dend[3], g_m[3], g_NR[3],g_AIS[3];
  				
			for (int i=0; i<3; i++)
	 		{
	    			g[i] =  g_bar_WB[i];
	    			g_NR[i] = g_bar_WB[i];
	    			g_AIS[i] = g_bar_WB[i];
	    			g_m[i] = 0.;
	    			g_dend[i] = 0.2*g_bar_WB[i];
				// unmyelined axon   g_m[i] = g_bar_WB[i];
	  		}

   		g_dend[1] *= 0.33;
		g_dend[2] *= 0.33;
		//g_m[0] = g_bar_WB[0]*0.4;
  	 	g_m[0] = g_bar_WB[0]*0.2;
 		g_NR[1]*= gFactor;
	 	g_NR[2] *= gFactor;
	  	g_AIS[1] *=gFactor;
	   	g_AIS[2] *=gFactor;
	  	/*g_AIS[0] = g_bar_WB[0];*/
		//g_NR[0] = g_bar_WB[0];
		//g[0] = g_bar_WB[0];
		double gl[3];
		gl[0]=g_bar_WB[0];
		gl[1]=gl[2]=0;
	  	
			for (int i=0; i<2; i++)
	  		{
	    			dV0[i] = 0.;
	    			dV1[i] = 0;//-V12;
	    			dV2[i] = 0.;
	  		}

	  	int parts = 0;
	  	int axon1_comp = 0;
	  	int axon2_comp = 0;

  			if(xa<1){
	    			parts = 26;
	    			axon2_comp =25;
	  			}
	  		else{
	    			parts = 27;
	    			axon1_comp = int(xa/2);
				// 	cout<<"#\t"<<axon1_comp<<endl;
	    			axon2_comp = int(25-xa/2);
	    			}

		double pi =0.1;
		double p0 = 0.;

 	  	neuron_part *neuronparts = new neuron_part[parts]; // name, #comp, length, a_in, a_out,Cm, r_L, ...
 	  	neuronparts[0].init("dendrite", 6, 50.,  0.5, 0.5, Cm, r_L, g_dend, dV0,p0);//g_dend 6
 	  	neuronparts[1].init("dendrite2", 10, 10., 0.5, 1., Cm, r_L,g_dend, dV0,p0);//g_dend 10
 	  	neuronparts[2].init("soma1", 10, 2.,  1., 10., Cm, r_L,g , dV0,p0);//g 10
 	  	neuronparts[3].init("soma2", 10,2.,  10., 1.5, Cm, r_L,g , dV0,p0);//g 10


 //	  	if(xa>1){
 	 		neuronparts[4].init("hillock",axon1_comp, 2.,  1.5,  1.5-axon1_comp*0.04,Cm, r_L, g, dV0,0);// g
// 	  	}

		  neuronparts[parts-22].init("proper",  axon2_comp, 2., 1.5-axon1_comp*0.04 , a_out, Cm, r_L, g_AIS, dV1,0);//g_AIS parts-22
	//   neuronparts[2].init("hillock0", 1, 2., a_out+0.04 , a_out , Cm, r_L, g_AIS, dV0,p0);
 	  	neuronparts[parts-21].init("myelin1",  5, 10., a_in, a_out, Cm_myelin, rL_myelin,g_m , dV0,p0);//
	  	neuronparts[parts-20].init("NR1",  1, 2., a_in, a_out, Cm, r_L, g_NR, dV2,p0);//g_NR
 	  	neuronparts[parts-19].init("myelin2",  5, 10., a_in, a_out, Cm_myelin, rL_myelin, g_m, dV0,p0);//
 	 	neuronparts[parts-18].init("NR2",  1, 2., a_in, a_out, Cm, r_L,g_NR, dV2,p0);//g_NR
	 	neuronparts[parts-17].init("myelin3",  5, 10., a_in, a_out, Cm_myelin, rL_myelin, g_m, dV0,p0);
 	 	neuronparts[parts-16].init("NR3",  1, 2., a_in, a_out, Cm, r_L,g_NR, dV2,p0);//g_NR
 	 	neuronparts[parts-15].init("myelin4",  5, 10., a_in, a_out, Cm_myelin, rL_myelin, g_m, dV0,p0);
 	 	neuronparts[parts-14].init("NR4",  1, 2., a_in, a_out, Cm, r_L,g_NR, dV2,p0);//g_NR
 	 	neuronparts[parts-13].init("myelin5",  5, 10., a_in, a_out, Cm_myelin, rL_myelin, g_m, dV0,p0);
 	 	neuronparts[parts-12].init("NR15",  1, 2., a_in,a_out, Cm, r_L,g_NR, dV2,p0);//g_NR
 	 	neuronparts[parts-11].init("myelin6",  5, 10., a_in, a_out, Cm_myelin, rL_myelin, g_m, dV0,p0);
 	 	neuronparts[parts-10].init("NR6",  1, 2., a_in, a_out, Cm, r_L,g_NR, dV2,p0);//g_NR
 	 	neuronparts[parts-9].init("myelin7",  5, 10., a_in, a_out, Cm_myelin, rL_myelin, g_m, dV0,p0);
 	 	neuronparts[parts-8].init("NR7",  1, 2., a_in, a_out, Cm, r_L,g_NR, dV2,p0);// g_NR
 	 	neuronparts[parts-7].init("myelin8",  5, 10., a_in, a_out, Cm_myelin, rL_myelin, g_m, dV0,p0);
 	 	neuronparts[parts-6].init("NR8",  1, 2., a_in, a_out, Cm, r_L,g_NR, dV2,p0);//g_NR
 	 	neuronparts[parts-5].init("myelin9",  5, 10., a_in, a_out, Cm_myelin, rL_myelin, g_m, dV0,p0);
 	 	neuronparts[parts-4].init("NR9",  1, 2., a_in, a_out, Cm, r_L,g_NR , dV2,p0);//g_NR
		neuronparts[parts-3].init("myelin10",  5, 10., a_in, a_out, Cm_myelin, rL_myelin, g_m, dV0,p0);
 	 	neuronparts[parts-2].init("bleb1",  3, 3.,  1., 3., Cm, r_L,g , dV0,p0);//g
 	 	neuronparts[parts-1].init("bleb2", 3, 3.,  3., 1., Cm, r_L,g , dV0,p0);//g

int comp_no = 126;
	  	neuron wholeneuron(parts, neuronparts, V_init, dt);

		//   ofstream fMO((prefix+"_MO.dat").c_str());
	//     neuron_part::print_param_names(fMO);

	  		for(int p=0; p<parts; p++)
	  		{
	   			neuronparts[p].print_paras();
		  //neuronparts[p].print_params(fMO);
	  		}
		

		  	unsigned long NtimeSteps = (unsigned long) (time/dt+1/*/2*/);
	  		cout << "maximal time steps are " << NtimeSteps << endl;
	  		cout << "dt [ms]: \t" <<  dt << endl << endl;
	  //create arrays to store
	  	double *timeStore = new double [NtimeSteps];

			for(unsigned long tStep=0; tStep < NtimeSteps; tStep++)
	    			timeStore[tStep] = 0.;
		

		double **mbStore = new double *[NtimeSteps];

			for(unsigned long tStep=0; tStep < NtimeSteps; tStep++){
				mbStore[tStep] = new double [comp_no];//wholeneuron.compartments()
				for(int c=0;c<comp_no;c++)
				{
					mbStore[tStep][c] = 0;
				}
			}
		double *NaI = new double [NtimeSteps];
		for(unsigned long tStep=0; tStep < NtimeSteps; tStep++)
			NaI[tStep] = 0; 

		double *SpikeTimeStore= new double [1001];
	         for(int x=0; x < 1001; x++)
			SpikeTimeStore[x] = 0; 
		 //*/
		// initialize the range for injected current
		double I_max =8e-3;
// 			cout<<"random\t"<<rand_mu<<"\t I_max \t"<<I_max<<endl;
	  	double I_min = 0.;
	  	double drate = 0.;

		double I0=0.5*(I_max+I_min); //0.015625  ;//0.0019;
	   	double sgI = 6e-3;
	 //	double sgI= 0.000085;  // nA
	//	double tau_c = 10;
		double tau_c=5;  //ms
		double fI = frequencies[f_ind];
	//	double f = 0.001; //kHz
		double I1;

		double expt  = exp(-dt/tau_c);
		double prefactor = sqrt(2*dt/tau_c);
	   	double Ix = 0;

		double Psi=-10;
	  //  	double I_inj
		double sc_rate;
		double flag_tick=0;
		double flag_tick2=0;
		double flag_out=0;
		unsigned long ticker;
		double check_pot;
		//double check_pot2;
		double mbcheck;
		double I_sin;

		srand(clock());

//   	      cout<<"Loop 1 Start"<<endl;
 	  		do
			{
				sc_rate=0;
				ticker=0;
	    			bool all_fired = false;
	    			wholeneuron.freq = 0.;
					double t=0.;
	    			unsigned long tStep = 0;

					for(unsigned long tStep=0; tStep < NtimeSteps; tStep++)
	      					timeStore[tStep] = 0.;
		      
		      			do
	      				{
  			      			t+=dt;
        
       
                        I1 = 0.4*I0;
 						I_sin = I1*cos(2*M_PI*fI*t);
    
						double randn = gaussrand();
						Ix = Ix*expt + prefactor*randn;//
   						double I_inj = I0 + sgI*Ix;// +I_sin;//

						neuronparts[3].current(I_inj, 1);
						NaI[tStep] = I_inj;

						  wholeneuron.timestep();

						check_pot = wholeneuron.pot[60];
					
 						
							
							if(flag_tick==0){

								if(check_pot>Psi){

									flag_tick++;
		                					SpikeTimeStore[ticker]=t;
									ticker++;
								}
							}

							if(flag_tick>0){

								if(check_pot<Psi){

									flag_tick=0;
								}
							}

 			    			timeStore[tStep] = t;

              				for(int i=0;i<comp_no;i++)
              						{
								mbStore[tStep][i] = wholeneuron.pot[i];
							}

		            			tStep++;


      					}while(t<time);

 						if(ticker>0)
							sc_rate =1000*((ticker-1)/t);
						else
							sc_rate = 0;
					cout << "ticker\t"<<ticker<<"\t time \t"<<t<<endl;
	      				cout << "spike frequency = " << sc_rate << "Hz" << endl;

				if(ticker>2){
					    if(sc_rate > 10)
					    {
						I_max = I0;
					    }else
						I_min=I0;

				 }
				if(ticker==1)
					    {
						I_max = I0;
					     }
				if(ticker<1){
			 			I_min = I0;
					}

				I0= 0.5*(I_max+I_min);
					cout<<I0<<endl;

				drate = sc_rate-10;
				drate = (drate >= 0) ? drate : -drate;

					cout << "timesteps = " << int(t/dt) << endl;
	   				cout << "sim time = " << t << " ms" << endl;

			} while(drate > 5);

	  	time1 = clock() - tstart;     // end
	  	time1 = time1/CLOCKS_PER_SEC;  // rescale to seconds
	  		cout << "cpu time = " << time1 << " s" << endl;
			cout<<"********************* run_no = "<< r+1 <<"\t frequency "<<fI<< " *********************"<<endl;
	        //save membrane potential to disk
	     if ( r+2>run_total){
	    	ofstream fMP((prefix+"_MP.dat").c_str());

	      		for(unsigned long tStep=0; tStep < NtimeSteps; tStep++)
	        	{
	          			if (timeStore[tStep] == 0.)
	            				break;

	          		fMP << timeStore[tStep];

	          			for(int i=0;i<comp_no;i++)//wholeneuron.compartments()
	              				fMP << " " << mbStore[tStep][i];

				fMP << endl;
		        }

 		fMP.close();  
		}


  			for(int i=0;i<ticker;i++){
				fspike <<SpikeTimeStore[i]<<endl;;
			}

		fspike.close();

	 
	   	delete[] timeStore;

				for(unsigned long tStep=0; tStep < NtimeSteps; tStep++)
	    			{
				delete[] mbStore[tStep];
			}

	  	delete[] mbStore;
		delete[] SpikeTimeStore;
	  	delete[] NaI;

		}

	}

return EXIT_SUCCESS;
}


double gaussrand()
{
  static double V1, V2, S;
  static int phase = 0;
  double X;

  	if(phase == 0) {

		do {
 		//srand(time(0));
			double U1 = (double)rand() / RAND_MAX;
      			double U2 = (double)rand() / RAND_MAX;

      			V1 = 2 * U1 - 1;
      			V2 = 2 * U2 - 1;
      			S = V1 * V1 + V2 * V2;

		} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
 	} else
    		X = V2 * sqrt(-2 * log(S) / S);

 phase = 1 - phase;
 return X;
}

