/**
  @author Pinar Oz <poz.neuro AT gmail.com>
  @author Michael Kreissl <mig80 AT gmx.net>


 */
 #include "neuron.h"
#include "WB.h"
#include<assert.h>


/*! \file neuron.cpp
    \brief All the individual potential estimations are written here.
    
    Everything is as in original if not mentioned otherwise.
*/ 

neuron_part::neuron_part()
{
 I_inj = 0.;
 comp_inj = 0;
}

neuron_part::~neuron_part()
{
 delete[] area;
 delete[] R_L1;
 delete[] R_L2;
 delete[] radius_in;
}

void neuron_part::init(string name, int comp, double l, double a_in, double a_out, double c_m, double r_l, double gi[3], double Vi[2],double pr) 
{
/*! The parameters that are defined in neuron.h are given their values that are called from multicomp.cpp.
\par Defining the neuron geometry
The model geometry is given in this part. There are two alternatives so far that were used. In the current model the octagon shaped soma is used. (M_PI is pi value.)
\code
		radius_in[x] = a_in+x*(a_out-a_in)/comp;        
      		radius_out = a_in+(x+1)*(a_out-a_in)/comp;
      		radius_add = radius_in[x]+radius_out;
      		area[x] = M_PI*radius_add*l;
      
\endcode 
	The second alternative is not used anymore.
\code
		double Lhalf = l*comp/2;
    		double xi = x-Lhalf;
    			 if(a_max>a_in)
     			{
      				radius_in[x] = sqrt(a_max*a_max-xi*xi*(a_max*a_max-a_in*a_in)/Lhalf/Lhalf);  
				radius_out = sqrt(a_max*a_max- (xi+1)*(xi+1)*(a_max*a_max-a_in*a_in)/Lhalf/Lhalf); 
     			}
     			else
     			{
       				radius_in[x] = a_in;
       				radius_out=a_in;  
     			}
\endcode	

\par
Calculating the radii enables us to estimate the longtidunal resistance per compartment : 
\code
		R_L1[x] = r_l*l/M_PI/radius_in[x]/radius_add;
     		R_L2[x] = r_l*l/M_PI/radius_out/radius_add;
\endcode
*/

 this -> name = name;
 this -> comp = comp;
 this -> l    = l;
 this -> c_m  = c_m;
 this -> r_l  = r_l;
  /*(channels.find("L")!=string::npos) ? flag[0] = true : flag[0] = false;
  (channels.find("N")!=string::npos) ? flag[1] = true : flag[1] = false;
  (channels.find("K")!=string::npos) ? flag[2] = true : flag[2] = false;
  (channels.find("A")!=string::npos) ? flag[3] = true : flag[3] = false;*/
 this -> gi  = gi;
 this -> Vi  = Vi;
 this -> pr  = pr;
	for(int i=0; i<3; i++)
    		Ei[i] = E_WB[i];
  
 area = new double [comp];
 R_L1 = new double [comp];
 R_L2 = new double [comp];
 radius_in = new double [comp];
 double radius_out;
 double radius_add;
  
	for (int x=0; x<comp; x++)
  	{    
//    		double Lhalf = l*comp/2;
//    		 double xi = x-Lhalf;
//    			 if(a_max>a_in)
//     			{
//      			radius_in[x] = sqrt(a_max*a_max- 		xi*xi*(a_max*a_max-a_in*a_in)/Lhalf/Lhalf);  
//				radius_out = sqrt(a_max*a_max- (xi+1)*(xi+1)*(a_max*a_max-a_in*a_in)/Lhalf/Lhalf); 
//     			}
//     			else
//     			{
//       			radius_in[x] = a_in;
//       			radius_out=a_in;  
//     			}
    
      		radius_in[x] = a_in+x*(a_out-a_in)/comp;        
      		radius_out = a_in+(x+1)*(a_out-a_in)/comp;
      		radius_add = radius_in[x]+radius_out;
      		area[x] = M_PI*radius_add*l;
      		R_L1[x] = r_l*l/M_PI/radius_in[x]/radius_add;
     		R_L2[x] = r_l*l/M_PI/radius_out/radius_add;
    //			cout<< "area: \t" << radius_in[x] <<"\tR_L\t" <<R_L1[x] << "\t" << R_L2[x] << endl;
     //		if(x>comp-2){cout<<radius_in[0]<<"\t"<<x<<"\t"<<radius_in[x]<<"\t"<<radius_out<<endl;}
   	}
   
}

void neuron_part::print_param_names(ostream& str)
{
  str << "#comp l radius_in area" << endl; 
}               
                
void neuron_part::print_params(ostream& str)
{
  	for (int x=0; x<comp; x++)
  		str << x+1 << " " << l  << " " <<  radius_in[x]  << " " << area[x] << endl; 
 }               
   
                
void neuron_part::print_paras()
{               
  	cout << "name: \t\t" << name << endl;
  /*	cout << "comp: \t\t" << comp << endl;
  	cout << "l: \t\t" << l << " mikron" << endl;
  	cout << "radius_in[0]: \t" << radius_in[0] << " mikron" << endl;
  	cout << "radius_in["<<comp-1<<"]: \t" << radius_in[comp-1] << " mikron" << endl;
	cout << "c_m: \t\t" << c_m << " nF/mikron^2" << endl;
  	cout << "r_l: \t\t" << r_l << " MOhm mikron" << endl;
  	cout << "gi: \t\t";
  		for (int i=0; i<3; i++)
    			cout << gi[i] << " ";
  	cout << endl; 
   	cout<< "cooperativity fraction: \t\t"<< pr << endl;
  	cout << "Ei: \t\t";
  		for (int i=0; i<3; i++)
  			cout << Ei[i] << " ";
 	cout << endl; 
  	cout << "dVi: \t\t 0 ";
  		for (int i=0; i<2; i++)
    			cout << Vi[i] << " ";
  	cout << endl;   
  	cout << "I_inj [nA]: \t" << I_inj << endl;
  	cout << "comp_inj: \t" << comp_inj << endl << endl;
  
  //	cout << "name\t #comps\t length of comp\t radius in\t radius out\t c_m\t r_l\tg_L\tg_Na" << endl;*/
  
}



void neuron_part::current (double I, int c)
{
  I_inj = I;
  comp_inj = c-1;               // -1 because it starts with 0
}


neuron::neuron(int partsinneuron, neuron_part *nxxx, double initialV, double tstep) : parts(partsinneuron), neuronparts(nxxx), V_init(initialV), dt(tstep)
{
  

comps = 0;


  	for (int p=0; p<parts; p++)
   	{
     		comps += neuronparts[p].comp;
   
     		/*	if(p == 3) 
       				comp_soma = comps;*/
			
     
      	}



	cout << "The neuron consist of " << parts << " parts and " <<  comps << " compartments." << endl ;
  	//cout << "The soma compment is " << comp_soma <<endl;
  	// If a compartment is cut!!
// 	cout<<"The cut is at\t" <<comp_cut<<endl;

  g_plus = new double [comps];
  g_minus = new double [comps];
  
  c_m = new double [comps];
  i_e = new double [comps];
  
  alpha = new double [comps];
  beta = new double [comps];
  gamma = new double [comps];
  delta = new double [comps];
    
  pot = new double [comps];  
  
  dudt = new double [comps];
  dudt_last = new double [comps];
  ap[0]=ap[1]=ap[2]=ap[3]=0.;
  freq = 0.;
  latency = new double[comps];
  
  lengths = new double[comps];
  i_Na = new double[comps];
  i_lateral = new double[comps];
  
  membrane = new WB[comps];
  
  comp_lat =0;
  

  double R[comps];//, area[comps];
  R[0] = neuronparts[0].R_L1[0];
  lengths[0] = neuronparts[0].l;
  c_m[0] = neuronparts[0].c_m/dt;

/*  int x = 1;						//compartments in whole neuron
  for(int p=0; p<parts; p++)				//parts in whole neuron
  {

    for(int px=1; px<neuronparts[p].comp; px++)		//compartments in parts
    {
      R[x] = neuronparts[p].R_L2[px-1] + neuronparts[p].R_L1[px];

      lengths[x] = neuronparts[p].l;
      c_m[x] = neuronparts[p].c_m/dt;

      x++;
    }
  
    if (p+1 < parts)
    {
      R[x] = neuronparts[p].R_L2[neuronparts[p].comp-1] + neuronparts[p+1].R_L1[0];

      lengths[x] = neuronparts[p].l;
      c_m[x] = neuronparts[p].c_m/dt;

      x++;
    }
    
  }*/

 

  	for(int p=0, px=0, x=0; x<comps-1; x++, px++) /// if there is only one compartment x should be smaller than "comps" otherwise comps-1
  	{
	    		assert(px>=0);
    			if (px<neuronparts[p].comp-1) ///next part of neuron not reached
      				R[x+1] = neuronparts[p].R_L2[px]+neuronparts[p].R_L1[px+1];
    			else{ 
      				R[x+1] = neuronparts[p].R_L2[px]+neuronparts[p+1].R_L1[0];
      				p++; px=-1;
    			}
    
///Note here that membrane capacitance is normalized by the time interval.
    			if (px>=neuronparts[p].comp) {p++; px=0;}
    //				cout<< "x=: \t" << x <<"\tp="<<p<<"\tpx="<<px<<"\tR[x+1]=" <<R[x+1] << endl;
    		lengths[x] = neuronparts[p].l;
    		c_m[x] = neuronparts[p].c_m/dt;
  	}


	//cummulative sum of the lengths
  	for(int x=1; x<comps; x++)
    		lengths[x] = lengths[x] + lengths[x-1];


/*  x = 0;						//compartments in whole neuron
  for(int p=0; p<parts; p++)				//parts in whole neuron
    for(int px=0; px<neuronparts[p].comp; px++)		//compartments in parts
    {
      if (x<comps-1)
      {
	g_minus[x] = 1./R[x]/neuronparts[p].area[px];
	g_plus[x] = 1./R[x+1]/neuronparts[p].area[px];
        x++;
      }
    }
*/
  //set values for first and last compartment
 /* g_minus[0] = 0.;
  g_plus[comps-1] = 0.;


    g_plus[0] = 1./R[1]/neuronparts[0].area[0];
    g_minus[comps-1] = 1./R[comps-2]/neuronparts[parts-1].area[neuronparts[parts-1].comp-1];
  */



  		for(int p=0, px=1, x=1; x<comps-1; x++, px++)
  	{  
    		if (px>=neuronparts[p].comp) {p++; px=0;}
    	
		g_plus[x] = 1./R[x+1]/neuronparts[p].area[px];
    		g_minus[x] = 1./R[x]/neuronparts[p].area[px];
  	}
	g_minus[0] = 0.;
	g_plus[0] = 1./R[1]/neuronparts[0].area[0];
	g_minus[comps-1] = 1./R[comps-2]/neuronparts[parts-1].area[neuronparts[parts-1].comp-1];
 	g_plus[comps-1] = 0.;



  	for(int x=0; x<comps; x++)
	{
		alpha[x] = -g_minus[x];
    		gamma[x] = -g_plus[x];
		
		pot[x] = V_init;
		membrane[x].init(V_init);
    		dudt_last[x] = 0.;
    //			cout<< "alpha: \t" << alpha[x] <<"\tgamma\t" <<gamma[x] << endl;
     	}


/* gamma[47] = 0.;
 alpha[48] = 0.;
 gamma[48] = 0.;
 alpha[49]=0.;*/


/*!  \par Total Resistance
First thing to calculate for each comaprtment one after another is the total resistance between two comaprtments
\code
 double R[comps];
   R[0] = neuronparts[0].R_L1[0];
  
  	for(int p=0, px=0, x=0; x<comps; x++, px++) // If there is only one compartment x should be smaller than "comps" otherwise comps-1
  	{
		if(comps>1){ // This piece is only for multi-compartmental cases. (see Warning)
  	  		assert(px>=0);
    			if (px<neuronparts[p].comp-1) // next part of neuron not reached
      				R[x+1] = neuronparts[p].R_L2[px]+neuronparts[p].R_L1[px+1];
    			else{ 
      				R[x+1] = neuronparts[p].R_L2[px]+neuronparts[p+1].R_L1[0];
      				p++; px=-1;
    			}
		}
		if (px>=neuronparts[p].comp) {p++; px=0;}
\endcode

\Note 
Note here that membrane capacitance is normalized by the time interval.
\code
     		c_m[x] = neuronparts[p].c_m/dt;
  	}
\endcode 
*/

/*!
\warning
For control purposes, the multi-compartmental model can be reduced to a single compartment one. In this case, pay attention to the sites that should be excluded, and that would lead to errors like segmentation faults otherwise. This type of code pieces are marked as "see Warning" in the boxes.  

\par Influx and Outflux
Here the diffusing conductances to the previous (g_minus) and next (g_plus) compartments are calculated 

\code
  	if(comps>1){     //(see Warning)
	for(int p=0, px=1, x=1; x<comps-1; x++, px++)
  	{  
    		if (px>=neuronparts[p].comp) {p++; px=0;}
    	
		g_plus[x] = 1./R[x+1]/neuronparts[p].area[px];
    		g_minus[x] = 1./R[x]/neuronparts[p].area[px];
  	}
\endcode

\par
The initial compartment shouldn't have backward diffusion and the final compartment shouldn't have the forward diffusion

\code
	 g_minus[0] = 0.;
	 g_plus[0] = 1./R[1]/neuronparts[0].area[0];
	 g_minus[comps-1] = 1./R[comps-2]/neuronparts[parts-1].area[neuronparts[parts-1].comp-1];
 	g_plus[comps-1] = 0.;}
\endcode

\par
Then the abbrevations for g_plus and g_minus are given. Note the minus sign. Also V_{mem} is initiated.

\code
for(int x=0; x<comps; x++)
	{
		alpha[x] = -g_minus[x];
    		gamma[x] = -g_plus[x];
		
		pot[x] = V_init;
    		membrane[x].init(V_init);
    		dudt_last[x] = 0.;
     	}

\endcode

*/
}

neuron::~neuron()
{
  delete[] g_plus;
  delete[] g_minus;
  delete[] i_e;
  delete[] c_m;
  delete[] alpha;
  delete[] beta;
  delete[] gamma;
  delete[] delta;
  delete[] pot;
  delete[] membrane;
  delete[] dudt;
  delete[] dudt_last;
  delete[] latency;
  delete[] lengths;
  delete[] i_Na;
  delete[] i_lateral;
}
 
inline void neuron::diffusion_implicit()
{/*! \fn inline void diffusion_implicit 
\brief This is the function which calculates the in
Notice that if comps=1, then the loops wont work. So there is no need of commenting out.
 */
  	for(int x=1; x<comps; x++)
  	{			
		beta[x] -= gamma[x-1] * alpha[x]/beta[x-1];
    		delta[x] -= delta[x-1] * alpha[x]/beta[x-1];
		i_lateral[x] = alpha[x]*(pot[x-1]-pot[x]);
  	}

  i_lateral[0] = 0;
  pot[comps-1] = delta[comps-1]/beta[comps-1];

	for(int x=comps-2; x>=0; x--)
  	{
    		pot[x] =(delta[x] - gamma[x]*pot[x+1])/beta[x];
	}
			

}



void neuron::timestep()
{/*!\brief This is the Core. For each time step, the results of probability calculations are called from WB.cpp, multiplied by the g_bar values and the reversal potential, summed and intermingled into the equations. The equations for diffusion estimation is called at the end to obtain the \DeltaV value (pot[x]). This is repeated for each and every compartment, and a vector "pot" of size "#_{comp}" provided obtained. This timestep() function is called in the main function in a loop that repeats for the given simulation time.  
  */
 /*! If any current is injected at a site, it is normalized by the area of the site (compartment). If not, it is set to 0. 
\code
for(int p=0, px=0, x=0; x<comps; x++, px++)
{
if (px>=neuronparts[p].comp) {p++; px=0;}    
 if ((neuronparts[p].comp_inj == px)&&(neuronparts[p].I_inj != 0.))
i_e[x] = neuronparts[p].I_inj / neuronparts[p].area[px];
else
i_e[x] =0.;
}
\endcode
*/


/*! \par Core loop for DeltaV estimation. 
The loop repeats for every compartment and stores the values. The gating probabilities are called from WB.cpp and multiplied by the \bar{g} values of the original WB model. All of the part here is same as original. Until cooperative part...
\code
  	for(int px=0, p=0, x=0; x < comps; x++, px++)
  	{
  			if (px>=neuronparts[p].comp) {p++; px=0;}

    			if (neuronparts[p].gi[0]>0.) 
				g[0] = neuronparts[p].gi[0];  else g[0] = 0.;
    			if (neuronparts[p].gi[1]>0.) 
				g[1] = neuronparts[p].gi[1]*membrane[x].calculate_Na(pot[x]-neuronparts[p].Vi[0], dt);   else g[1] = 0.;
   			if (neuronparts[p].gi[2]>0.) 
				g[2] = neuronparts[p].gi[2]*membrane[x].calculate_K(pot[x]-neuronparts[p].Vi[1], dt);    else g[2] = 0.;
\endcode

\par
The same procedure repeated for cooperative gating. 
\note 
If the percentage of cooperativity is given as 0 in the main loop (see multicomp.cpp), then this part will not be active and be skipped. Therefore setting pr = 0 naturally means the classical model, which should behave same as the original model.
  \code
  			if (neuronparts[p].pr>0.)
			{
				gj = neuronparts[p].gi[1]*membrane[x].calculate_Naj(pot[x]-neuronparts[p].Vi[0], dt);
			}
\endcode

\par
The currents and the conductances are summed seperately and combined in the equation below: 
\code
  		sum_L =sum_K =sum_Na=0;
  		sum_L = g[0]*neuronparts[p].Ei[0];
  		sum_K = g[2]*neuronparts[p].Ei[2];
  		sum_Na = (((1-neuronparts[p].pr)*g[1]) + (neuronparts[p].pr*gj))*neuronparts[p].Ei[1];
  		i_Na[x] = sum_Na;
  		sum_gi = sum_giEi = 0;
 
		for(int i=0; i<3; i++) 
     			{
      				sum_gi += g[i];
    			}

  		sum_gi +=gj;
  		sum_giEi = sum_Na+sum_K+sum_L;
  
  		beta[x] = c_m[x] + g_minus[x] + g_plus[x] + sum_gi;
  		delta[x] = pot[x]*c_m[x] + i_e[x] + sum_giEi;

   		dudt[x] = pot[x];
\endcode
\par
After this point, the loop ends and the pot[x] vector is sent to the  diffusion_implicit.
\code
  	}
 
   diffusion_implicit();
}
\endcode
*/

// THE ACTUAL CODE STARTS HERE
  	for(int p=0, px=0, x=0; x<comps; x++, px++)
  	{
    		if (px>=neuronparts[p].comp) {p++; px=0;}    
   
    		if ((neuronparts[p].comp_inj == px)&&(neuronparts[p].I_inj != 0.))
      			i_e[x] = neuronparts[p].I_inj / neuronparts[p].area[px];
    		else
      			i_e[x] = 0;
	}
	double gj = 0; /**< NEW. Defines the cooperative sodium conductance.*/

  	for(int px=0, p=0, x=0; x<comps; x++, px++)
  	{
		if (px>=neuronparts[p].comp) {p++; px=0;}
           
    			if (neuronparts[p].gi[0]>0.) 
				g[0] = neuronparts[p].gi[0];  else g[0] = 0.;
    			if (neuronparts[p].gi[1]>0.) 
				g[1] = neuronparts[p].gi[1]*membrane[x].calculate_Na(pot[x]-neuronparts[p].Vi[0], dt);   else g[1] = 0.;
   			if (neuronparts[p].gi[2]>0.) 
				g[2] = neuronparts[p].gi[2]*membrane[x].calculate_K(pot[x]-neuronparts[p].Vi[1], dt);    else g[2] = 0.;
    			//if (neuronparts[p].pr>0.)
			//{
			  //if (neuronparts[p].gi[1]>0.)
			  //gj = neuronparts[p].gi[1]*membrane[x].calculate_Na(pot[x]-neuronparts[p].Vi[0], dt);   else gj = 0.;
 			gj = neuronparts[p].gi[1]*membrane[x].calculate_Naj(pot[x]-neuronparts[p].Vi[0], dt);
			//}
		
// 		 cout<<"check x"<<endl; 
  		sum_L =sum_K =sum_Na=0;
  		sum_L = g[0]*neuronparts[p].Ei[0];
  		sum_K = g[2]*neuronparts[p].Ei[2];
  	         sum_Na = (((1-(neuronparts[p].pr))*g[1]) + ((neuronparts[p].pr)*gj))*neuronparts[p].Ei[1];
  		
	     
		i_Na[x] = sum_Na;
 //		i_Na[x] = g[1]*(pot[x]-neuronparts[p].Ei[1]);
  		sum_gi = sum_giEi = 0;
 
   // 			cout << "Ie: \t\t" << i_e[16] << "\t" << i_e[17]<< endl;

			for(int i=0; i<3; i++) 
     			{
      				sum_gi += g[i];
      //			sum_giEi += g[i]*neuronparts[p].Ei[i];
    			}

  		sum_gi +=gj;
  		sum_giEi = sum_Na+sum_K+sum_L;
                                                                                                                                                                                                                                         
  		beta[x] = c_m[x] + g_minus[x] + g_plus[x] + sum_gi;
		delta[x] = pot[x]*c_m[x] + i_e[x] + sum_giEi;
		dudt[x] = pot[x];
  	}
 
   diffusion_implicit();
}

int neuron::compartments()
{
  return comps;
}

bool neuron::check_latency(double dudt_ap, double t)
{
  
  for(int x=0; x<comps; x++)
    dudt[x] = (pot[x]-dudt[x])/dt;         //notwendig, dass vorher in timestep dudt=pot vor gaussreduktion
    //                          /dt hier weg -> schneller ?
  
  //find the AP initiation site
  bool flag = 0;
  dudt_max=0.;
  if(ap[0])
    flag = 1;
  
  for(int x=0; x<comps; x++)
    if ((dudt_max < dudt[x])&& (flag == 0))
    {
      dudt_max = dudt[x];
      max = x;
     }
    
     // ap onset ?
     
  if ((dudt[max]>=dudt_ap)&&(dudt_last[max]<dudt_ap))
  {
    for (int i=0; i<4; i++)
      if (!ap[i])
      {
        ap[i] = t;
        cout << i+1 << ". ap: " << ap[i] << " ms" << endl;
        
      if (ap[3])
      {cout<<"max:"<< max <<endl;
           freq = 1000./(ap[2]-ap[1]);
      }   
        
        break;
      }
  }
     
      // third ap -> calc onset latencies
    if (ap[2]&&!ap[3])
    
      for(int x=0; x<comps; x++)
       if ((dudt[x]>=dudt_ap)&&(dudt_last[x]<dudt_ap))
        {
          latency[x] = t-ap[2];
          comp_lat++;
        }
        
        
        
    for(int x=0; x<comps; x++)
      dudt_last[x] = dudt[x];
  
    if (comp_lat==comps) return true;
    else return false;
}


