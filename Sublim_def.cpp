#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
//#include <TFile.h>
//#include <TCanvas.h>
//#include <TF1.h>
//#include <TGraph.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
//#include <TMath.h>
#include <vector>
#include <iomanip>

using namespace std;

const double pi = 3.1415926535897932384626433832795; 


//---------------------------------Functions Begining---------------------------------------------------------
double ay(double t, double x, double y, double z, double F, double phi0, double m, double GM, double m_nuc, double zetta) 
{
	double r = sqrt(pow(x, 2) + pow(y, 2)+pow(z, 2));
	if (m > m_nuc) return -GM*y / pow(r, 3) + F*sin(phi0)*cos(zetta) / m;
	else return -GM*y / pow(r, 3);
	
	
}

double ax(double t, double x, double y, double z, double F, double phi0, double m, double GM, double m_nuc, double zetta)
{
	double r = sqrt(pow(x, 2) + pow(y, 2)+pow(z, 2));
	if (m > m_nuc) return -GM*x /pow(r, 3) + F*cos(phi0)*cos(zetta) / m;
	else return -GM*x / pow(r, 3);
													 
	
}
double az(double t, double x, double y, double z, double F, double phi0, double m, double GM, double m_nuc, double zetta)
{
	double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	if (m > m_nuc) return -GM*z / pow(r, 3) + F*sin(zetta) / m;
	else return -GM*z / pow(r, 3);
}



//------------------------------ Functions End ------------------------------------------------------


// -----------------------------  Main Function Start --------------------------------------------------

void Sublim_def(int nrev)					// nrev - number of revolutions
{
	// ----------------- Core inputs and calculation of pebble numbers for one reactor block ----------------
	double R_c =205;  // core radius in cm
	double core_levels = 5; //number of pebble layers per block
	cout << "core levels = " << core_levels << endl;
	int count_in = 0; //variable for inner pebble number calculations
	int N_cell = 0; //cell number
	int count_near = 0; //variable for lateral pebble number calculations
	double N_level_spring = 0;
	double Np_diam = 2*int (R_c /9) + 1;
	int c_c = 0; // variable for cell number calculations


	for (int i = 9 * int(R_c / 9); i >= 0; i = i - 9)
	{
		c_c = count_in;

		for (int j = -9 * int(sqrt(R_c*R_c - i*i) / 9); j <= 0; j = j + 9)
		{
			if ((((j-9)*(j-9) + i*i )> R_c*R_c)||(((i+9)*(i+9) + j*j )> R_c*R_c))
			{
				count_in++;
				count_near++;

			}
			else
			{
				count_in++;
			}




		}
		if ((count_in - c_c) == 1)
		{
			N_level_spring++;
		}
		if ((i - 9) >= 0)
		{
			N_cell = N_cell + (count_in - c_c - 1);
		}

	}
	cout << "ccvac  = " << N_level_spring << endl; // this will show, if there are any pebbles that are connected with their layer by a single spring

	double Np_near_level =(4 * count_near )- 4; //-4 because 4 of the pebbles are calculated 2 times
	cout << "Np_outer_layer = " << Np_near_level << endl; // number of outer pebbles for 1 layer
	cout << "Np_diam = " << Np_diam << endl; // number of pebbles on diameter


	double Np_base = (count_in - Np_diam) * 4 + (2 * Np_diam)-1;   //number for inner pebbles on each layer
	cout << "Np base = " << Np_base << endl;
	double Np = Np_base * core_levels; // number of inner pebbles
	double Np_lat = (Np_near_level*(core_levels - 2)); //number of lateral pebbles
	N_level_spring = N_level_spring * 4;
	N_cell = 4 * N_cell;
	N_level_spring = N_level_spring + (Np_base + (N_cell + 1) - 2);     //Eyler theorem
	double N_spring = (N_level_spring * core_levels) + (Np_base * (core_levels - 1)); //final number of springs

	cout <<"N_spring = "<< N_spring << endl;

						

									   // ------------------ Initial Parametres --------------------------------------------------------
	double t_ = 0;
	double dt = 10;		 	    // integration step
	double dt_print = 86400;		// printing step
	double incline_deg = 0;         // NEO-s planes inclination from ecliptic plane
	double incline = incline_deg * pi / 180; //the same but in radians
	double zetta_grad =0;              //F_jet angle to orbit plane (-90 to 90)
	double zetta = zetta_grad*pi / 180;    //the same but radians
	double phi0_grad =90;			// initial angle of vector F0 relative to X axis (deg)
	double phi0 = phi0_grad*pi / 180;  // the same but radians
	double r_ast = 500;             // asteroid's radius (m)
	cout << "r_ast =" << r_ast << endl;
	double ro_ast = 1000;		    // asteroid's density (kg/m^3)
	double F;						// jet thrust
	double P;						// pressure of outflowing water vapor
	double k_F = 1;					// thrust coefficient        
	double t_year = 365.24 * 86400;	// (s)
	double Tp_PBMR = 900;			// pebble surface temperature (degree C)
	double Tp = Tp_PBMR + 273;      // the same but K
	double Rp = 0.03;				// pebble radius
	double a = 0.09;				// size of one lattice unit
	double k_p = pow(a / (2 * Rp), 2);	// neutron flux scale factor - uranium enrichment factor
	double m_u_1 = 0.009;           // mass of uranium in one pebble (kg)
	double dp_gap = a - 2 * Rp;		// gap between neighboring pebbles 
	double V_1u = pow(a, 3);			// volume of one lattice unit
	double S_u_1 = a*a;              // area of one cell
	double N_block = 8;				//number of reactor <<slices>> *********************************************************************************************************************************//
	cout << "N block = " << N_block << endl;
	double Sp_1 = 4 * pi*pow(Rp, 2);	// surface area of one pebble
	double Sp = Np*Sp_1;	        // total surface area of pebbles
	double Vp_1 = (4.0 / 3.0)*pi*pow(Rp, 3);// volume of one pebble
	double d_boron = 0.0002;		// Thickness of boron carbid layer on a pebble
	double ro_BC = 2520;			// density of boron carbide kg/m^3
	double m_boron_1 = ro_BC*Sp_1*d_boron; // mass of absorbong layer on a pebble
	double m_pebble_1 = Vp_1 * 2000 + k_p*m_u_1; // 0.246 kg;	 mass of one pebble
	double m_pebbles = Np*(m_pebble_1 + m_boron_1);  // total mass of pebbles
	double d_shell = 0.01;		    // thickness of reactor wall    cm
	double w_p_gap = (Np_diam - 1)*0.009/2;  // the gap between the lateral wall and core taking into account expansion of core at 900 deg C
	double h_core = (core_levels - 1)*a;  // height of the reactor core
	double r_core = R_c / 100;  // reactor core radius   m
	double h_reactor = h_core + 2 * (d_shell + Rp + w_p_gap);  // reactor core height
	double r_reactor = r_core + (d_shell + Rp + w_p_gap);  // radius of the reactor   
	double m_spring_1 = 0.037617;			 // mass of one spring
	double m_springs = m_spring_1 * N_spring; //total mass of springs
	cout << " mass of springs = " << m_springs/1000 <<" tone" << endl;				
	double s_shell_h = 2 * pi*(r_reactor - d_shell / 2)*(h_reactor - 2 * d_shell); // area of the reactor's lateral surface
	double s_shell_r = 2 * pi*(pow(r_reactor, 2));   // area of the reactor's bases
	double s_shell = s_shell_h + s_shell_r;		 // total area area of the reactor's surface
	double m_shell = 1600 * d_shell*s_shell;       // mass of the reactor wall; 1600 - graphite-epoxy density
	double m_nuc = N_block*(m_pebbles + m_springs + m_shell);    	// mass of the reactor (kg)
	double k_W = N_block*(Np - 0.5*Np_base - 0.1666*Np_lat) * 400 / 452000; // reactor power scale-up  
	double W0 = k_F*k_W * 1000000;				 // thermal power of the nuclear reactor (MW)
	double W = W0;                               
	double Ts;							         // ice sublimation temperature
	double m0 = 4 * pi*pow(r_ast, 3)*ro_ast / 3 + m_nuc; // asteroid's mass (kg)
	double m = m0;                  // current mass of NEO
	double m_loss = 0;              // the NEO mass loss
	double h_well = 0;              // borehole height
	double Lambda = 335000;	        // enthalpy of fusion - ice  (Si)
	double L_boil = 2256000; 	    // specific heat of evaporation - water (Si)
	double C_ice = 1300;	        // specific heat - ice (Si)
	double C_steam = 1850;	        // specific heat - water vapor at ~200K (Si)
	double delta_T = 412;			// 412 increase of water vapor temperature (K)
	double T0 = 40;			        // temperature of subsurface ice (K) :  0 < T0 < Ts0
	double T;                       // tempreture of gas outflowing from the core
	double R = 8.3144598;	 	    // gas constant
	double x0 = 149597870700;       // radius of Earths orbit
	double x = x0;					// initial x
	double y = 0;                   // initial y
	double z = 0;					// initial z
	double x_e;						// x coordinate of point E
	double y_e;						// y coordinate of point E
	double ddx;						// x component of ddr
	double ddy;						// y component of ddr
	double ddz;						// z component of ddr
	double ddr = 0;					// minimal bypassing distance
	double ddr_;				    // distance from point E at the previous integration step 
	double alpha = 0, alpha_deg = 0;// current angle of NEO from x axis
	double t = 0;					// time passed from the start of operation
	double t_c = 0;					// the time of reaching the NEO's center
	double dVy = 0;                 // variable to calculate mass of rocket fuel in case of rocket propulsion method. 0.207 ;
	double Vp = 39617.8298017968;   // 9 year orbit - periheluim velocity  m/s
	double Vy = Vp * cos(incline) + dVy; // y component of Vp
	double Vx = 0;                  // x component of Vx
	double Vz = Vp * sin(incline);  // z component of Vz
	double n = 0;			 	    // number of revolutions around the Sun   
	double G = 6.67545E-11;	 	    // gravitational constant (Si)
	double M = 1.9885E30;	 	    // mass of Sun  (kg)
	double GM = G * M;
	double r = sqrt(pow(x, 2) + pow(y, 2)+ pow(z, 2));  //  distance from Sun (m)
	double V;                       // spacial velocity of NEO
	double m_min = 14137167;		// minimum mass - can burn out in the atmosphere: corresponds to 15 m radius

									// ---------------------------- solving equation to find Ts ---------------->--------
	double a_Ts = 3.41E+12;         // coefficient of water ice sublimation curve
	double b_Ts = 6130;             // coefficient of water ice sublimation curve
	double r_jet = 1.5*r_reactor;   // borehole/jet radius
	double S_jet = pi*pow(r_jet, 2);// area of borehole
	double Tmax;                    // variable for solving the equation for Ts 
	double Tmin;                    // variable for solving the equation for Ts 
	double Ts0;						// variable for solving the equation for Ts 
	double diff;					// variable for solving the equation for Ts 
	double ms;                      // NEO mass loss rate
	double C;                       // jet outflow velocity         
	Tmax = 273;
	Tmin = 40;				        // initial subsurface temperature
	Ts0 = Tmax;
	Ts = (Tmax + Tmin) / 2;
	while (abs(Ts - Ts0) > 1)
	{

		C = sqrt(1000 * 1.38*(Ts + delta_T) / 2.7);
		ms = -W / (C_ice*(Ts - T0) + Lambda + L_boil + C_steam*delta_T);
		diff = a_Ts*exp(-b_Ts / Ts) - abs(ms)*R*(Ts + delta_T) / (0.018*S_jet*C);
		Ts0 = Ts;
		if (diff < 0)
		{
			Tmin = Ts;
			Ts = (Ts + Tmax) / 2;
		}
		else
		{
			Tmax = Ts;
			Ts = (Ts + Tmin) / 2;
		}
	}
	// ----------------------------- solving equation to find Ts ------------------------
	T = Ts + delta_T;
	double F0 = abs(ms)*(C + R*T / C / 0.018);
	F = F0;
	P = abs(ms)*R*T / (0.018*S_jet*C);

	cout << "Np = " << Np << endl;
	cout << "W = " << k_W << " MWt" << endl;
	cout << "phi0_grad = " << phi0_grad << " degree" << endl;
	cout << "zetta_grad = " << zetta_grad << " degree" << endl;
	cout << "r_reactor = " << r_reactor << " m" << endl;
	cout << "d_reactor = " << r_reactor * 2 << " m" << endl;
	cout << "h_reactor = " << h_reactor << " m" << endl;
	cout << "m_pebble_1 = " << m_pebble_1 << " kg" << endl;
	cout << "m_nuc = " << m_nuc / 1000 << " tone" << endl;
	cout << "Tsub = " << Ts << " K" << endl;
	cout << "ms = " << ms << " kg/s" << endl;
	cout << "P = " << P << " Pa" << endl;
	cout << "F = " << F << " N" << endl;


	//----------- Runge-Kutta ----------------- //

	double k_vy[5] = { 0 };
	double k_vx[5] = { 0 };
	double k_vz[5] = { 0 }; 

	double k_y[5] = { 0 };
	double k_x[5] = { 0 };
	double k_z[5] = { 0 };
	double u[5] = { 0, 0, 1./2, 1./2, 1 };
	double Vx_temp[5] = { 0 };
	double Vy_temp[5] = { 0 };
	double Vz_temp[5] = { 0 };
	double t_temp[5] = { 0 };
	double y_prev = 0; 
	double x_temp[5] = { 0 };
	double y_temp[5] = { 0 };
	double z_temp[5] = { 0 }; 
	double m_temp[5] = { 0 };
	double ms_temp = 0;
	double W_temp =  0 ;
	double F_temp =  0;
	double k_Imp = 0 ;

	//----------- File I/O -------------------------------

	
	ofstream myfile;

	{
		char buffer[100];
		sprintf(buffer, "Sublim_Def _%d_%d.txt", int(phi0_grad), int(r_ast));
		myfile.open(buffer);
	}

	myfile << "Np\tr_nuc\th_nuc\tm_nuc\tW0\n";

	myfile << setprecision(15) << Np << "\t";
	myfile << r_reactor << "\t";
	myfile << h_reactor << "\t";
	myfile << m_nuc << "\t";
	myfile << W0 << endl;
	myfile << "t\tx\ty\tz\tr\tddr\talpha\tVx\tVy\tVz\tms\tm\tW\tF\n";
	myfile << setprecision(15) << t << "\t";
	myfile << x << "\t";
	myfile << y << "\t";
	myfile << z << "\t";
	myfile << r << "\t";
	myfile << ddr << "\t";
	myfile << alpha_deg << "\t";
	myfile << Vx << "\t";
	myfile << Vy << "\t";
	myfile << Vz << "\t";
	myfile << ms << "\t";
	myfile << m << "\t";
	myfile << W << "\t";
	myfile << F << endl;

 // ----------------loop start--------------------------------------
	while (n < 1) // n is number of revolutions around the Sun
	{
		W_temp = k_F * W0; 
		ms_temp = -W_temp / (C_ice*(Ts - T0) + Lambda + L_boil + C_steam*delta_T);
		F_temp = abs(ms_temp)*(C + R*(Ts + delta_T) / C / 0.018);
		k_Imp = F_temp * dt;
		//------------------ using Runge Kutta formulas -----------
		for (int i = 1; i <= 4; i++)
		{

			Vy_temp[i] = Vy + u[i] * k_vy[i - 1];
			Vx_temp[i] = Vx + u[i] * k_vx[i - 1];
			Vz_temp[i] = Vz + u[i] * k_vz[i - 1]; 
			t_temp[i] = t + u[i] * dt;
			m_temp[i] = m + u[i] * ms_temp * dt;
			if (m_temp[i] < m_nuc) m_temp[i] = m_nuc;
			k_y[i] = Vy_temp[i] * dt;
			y_temp[i] = y + u[i] * k_y[i - 1];
			k_x[i] = Vx_temp[i] * dt;
			x_temp[i] = x + u[i] * k_x[i - 1];
			k_z[i] = Vz_temp[i] * dt; 
			z_temp[i] = z + u[i] * k_z[i - 1]; 
			k_vy[i] = dt * ay(t_temp[i], x_temp[i], y_temp[i], z_temp[i], F_temp, phi0, m_temp[i], GM, m_nuc, zetta);
			k_vx[i] = dt * ax(t_temp[i], x_temp[i], y_temp[i], z_temp[i], F_temp, phi0, m_temp[i], GM, m_nuc, zetta);
			k_vz[i] = dt * az(t_temp[i], x_temp[i], y_temp[i], z_temp[i], F_temp, phi0, m_temp[i], GM, m_nuc, zetta);
		}
		//-------------- calculating parametres --------------------
		t = t + dt;
		m_loss = m0 - m;
		h_well = m_loss / N_block / S_jet / ro_ast;
		if (h_well >= r_ast)
		{
			P = 0;
			F = 0;
			ms = 0;
			W = 0;
			k_F = 0; 
		}
		else t_c = t; 
		m = m + ms * dt;
		if (m <= m_nuc) m = m_nuc; // <= reactor stop criteria
		x = x + (k_x[1] + 2 * (k_x[2] + k_x[3]) + k_x[4]) / 6;
		y = y + (k_y[1] + 2 * (k_y[2] + k_y[3]) + k_y[4]) / 6;
		z = z + (k_z[1] + 2 * (k_z[2] + k_z[3]) + k_z[4]) / 6; 
		r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		Vy = Vy + (k_vy[1] + 2 * (k_vy[2] + k_vy[3]) + k_vy[4]) / 6;
		Vx = Vx + (k_vx[1] + 2 * (k_vx[2] + k_vx[3]) + k_vx[4]) / 6;
		Vz = Vz + (k_vz[1] + 2 * (k_vz[2] + k_vz[3]) + k_vz[4]) / 6;
		V = sqrt(pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2));
		alpha = atan2(y, x); 
		if (y < 0) alpha = 2 * pi + alpha;
		alpha_deg = alpha * 180 / pi;
		if (m < m_min) cout << "Minimum mass reached at alpha=" << alpha << endl;
	
		x_e = x0 * cos(2 * pi*t / t_year);
		y_e = x0 * sin(2 * pi*t / t_year);
		ddx = (x - x_e) / 1000;
		ddy = (y - y_e) / 1000;
		ddz = z / 1000; 
		ddr_ = ddr;
		ddr = sqrt(pow(ddx, 2) + pow(ddy, 2) + pow(ddz, 2));

		//---------------Writing into file --------------------------------------------

		/*if ((t - t_) > dt_print)
		{
			t_ = t;
		    myfile << setprecision(15) << t << "\t";
			myfile << x << "\t";
			myfile << y << "\t";
			myfile << z << "\t";
			myfile << r << "\t";
			myfile << ddr << "\t";
			myfile << alpha_deg << "\t";
			myfile << Vx << "\t";
			myfile << Vy << "\t";
			myfile << Vz << "\t";
			myfile << ms << "\t";
			myfile << m << "\t";
			myfile << W << "\t";
			myfile << F << endl;  
		}
		*/
		if ((t > 8.5*t_year) && (ddr > ddr_))  n++;  // t > 8.5*t_year  -  the rendezvous should be not earlear than the last half year 
	}
	myfile << setprecision(15) << t << "\t";
	myfile << x << "\t";
	myfile << y << "\t";
	myfile << r << "\t";
	myfile << ddr << "\t";
	myfile << alpha_deg << "\t";
	myfile << Vx << "\t";
	myfile << Vy << "\t";
	myfile << ms << "\t";
	myfile << m << "\t";
	myfile << W << "\t";
	myfile << F << endl;

	cout << "x = " << x/1000 << " km" << endl;
	cout << "y = " << y/1000 << " km" << endl;
	cout << "z = " << z / 1000 << " km" << endl;
	cout << "ddx = " << ddx << " km" << endl;
	cout << "ddy = " << ddy << " km" << endl;
	cout << "ddz = " << ddz << " km" << endl;
	cout << "ddr = " << ddr << " km" << endl;
	cout << "ddr_ = " << ddr_ << " km" << endl;
	cout << "alpha_deg = " << alpha_deg << " deg" << endl;
	cout << "t = " << t<< endl;
	cout << "dt = " << (t - 9 * 365.24 * 86400) << " sec" << endl;
	cout << "m0     = " << m0 << endl;
	cout << "m      = " << m << endl;
	m_loss = m0 - m;
	cout << "m_loss = " << m_loss << endl;
	cout << "r_jet = " << r_jet << "\t" << endl;
	cout << "h_well = " << h_well << endl;
	cout << "t_c = " << t_c / 365.24/ 86400 << " year" << endl;

	double m_fuel = m0*dVy / 3000 / 1000;
	cout << "m_fuel = " << m_fuel << "ton" << endl;
	cout << "ax = " << ax(t, x, y, z, F, phi0, m, GM, m_nuc, zetta) << " m/s^2" << endl;
	cout << "ay = " << ay(t, x, y, z, F, phi0, m, GM, m_nuc, zetta) << " m/s^2" << endl;
	cout << "az = " << az(t, x, y, z, F, phi0, m, GM, m_nuc, zetta) << " m/s^2" << endl;


	myfile << "\n\n\n\nphi0_grad\tr_ast\tW0\tdx\tdy\tdr\tF0\tm_nuc\n";

	myfile << phi0_grad << "\t";
	myfile << r_ast << "\t";
	myfile << W0 << "\t";
	myfile << ddx << "\t";
	myfile << ddy << "\t";
	myfile << ddr << "\t";
	myfile << F0 << "\t";
	myfile << m_nuc / 1000 << "\n";
	myfile.close();

	//------------------------main Function End------------------
}

int main()
{
	Sublim_def(0);
	return 0;
}
