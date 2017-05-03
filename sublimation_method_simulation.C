#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;
const long double pi = 3.1415926535897932384626433832795;

//--------------------------------- Functions Begin ---------------------------------------------------------
long double ay(long double t, long double x, long double y, long double F, long double phi0, long double m, long double w, long double GM, long double m_nuc)
{
  long double r=sqrt(pow(x,2)+pow(y,2));
  if(m > m_nuc) return -GM*y/pow(r,3)+F*sin(w*t+phi0)/m;
  else return -GM*y/pow(r,3) ;
}

long double ax(long double t, long double x, long double y, long double F, long double phi0, long double m, long double w,long double GM, long double m_nuc)
{
  long double r=sqrt(pow(x,2)+pow(y,2));
  if(m > m_nuc) return -GM*x/pow(r,3)+F*cos(w*t+phi0)/m;
  else return -GM*x/pow(r,3) ;
}

bool ChangeRev(long double y_next,long double y_prev)
{
   if (y_prev < 0 && y_next >= 0)
      return true;
  return false;
}
//------------------------------ Functions End ------------------------------------------------------
// ----------------------------- Main Function Begin ----------------------------------------------
void Kepler(int nrev)			   // nrev - number of revolutions
{
// ------------------ initial parameters -------------

  long double nw_rev = 0;		   // number of revolutions of asteroid around its axis per day
  long double w  = nw_rev * 2 * pi / 86400; 	   // angular velocity of asteroid's rotation around its axis
  long double t_  = 0;
  long double dt  = 100 ; 		   // integration step
  long double dt_print = 86400;		   // printing step
  long double phi0_grad = 180;
  long double phi0  = phi0_grad * pi / 180;      // initial angle of vector F0 relative to X axis
  long double r_ast  = 50;  	  	   // asteroid's radius
  long double ro_ast = 1000;		   // asteroid's density
  long double k_F = 47.5 ;		   // Thrust factor
  long double W0 = k_F * 1000 ;	   	   // the heat generation power of nuclear power unit
  long double W;
  long double F;
  long double m_nuc = 0.0001064*W0; 		   // mass of nuclear power unit
  long double m = 4 * pi * pow(r_ast,3) * ro_ast / 3 + m_nuc;   // asteroid's mass
  long double t_year = 284010624 ;
  long double T12 = 1* t_year  ;		   // half-life
  long double Lambda  =  335000 ;	   // enthalpy of fusion - ice
  long double L_boil  = 2256000 ; 	   // specific heat of evaporation - water
  long double C_ice   = 4190000 ;	   // specific heat - ice
  long double C_steam = 2000000 ;	   // specific heat - steam
  long double delta_T = 0 ;
  long double T0 = 173 ;		   //  K :  0 < T0 < 273
  long double R  = 8.3144598;	 	   //  gas constant
  long double ms = -W0/(C_ice * (273-T0)+Lambda+L_boil+C_steam * delta_T);		   // mass loss rate
  long double C  = sqrt(1000 * 3 * 1.38 * (T0+delta_T) / 2.7) ;          //  outflow velocity
  long double F0 = abs(ms) * (C+R*(T0+delta_T) / C / 0.018);         // jet thrust
  long double x0 = 149597870700;
  long double x  = x0 ;   // initial x
  long double y = 0;	   // initial y
  long double alpha=0, t = 0;
  long double Vy = 39617.8298017968, Vx = 0;
  long double n  = 0;			   // number of revolutions around the Sun
  long double G  = 6.67545E-11;	 	   // gravitational constant
  long double M  = 1.9885E30;	 	   // mass of Sun
  long double GM = G*M;
  long double r  = sqrt(pow(x,2)+pow(y,2));    // distance from Sun
  long double V  = Vy;
  long double Impulse = 0;	    	   // Imparted impulse

//-----------coifficients Runge-Kutta -----------------
  long double k_vy[5] = {0};  	               	   // k, l - coefficients for Vx, Vy
  long double k_vx[5] = {0};
  long double k_y[5]  = {0}; 	               	   // s,  j - coeficeints for X, Y
  long double k_x[5]  = {0};
  long double u[5]    = {0,0,1./2,1./2,1};
  long double Vx_temp[5] = {0};               	   //  Vx_temp =Vx(1), Vx(2) Vx(3), Vx(4)
  long double Vy_temp[5] = {0};
  long double t_temp[5]  = {0};
  long double x_temp[5]  = {0};
  long double y_prev     =  0;
  long double y_temp[5]  = {0};
  long double m_temp[5]  = {0};
  long double ms_temp[5] = {0};
  long double W_temp[5]  = {0};
  long double F_temp[5]  = {0};
  long double k_Imp[5]   = {0};

  //-----------File I/O-------------------------------
  int n_lines=0;
  ofstream myfile;
  myfile.open ("Kepler.txt");
  myfile << "W0\tms\tF0\tphi0_grad\tnw_rev\tr_ast\tm_nuc\n";
  myfile << setprecision(15) << W0 <<"\t";
  myfile << ms <<"\t";
  myfile << F0 <<"\t";
  myfile << phi0_grad <<"\t";
  myfile << nw_rev <<"\t";
  myfile << r_ast <<"\t";
  myfile << m_nuc <<"\n";
  myfile << "t\tx\ty\tr\talpha\tVx\tVy\tms\tm\tW\tF\n";

  myfile << setprecision(15) << t  <<"\t";
  myfile << x  <<"\t";
  myfile << y  <<"\t";
  myfile << r  <<"\t";
  myfile << alpha <<"\t";
  myfile << Vx <<"\t";
  myfile << Vy <<"\t";
  myfile << ms  <<"\t";
  myfile << m  <<"\t";
  myfile << W0  <<"\t";
  myfile << F0 <<endl;

// ----------------starting loop,--------------------------------------
  while (n < nrev)     // n is number of revolutions around the Sun
   {
//------------------calculating coificinets of  Runge Kutta -----------
     for(int i=1; i<=4; i++)
       {
        Vy_temp[i]  =  Vy+u[i] * k_vy[i-1];
        Vx_temp[i]  =  Vx+u[i] * k_vx[i-1];
        t_temp[i]   =  t + u[i] * dt;
        W_temp[i]   =  W0 * pow(2,-t_temp[i] / T12) ;
        ms_temp[i]  =  -W_temp[i] / (C_ice * (273-T0) + Lambda + L_boil + C_steam * delta_T);
F_temp[i]   =  abs(ms_temp[i]) * (C + R * (T0+delta_T) / C / 0.018);
k_Imp[i]    =  F_temp[i] * dt ;
m_temp[i]   =  m + u[i] * ms_temp[i-1] * dt;
if(m_temp[i] < m_nuc) m_temp[i]=m_nuc ;
k_y[i]      =  Vy_temp[i] * dt;
        y_temp[i]   =  y + u[i] * k_y[i-1];
        k_x[i]      =  Vx_temp[i] * dt;
        x_temp[i]   =  x+u[i] * k_x[i-1]; 	//  x_temp[i] x of small cicle
        k_vy[i]     =  dt * ay(t_temp[i], x_temp[i], y_temp[i], F_temp[i], phi0, m_temp[i], w, GM, m_nuc);
        k_vx[i]     =  dt * ax(t_temp[i], x_temp[i], y_temp[i], F_temp[i], phi0, m_temp[i], w, GM, m_nuc);
       }

//-------------- calculating parametres---------------------
        t      =  t + dt;
        W      =  W0 * pow(2,-t / T12) ;
        ms     =  -W/(C_ice * (273 - T0) + Lambda + L_boil + C_steam*delta_T);
        F      =  abs(ms) * (C + R * (T0+delta_T) / C / 0.018);
        m      =  m + ms * dt;
        if(m < m_nuc) m = m_nuc ;
        y_prev =  y;
        x      =  x + ( k_x[1] + 2 * k_x[2] + 2 * k_x[3] + k_x[4] ) / 6;
        y      =  y + ( k_y[1] + 2 * (k_y[2] + k_y[3] ) + k_y[4] ) / 6;
        r      =  sqrt( pow(x,2) + pow(y,2) );
        Vy     =  Vy + ( k_vy[1] + 2 * ( k_vy[2] + k_vy[3] ) + k_vy[4] ) / 6;
        Vx     =  Vx + ( k_vx[1] + 2 * ( k_vx[2] + k_vx[3] ) + k_vx[4] ) / 6;
        alpha  =  atan2 (y,x) * 180 / pi;
        V      =  sqrt( pow(Vx,2) + pow(Vy,2) );
        Impulse = Impulse + (k_Imp[1] + 2 * k_Imp[2] + 2 * k_Imp[3] + k_Imp[4]) / 6;
        n_lines ++ ;
//-----------------condition for stopping while cicle--------
        if ( ChangeRev(y,y_prev) )
               n++;
//---------------Writing into file --------------------------------------------

if((t-t_) > dt_print)
{
t_ =  t;
myfile << setprecision(15) << t  <<"\t";
myfile << x  <<"\t";
myfile << y  <<"\t";
myfile << r  <<"\t";
myfile << alpha <<"\t";
myfile << Vx <<"\t";
myfile << Vy <<"\t";
myfile << ms  <<"\t";
myfile << m  <<"\t";
myfile << W  <<"\t";
myfile << F <<endl;
}
}
 myfile << setprecision(15) << t  <<"\t";
myfile << x  <<"\t";
myfile << y  <<"\t";
myfile << r  <<"\t";
myfile << alpha <<"\t";
myfile << Vx <<"\t";
myfile << Vy <<"\t";
myfile << ms  <<"\t";
myfile << m  <<"\t";
myfile << W  <<"\t";
myfile << F <<endl;
long double ddx = (x - x0) / 1000 ;
long double ddy = (t - t_year) * 30 ;  // 284010624 � 1 year
long double ddr = sqrt(pow(ddx,2) + pow(ddy,2));
cout<< "r_ast="<<r_ast<<endl;
cout<< "W= "<<W0/1.E+06<<endl;
cout<<"F= "<<F0<<endl;
cout<<"ddx = "<<ddx<<" km"<<endl ;
cout<<"ddy = "<<ddy<<" km"<<endl ;
cout<<"ddr = "<<ddr<<" km"<<endl ;
cout<<"t = "<<t<<endl;
cout<<"Impulse  = "<<Impulse<<endl;
  //------------------ploting graph-----------------------------------------
    myfile.close();

//------------------------main Function End------------------
 }

 int main()
 {
	 Kepler(1);
	 return 0;
 }
