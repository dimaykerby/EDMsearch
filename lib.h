#ifndef __lib_h__
#define __lib_h__

#include"random.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
//#include <TVector3.h>
//#include <TMath.h>
//#include <TF1.h>

using namespace std;

double dot_product(vector <double> vector_a, vector <double> vector_b) {
   double product = 0.;
   for (int i = 0; i < 3; i++)
   product = product + vector_a[i] * vector_b[i];
   return product;
}

vector <double> cross_product(vector <double> v_A, vector <double> v_B) {
   vector <double> temp = {0,0,0};
   temp[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   temp[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   temp[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];

   return temp;
}


//this function prepares the random numbers
void set_random_gen( Random &ran_gen ){
    int seed[9];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             ran_gen.SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

//angular distribution of the decay
//p[0]=alpha_psi,p[1]=alpha_1,p[2]=alpha_2, p[3]=eta
//x[0]=n1x,x[1]=n2x,x[2]=n1y, x[3]=n2y, x[4]=n1z, x[5]=n2z
double angular_distribution(vector <double> x, vector <double> p, string s){
  if(s=="direction")
    return 3+p[0]+(p[1]*p[2])*((1-p[0])*(x[0]*x[1]-x[2]*x[3])-(1+3*p[0])*x[4]*x[5]);

  else {
    double n1x = sin(x[0])*cos(x[1]);
    double n1y = sin(x[0])*sin(x[1]);
    double n1z = cos(x[0]);
    double n2x = sin(x[2])*cos(x[3]);
    double n2y = sin(x[2])*sin(x[3]);
    double n2z = cos(x[2]);
    return 3+p[0]+(p[1]*p[2])*((1-p[0])*(n1x*n2x-n1y*n2y)-(1+3*p[0])*n1z*n2z);
  }
}


//this function generates randomly points in the unit circle
vector <double> angle_unif(Random& Genny, string s) {
  vector <double> dir {Genny.Gauss(0,1),Genny.Gauss(0,1),Genny.Gauss(0,1) };
  double norm = pow(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2],0.5);
  dir[0] = dir[0]/norm;
  dir[1] = dir[1]/norm;
  dir[2] = dir[2]/norm;
  if(s=="direction")
    return dir;

  vector <double> angle {0, 0};
  angle[0] = acos(dir[2]);
  if(atan2(dir[1],dir[0])>0)
    angle[1] = atan2(dir[1],dir[0]);
  else angle[1] = atan2(dir[1],dir[0])+2*M_PI;
  return angle;
}

vector <double> accept_reject_indirect(Random& Genny, vector <double> parameters) {
  vector <double> direction1 = angle_unif(Genny,"angles");
  vector <double> direction2 = angle_unif(Genny,"angles");

  vector <double> x = {0,0,0,0};
  for(int i=0; i<2;i++) {
    x[i] = direction1[i];
    x[i+2] = direction2[i];
  }

  double temp = Genny.Rannyu(0.,4.);
  if(angular_distribution( x, parameters,"angles") < temp)
    return x;
  while(angular_distribution( x, parameters,"angles") > temp) {
    temp = Genny.Rannyu(0.,4.);
    direction1 = angle_unif(Genny,"angles");
    direction2 = angle_unif(Genny,"angles");
    for(int i=0; i<2;i++) {
      x[i] = direction1[i];
      x[i+2] = direction2[i];
    }
  }

  return x;

}


//this calculates the normalization factor in the lab FoR
double calc_norm_LAB(vector <double> dir) {
  return pow(pow(dir[0],2)+pow(dir[1],2)+pow(31.25*(dir[2]+1),2),0.5);
}

//this calculates the norm of a vector
double calc_norm(vector <double> dir) {
  return pow(pow(dir[0],2)+pow(dir[1],2)+pow(dir[2],2),0.5);
}

//this calculates the directions given the angles [theta,phi]
vector <double> calc_dir(vector <double> angles_LAB) {
  vector <double> temp {0,0,0};
  temp[0] = sin(angles_LAB[0])*cos(angles_LAB[1]);
  temp[1] = sin(angles_LAB[0])*sin(angles_LAB[1]);
  temp[2] = cos(angles_LAB[0]);
  return temp;
}

//this calculates the spin polarisation vector given the angles
vector <double> calc_spin(vector <double> angles, double alpha_1, double alpha_psi) {

  vector <double>  temp = {0,0,0};
  //directions of lambda_1
  double n1x = sin(angles[0])*cos(angles[1]);
  double n1y = sin(angles[0])*sin(angles[1]);
  double n1z = cos(angles[0]);

  //we can calculat s_2 with the directions of lambda_1
  temp[0]=alpha_1*(1-alpha_psi)/(3+alpha_psi)*n1x;
  temp[1]=-alpha_1*(1-alpha_psi)/(3+alpha_psi)*n1y;
  temp[2]=-alpha_1*(1+3*alpha_psi)/(3+alpha_psi)*n1z;

  return temp;
}

//this calculates the angular velocity of the precession given the direction where the particle is flying
vector <double> calc_precession(vector <double> direction,double g0) {
  vector <double>  temp = {0,0,0};
  double c =  0.99948783186249;
  long double muB_h=4.02782*pow(10,7);
  double g = g0;
  double Omega_x = muB_h*(g*(0-0.969*c*direction[1]*1*direction[0]));
  double Omega_y = muB_h*(g*(1-0.969*c*direction[1]*1*direction[1]));
  double Omega_z = muB_h*(g*(0-0.969*c*direction[1]*1*direction[2]));
  temp[0] = Omega_x;
  temp[1] = Omega_y;
  temp[2] = Omega_z;

  return temp;
}

//we need to rotate angular velocity from Lab frame to helicity frame
vector <double> rotate_omega_to_helicity(vector <double> omega, vector <double> angles_helicity_frame) {
  vector <double>  temp = {0,0,0};
  double beta = M_PI-angles_helicity_frame[0];
  double alfa = M_PI+angles_helicity_frame[1];
  vector <double>  xr = {cos(beta)*cos(alfa),sin(alfa)*cos(beta),-sin(beta)};
  vector <double>  yr = {-sin(alfa),cos(alfa),0};
  vector <double>  zr = {sin(beta)*cos(alfa),sin(alfa)*sin(beta),cos(beta)};
  temp[0] = dot_product(xr,omega);
  temp[1] = dot_product(yr,omega);
  temp[2] = dot_product(zr,omega);

  return temp;
}

//this calculates the spin after the precession in the magnetic field
vector <double> calc_spin_precession(vector <double> spin0, vector <double> omega, double Omega_norm, double t , ofstream& out_angles_precession) {
  vector <double>  temp = {0,0,0};
  double dot = dot_product(spin0, omega);
  vector <double>  cross = cross_product(spin0, omega);
  temp[0]=dot*omega[0]+(spin0[0]-dot*omega[0])*cos(Omega_norm*t)+cross[0]*sin(Omega_norm*t);
  temp[1]=dot*omega[1]+(spin0[1]-dot*omega[1])*cos(Omega_norm*t)+cross[1]*sin(Omega_norm*t);
  temp[2]=dot*omega[2]+(spin0[2]-dot*omega[2])*cos(Omega_norm*t)+cross[2]*sin(Omega_norm*t);
  out_angles_precession << calc_norm(temp) << endl;
  return temp;
}


//this is the angular distribution after the precession
double decay_spin_precess(vector <double> Omega,vector <double> Omega1,vector <double> Omega2,double g0, double alpha_psi, double alpha_1, double alpha_2, ofstream& out_angles_precession) {
  vector <double> angles  {Omega1[0],Omega1[1]}; //{theta_1, phi_1}
  vector <double> spin0 = calc_spin( angles, alpha_1, alpha_psi ); //initial spin value

  vector <double> direction = {0,0,0};
  direction = calc_dir(Omega);

  vector <double> direction_LAB  {0,0,0}; //calculate the direction of decaying lambda in the lab frame
  direction_LAB[0] = direction[0];
  direction_LAB[1] = direction[1];
  direction_LAB[2] = 31.25*(direction[2]+1);
  double norm = calc_norm(direction_LAB);
  direction_LAB[0] = direction_LAB[0]/norm;
  direction_LAB[1] = direction_LAB[1]/norm;
  direction_LAB[2] = direction_LAB[2]/norm;
  double t = 4/(299792458*0.9994878318624929*direction_LAB[2]);
  //cout << direction_LAB[0] << " " << direction_LAB[1] << " " << direction_LAB[2] << endl;
  //cout << calc_norm(direction_LAB) << endl;

  vector <double> Omega_velocity = {0,0,0};
  Omega_velocity=calc_precession(direction_LAB,g0);
  vector <double> omega = {0,0,0};
  omega[0]=Omega_velocity[0]/calc_norm(Omega_velocity);
  omega[1]=Omega_velocity[1]/calc_norm(Omega_velocity);
  omega[2]=Omega_velocity[2]/calc_norm(Omega_velocity);
  omega = rotate_omega_to_helicity(omega, Omega);
  //cout << calc_norm(omega) << endl;

  vector <double> spin_precess = calc_spin_precession(spin0,omega,calc_norm(Omega_velocity),t, out_angles_precession);
  vector <double> n2 = {sin(Omega2[0])*cos(Omega2[1]),sin(Omega2[0])*sin(Omega2[1]),cos(Omega2[0]) };
  return 1-alpha_2*dot_product(spin_precess,n2);
}

vector <double> accept_reject_direct(Random& Genny, vector <double> par, ofstream& out_angles_precession) {

  double alpha_psi = par[0];
  double alpha_1 = par[1];
  double alpha_2 = par[2];
  double g = par[3];

  vector <double> angles = {0,0,0,0,0,0,0};
  vector <double>  Omega = angle_unif(Genny,"angles");
  angles[0] = Omega[0];
  angles[1] = Omega[1];
  vector <double>  Omega1 = angle_unif(Genny,"angles");
  angles[2] = Omega1[0];
  angles[3] = Omega1[1];
  vector <double>  Omega2 = angle_unif(Genny,"angles");
  angles[4] = Omega2[0];
  angles[5] = Omega2[1];
  double maxi = Genny.Rannyu()*1.5;
  double y_val = decay_spin_precess(Omega, Omega1, Omega2,g,alpha_psi,alpha_1,alpha_2,out_angles_precession);
  if (maxi<=y_val) {
      angles[6] = y_val;
      return angles;
    }
  else {
    while(maxi>y_val) {
      Omega = angle_unif(Genny,"angles");
      angles[0] = Omega[0];
      angles[1] = Omega[1];
      Omega1 = angle_unif(Genny,"angles");
      angles[2] = Omega1[0];
      angles[3] = Omega1[1];
      Omega2 = angle_unif(Genny,"angles");
      angles[4] = Omega2[0];
      angles[5] = Omega2[1];
      maxi = Genny.Rannyu()*1.5;
      y_val = decay_spin_precess(Omega, Omega1, Omega2,g,alpha_psi,alpha_1,alpha_2,out_angles_precession);
    }
    angles[6] = y_val;
    return angles;
  }
}


#endif
