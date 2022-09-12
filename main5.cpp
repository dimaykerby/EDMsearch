#include"random.h"
#include"lib.h"

using namespace std;

int main() {

//random number generator
Random ran_gen ;
set_random_gen(ran_gen);

//define the angular distribution
//TF1 *ang_dist = new TF1("ang_dist", angular_distribution, xmin, xmax, npar );

string method;
cout << "direct or indirect method?" << endl;
cin >> method;

if( method == "indirect") {

  ofstream out_directions_y;
  //parameters
  vector <double> par = {0,0,0};
  int n_events = 45*pow(10,5);
  int n_pseudoex = 1;
  vector <double> extracted;

  double y_value = 0.;

  for( int i = 0; i < n_pseudoex ; i++) {
    cout << endl;
    cout << i+1 << "/" << n_pseudoex << endl;
    out_directions_y.open("angles/"+to_string(i+1+100)+".csv");
    par = {ran_gen.Gauss(0.461,0.013), ran_gen.Gauss(0.750,0.013),-ran_gen.Gauss(0.750,0.013)}; //{alpha_psi,alpha_1,alpha_2}
    for( int j = 0; j < n_events; j++) {

      if(j%(int(45*pow(10,4)))==0)
        cout << j/(45*pow(10,3)) << "%" << endl;

      extracted = accept_reject_indirect(ran_gen,par);
      y_value = angular_distribution(extracted,par,"angles");

      for( int k=0; k < 4; k++)
        out_directions_y << extracted[k] << ",";
      out_directions_y << y_value << endl;

    }
    out_directions_y.close();
  }

}

else {
  //ofstream out_directions_y;
  ofstream out;
  out.open("spin_precessed.csv");
  //parameters
  vector <double> par = {0,0,0,0}; //{alpha_psi,alpha_1,alpha_2,g}
  int n_events = 292225;
  int n_pseudoex = 1;
  vector <double> extracted;

  for( int i = 0; i < n_pseudoex ; i++) {
    cout << endl;
    cout << i+1 << "/" << n_pseudoex << endl;
    //out_directions_y.open("angles_direct/"+to_string(i+1/*+100*/)+".csv");
      //out_directions_y.open("s_0&s_precess.csv");
    par = {ran_gen.Gauss(0.461,0.013), ran_gen.Gauss(-0.750,0.013),ran_gen.Gauss(0.750,0.013),-1.458};
    //par = {0.461, -0.750, 0.750,-1.458};
    for( int j = 0; j < n_events; j++) {

      if(j%(int(n_events/10))==0)
        cout << j/(n_events/100) << "%" << endl;

      extracted = accept_reject_direct(ran_gen,par,out);

      /*for( int k=0; k < 6; k++)
        out_directions_y << extracted[k] << ",";
      out_directions_y << extracted[6] << endl;*/

    }
    //out_directions_y.close();
  }
  }
return 0 ;
}
