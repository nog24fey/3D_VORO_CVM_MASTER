/*
msdcoeffscaled.cc  
*/
#include "./regression1st.h"
#include "./regression2nd.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

int main(int argc, char **argv){

  const string imput1 = "/home/dein/Desktop/TD0419/MSDar";
  const string imput2 = "sm2df1vl";
  const string imput3 = "et22000ts10.dat";

  for (int vl = 2; vl <= 8; vl+=2) {
    for (int ar = 510; ar <= 560; ar+=10) {      
      const string imput = imput1+to_string(ar)+imput2+to_string(vl)+imput3;
      ifstream ifs;
      ifs.open(imput.c_str());

      const double vlreal = 0.001*(double)vl;
      const int D = 1;
      const double Dreal = 0.1*(double)D;
      double scalefactor = vlreal*vlreal/(3*Dreal);
      
      Regression1st *z1 = new Regression1st();
      //Regression2nd *z2 = new Regression2nd();

      //double x0;
      double x = 0.0; double y = 0.0;
      const double min_y = 0.01;
      const double max_y = 1.0;

      //bool flagx0 = false;

      while ( !ifs.eof() ) {
        ifs>>x; ifs>>y;
        //cout<<x<<" "<<y<<endl;
        if ( (y > min_y) && (y < max_y) ) {
          //cout<<x<<" "<<y<<endl;
          z1->pushData(x, y);
          //z2->pushData(x, 1000*y);
        }
      }
      /*
        while ( !ifs.eof() ) {
        ifs>>x; ifs>>y;
        if ( flagx0 == false && y > min_y) {
        x0 = x;
        flagx0 = true;
        }
        if ( (y > min_y) && (y < max_y) ) {
        z1->pushData(log(x-x0+1), log(y));
        z2->pushData(log(x-x0+1), log(y));
        }
        }
      */
      double a, b;
      //double a, b, c;
      z1->getResult(a, b);
      double Deff = a/(6.0*scalefactor);
      cout<<vl<<" "<<ar<<" "<<Deff<<" "<<b<<" ";
      if ( Deff > 0.025 ) {
        cout<<1;
      }else {
        cout<<0;
      }
      cout<<endl;
      //z2->getResult(a, b, c);
      //cout<<"2nd order model... "<<a<<" "<<b<<" "<<c<<endl;
      delete z1;
      //delete z2;
      ifs.close();
    }
    cout<<endl;
  }
  return EXIT_SUCCESS;	
}
