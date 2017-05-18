/*
intermediatescattering.cc
*/

#include "../src/spaceutility.h"
#include "../src/voro_func.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cassert>

using namespace std;

const int particles=6*6*6;

int main(int argc, char **argv){

  const string imput1 = "/home/dein/Desktop/TD0419/res/";
  const string imput2 = "ar";
  const string imput3 = "sm1df1vl";
  const string imput4 = "et22000ts10point";
  const string imput5 = ".dat";

  Boundary* bd = new Boundary(3.0);

  const int tend = 4399;
  const int tspan = 50;


  for (int vl = 2; vl <=8; vl+=2) {
    for (int ar = 510; ar <= 560; ar+=10) {
      vector< vector<double> > pos(particles, vector<double>(3*tspan));//pos[id][3*t+j]j=0->x,j=1->y,j=2->z

      int bin = 200;
      vector<double> gr(bin,0.0); //spacetic=0.01
      
      for (int ti = tend-tspan+1; ti <= tend; ++ti) {

        const string imput = imput1+to_string(ti)+imput2+to_string(ar)+imput3+to_string(vl)+imput4+imput5;
        //cout<<imput<<endl;
        ifstream ifs;
        ifs.open(imput.c_str());
        double x,y,z,g;
        int tid = 3*(ti-(tend-tspan+1));
        
        for( int id = 0; id !=particles; ++id) {

          ifs>>x;
          ifs>>y;
          ifs>>z;

          //cout<<id<<" "<<tid<<endl;
          pos[id][tid] = x;
          pos[id][tid+1] = y;
          pos[id][tid+2] = z;
          
          ifs>>g; ifs>>g; ifs>>g;
        }
        ifs.close();

        //cout<<"imput end"<<endl;

        for (int id = 0; id != particles-1; ++id) {
          for (int jd = id+1; jd != particles; ++jd) {
            double r = sqrt(squaredDistanceForMSD(pos[id][tid], pos[id][tid+1], pos[id][tid+2],pos[jd][tid], pos[jd][tid+1], pos[jd][tid+2],bd->x_axe_leng_,bd->y_axe_leng_,bd->z_axe_leng_));//suppose [0-10]
            //cout<<id<<" "<<jd<<" "<<r<<endl;
            assert(r<10.0);
            gr[(int)( floor( (r*(double)bin/10.0) ) )] += 1.0;
          }
        }
        
      }

      /*
      auto iter = std::max_element(gr.begin(), gr.end());
      auto index = std::distance(gr.begin(), iter);
      double r0 = (double)( (int)index)*10.0/(double)bin;
      */
      
      //const string output1 = ".dat";
      //const string output = imput1+"tra"+imput2+to_string(ar)+imput3+to_string(vl)+imput4+output1;

      int ii = 0;
      for( auto &g : gr) {
        double rr = (double)ii*10.0/(double)bin;
        g = 4.0*3.141592*rr*g;
        ++ii;
        //cout<<rr<<" "<<g<<endl;
      }
      //cout<<endl;

      double r0 = 0.0;
      //for ( int bi = (int)((double)bin*sqrt(3.0)/2.0); bi != bin; ++bi) gr[bi] = 0.0;
      for ( int i = 0; i != bin-1; ++i) {
        if (gr[i] > gr[i+1]) {
          r0 = (double)i*10.0/(double)bin;
          break;
        }
      }
      
      cout<<"result "<<ar<<" "<<vl<<" "<<r0<<endl;
    }
  } 

  return 0;
}

