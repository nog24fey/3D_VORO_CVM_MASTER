/*
hstfaceperimeter.cc
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

using namespace std;

const int particles=6*6*6;

int main(int argc, char **argv){

  const string imput1 = "/home/dein/Desktop/TD0419/res/";
  const string imput2 = "ar";
  const string imput3 = "sm1df1vl";
  const string imput4 = "et22000ts10point";
  const string imput5 = ".dat";
  const string output1 = "/home/dein/Desktop/TD0419/histfaceperimeter";

  Boundary* bd = new Boundary(3.0);

  const int tend = 4399;
  const int tspan = 100;


  const int numFaces = 40;
  const int numBinPeri = 500;
  const double maxPeri = 50.0;

  const double tic = 1.0/( (double)(tspan*particles) );
  for (int vl = 2; vl <=8; vl+=2) {
    for (int ar = 510; ar <= 560; ar+=10) {
      vector< vector<double> > density(numFaces, vector<double>(numBinPeri, 0.0));
      
      for (int ti = tend-tspan+1; ti <= tend; ++ti) {

        const string imput = imput1+to_string(ti)+imput2+to_string(ar)+imput3+to_string(vl)+imput4+imput5;
        //cout<<imput<<endl;
        ifstream ifs;
        ifs.open(imput.c_str());
        double x,y,z,g;

        container con(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
                      true,true,true,8);
        
        for( int id = 0; id !=particles; ++id) {

          ifs>>x;
          ifs>>y;
          ifs>>z;

          //cout<<id<<" "<<tid<<endl;
          con.put(id,x,y,z);
          
          ifs>>g; ifs>>g; ifs>>g;
        }
        ifs.close();

        c_loop_all cm(con);
        voronoicell_neighbor cl;
        if(cm.start()) do if(con.compute_cell(cl,cm)) {
              int faces = cl.number_of_faces();
              double peri = cl.total_edge_distance();
              density[faces][(int)( floor(peri*(double)numBinPeri/maxPeri) )] += tic;
            } while (cm.inc());
        
      }
      const string output = output1+imput2+to_string(ar)+imput3+to_string(vl)+imput4+imput5;
      ofstream ofs;
      ofs.open(output.c_str());
      
      for ( int i = 0; i != numFaces; ++i) {
        for ( int j = 0; j != numBinPeri; ++j) { 
          ofs<< i <<" "<< (double)j*maxPeri/(double)numBinPeri <<" "<<density[i][j]<<endl;
        }
        ofs<<endl;
      }
      ofs<<endl;
      cout<<"stats end "<< ar<<" "<<vl<<endl;
    }
  }

  return 0;
}

