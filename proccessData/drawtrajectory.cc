/*
drawtrajectory.cc
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

  Boundary* bd = new Boundary(3.0);

  const int tend = 4399;
  const int tspan = 50;


  for (int vl = 2; vl <=8; vl+=2) {
    for (int ar = 510; ar <= 560; ar+=10) {
      vector< vector<double> > pos(particles, vector<double>(3*tspan));//pos[id][3*t+j]j=0->x,j=1->y,j=2->z

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
      }
      //cout<<"imput end"<<ar<<" "<<vl<<endl;

      const string output1 = ".png";
      const string output = imput1+"tra"+imput2+to_string(ar)+imput3+to_string(vl)+imput4+output1;
      FILE *gp;

      gp = popen("gnuplot -persist", "w");
      fprintf(gp, "load \"~/Dropbox/GnuPlotPalette/gnuplot-palettes/rdpu.pal\" \n");
      fprintf(gp, "unset colorbox\n");
      fprintf(gp, "set term png size 1200, 1200\n");
      fprintf(gp, "set output \"%s\" \n", output.c_str());

      fprintf(gp, "sp[%f:%f][%f:%f][%f:%f]\'-\' u 1:2:3:4 w p pt 7 ps 1 lc pal notitle \n", bd->xmin_, bd->xmax_, bd->ymin_, bd->ymax_, bd->zmin_, bd->zmax_);
      //fprintf(gp, "sp[%f:%f][%f:%f][%f:%f]\'-\' u 1:2:3:4 w l lw 2 lc pal notitle \n", bd->xmin_, bd->xmax_, bd->ymin_, bd->ymax_, bd->zmin_, bd->zmax_);
      
      for ( int id = 0; id != particles; ++id) {
        for (int ti = 0; ti < tspan-1; ++ti) {
          double color = (double)id/(double)particles;

          double r = ( (pos[id][3*ti]-pos[id][3*(ti+1)])*(pos[id][3*ti]-pos[id][3*(ti+1)]) + (pos[id][3*ti+1]-pos[id][3*(ti+1)+1])*(pos[id][3*ti+1]-pos[id][3*(ti+1)+1]) + (pos[id][3*ti+2]-pos[id][3*(ti+1)+2])*(pos[id][3*ti+2]-pos[id][3*(ti+1)+2]) );
          if (r < 1.0){ 
            fprintf(gp,"%f %f %f %f\n", pos[id][3*ti], pos[id][3*ti+1], pos[id][3*ti+2], color);
            fprintf(gp,"%f %f %f %f\n", pos[id][3*(ti+1)], pos[id][3*(ti+1)+1], pos[id][3*(ti+1)+2], color);
          } else{
            fprintf(gp, "\n");
          }
          fprintf(gp, "\n");
        }
      }
      
      fprintf(gp,"e\n");
      fprintf(gp,"set output\n");
      fclose(gp);
      cout<<"drawing end "<<ar<<" "<<vl<<endl;
    }
  }

  return 0;
}

