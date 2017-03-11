/*
argv[1] directory name packing data files !!must be ended with "/" !!
argv[2] directory name packing image files !!must be ended with "/" !!
argv[3] target value of area; timed 0.01 
argv[4] sample id; also used for random seed
argv[5] 
*/

// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011
  
#include "./voro++/voro++.hh"
#include "./src/directorymake.h"
#include "./src/spaceutility.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#define Pi 3.14159265358979323846

using std::vector;
using std::string;
using std::cout;
using std::ofstream;
using std::endl;
using std::to_string;
using std::sin;
using std::cos;
using std::sqrt;

using std::mt19937;
using std::uniform_real_distribution;
mt19937 mt;

using namespace voro;
  
// Set up constants for the container geometry
const double x_min=-2.0,x_max=2.0;
const double y_min=-2.0,y_max=2.0;
const double z_min=-2.0,z_max=2.0;

const double x_axe_leng = x_max-x_min;
const double y_axe_leng = y_max-y_min;
const double z_axe_leng = z_max-z_min;

const double cvol=x_axe_leng*y_axe_leng*z_axe_leng;
 
// Set up the number of blocks that the container is divided into
const int n_x=4,n_y=4,n_z=4;
 
// Set the number of particles that are going to be randomly introduced
const int particles=4*4*4;
 

int main(int argc, char **argv) {
  //system parameters


  //data packing directories
  directoryMake(argv[1]);
  directoryMake(argv[2]);
  string datdirectoryname = argv[1];
  string imagedirectoryname = argv[2];

  //filename settings
  
  //initialize random function
  mt.seed(0);
  uniform_real_distribution<double> rnd(0.0,1.0);

  //spacetics
  //internal variables
 

  for ( int si = 0; si != 50; ++si) {
    
    container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                  true,true,true,8);
    
    for (int i = 0; i != particles; ++i) {
      double x = x_min+rnd(mt)*(x_max-x_min);
      double y = y_min+rnd(mt)*(y_max-y_min);
      double z = z_min+rnd(mt)*(z_max-z_min);

      con.put(i,x,y,z);
    }
    string pack = datdirectoryname+"packing"+to_string(si);
    con.print_custom("%i %v",pack.c_str());

    c_loop_all cl(con);
    voronoicell_neighbor c;
    if(cl.start()) do if(con.compute_cell(c,cl)) {
          int numfaces = c.number_of_faces();
          double vl = c.volume();
          double ar = c.surface_area();
          cout<<si<<" "<<numfaces<<" "<<vl<<" "<<ar<<endl;
        }while (cl.inc());


    string str_p = datdirectoryname+to_string(si)+"point.dat";
    string str_v = datdirectoryname+to_string(si)+"edges.dat";
    con.draw_particles(str_p.c_str());
    con.draw_cells_gnuplot(str_v.c_str());
    string str_pvj = imagedirectoryname+to_string(si)+"pointedges.png";
    FILE* gp;
    gp = popen("gnuplot -persist","w");
    fprintf(gp, "set term png size 1200, 1200\n");

    fprintf(gp, "set output \"%s\" \n", str_pvj.c_str());
    fprintf(gp, "sp [%f:%f][%f:%f][%f:%f]\"%s\" u 2:3:4 w p ps 2 pt 7 lc rgbcolor \"dark-green\" ti \"vpos\", \"%s\" u 1:2:3 w l lw 0.2 lc rgbcolor \"gray50\" ti \"edge\" \n", x_min, x_max, y_min, y_max, z_min, z_max, str_p.c_str(), str_v.c_str());
    fprintf(gp, "set output\n");
    pclose(gp);
  }
    


  return 1;
 }
