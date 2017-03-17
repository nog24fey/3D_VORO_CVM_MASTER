/*
argv[1] directory name packing data files !!must be ended with "/" !!
argv[2] directory name packing image files !!must be ended with "/" !!
argv[3] target value of area; timed 0.01 
argv[4] sample id; also used for random seed
argv[5] (int)time step
*/

// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011
#include "./src/directorymake.h"
#include "./src/spaceutility.h"
#include "./src/voro_func.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

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
  

// Set the number of particles that are going to be randomly introduced
const int particles=6*6*6;
 

int main(int argc, char **argv) {
  //system parameters
  const int i_tarea = atoi(argv[3]);
  const int i_smpl = atoi(argv[4]);
  const double tarea = 0.01*(double)i_tarea;
  const int timestep = atoi(argv[5]);
  const double timetics = 1.0/(double)timestep;

  //data packing directories
  directoryMake(argv[1]);
  directoryMake(argv[2]);
  string datdirectoryname = argv[1];
  string imagedirectoryname = argv[2];

  //filename settings
  const string str_tarea = "ar"+to_string(i_tarea);
  const string str_smpl = "sm"+to_string(i_smpl);
  const string str_param = str_tarea+str_smpl;

  //open msd dat file
  string str_msd = datdirectoryname+"MSD"+str_param+".dat";
  ofstream msd;
  msd.open(str_msd);

  //open histgram dat file
  string str_hst = datdirectoryname+"HST"+str_param+".dat";
  ofstream hst;
  hst.open(str_hst);
  vector<double> hstdata(40,0.0);
  
  //initialize random function
  mt.seed(i_smpl);
  uniform_real_distribution<double> rnd(0.0,1.0);

  //spacetics
  double dpos = 0.01;

  //boundary
  Boundary* bd = new Boundary(3.0);
  //internal variables
  vector<VoronoiPoint> vp;
  for( int i = 0; i != particles; ++i) vp.push_back(VoronoiPoint());

  
  container con(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
		true,true,true,8);
    
  // Randomly add particles into the container
  setInitialConfiguration(con, vp, bd, mt);

  con.print_custom("%i %v","packing.custom3");

  int endtime = 2000;
  int msdstarttime = endtime/20;
  double hsttics = (double)(particles*(endtime-msdstarttime));
  hsttics = 1.0/hsttics;
  for (int time = 0; time != endtime; ++time) {

     for( int ts = 0; ts != timestep; ++ ts) {
        container eon(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
		  true,true,true,8);
        execLangevinStep(eon, vp, bd, mt, tarea, dpos, timetics);
     }

     container don(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
                   true,true,true,8);
     setContainer(don, vp);

     if (time > msdstarttime) {
        msd<<time-msdstarttime<<" ";
        double meandiff = 1.0;
        for ( int i=0; i<particles; ++i) {
           meandiff *= pow(squaredDistanceForMSD(vp[i].x_,vp[i].y_,vp[i].z_,vp[i].x0_,vp[i].y0_,vp[i].z0_,bd->x_axe_leng_,bd->y_axe_leng_,bd->z_axe_leng_),1.0/(double)particles);
           //meandiff += squaredDistanceForMSD(vp[i].x_,vp[i].y_,vp[i].z_,vp[i].x0_,vp[i].y0_,vp[i].z0_,bd->x_axe_leng_,bd->y_axe_leng_,bd->z_axe_leng_);
           //meandiff = min(meandiff, squaredDistanceForMSD(vp[i].x_,vp[i].y_,vp[i].z_,x0[i],y0[i],z0[i],bd->x_axe_leng_,bd->y_axe_leng_,bd->z_axe_leng_));
        }
           //meandiff /= (double)particles;
        msd<<meandiff<<endl;
     }
     
     c_loop_all cmh(don);
     voronoicell_neighbor ch;
     if(cmh.start()) do if(don.compute_cell(ch,cmh)) {
              hstdata[ch.number_of_faces()] += hsttics;
           }while (cmh.inc());
 
     //start calculation msd;
     if (time == msdstarttime) {
        for (int i = 0; i < particles; ++i) {
           vp[i].x0_ = vp[i].x_;
           vp[i].y0_ = vp[i].y_;
           vp[i].z0_ = vp[i].z_;
        }
     }
    
     if (time%5 == 0) {
        string stime(to_string(time/5));
        string str_p = datdirectoryname+stime+str_param+"point.dat";
        string str_v = datdirectoryname+stime+str_param+"edges.dat";
        don.draw_particles(str_p.c_str());
        don.draw_cells_gnuplot(str_v.c_str());

        string str_pvj = imagedirectoryname+stime+str_param+"pointedges.png";
      
        FILE* gp;
        gp = popen("gnuplot -persist","w");
        fprintf(gp, "set term png size 1200, 1200\n");

        fprintf(gp, "set output \"%s\" \n", str_pvj.c_str());
        fprintf(gp, "sp Type& operator[](int index);%f:%f][%f:%f][%f:%f]\"%s\" u 2:3:4 w p ps 2 pt 7 lc rgbcolor \"dark-green\" ti \"vpos\", \"%s\" u 1:2:3 w l lw 0.2 lc rgbcolor \"gray50\" ti \"edge\" \n", bd->xmin_, bd->xmax_, bd->ymin_, bd->ymax_, bd->zmin_, bd->zmax_, str_p.c_str(), str_v.c_str());
        fprintf(gp, "set output\n");
        pclose(gp);
      
     }
  }
  msd.close();

  int hi = 0;
  for (auto & h : hstdata) {
    hst<<hi<<" "<<h<<endl;
    ++hi;
  }
  hst<<endl;
  hst.close();

  delete bd;
  return 1;
 }
