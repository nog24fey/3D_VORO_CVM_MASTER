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
  
// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;
// Set the number of particles that are going to be randomly introduced
const int particles=6*6*6;
 

int main(int argc, char **argv) {
  //system parameters
  const int i_tarea = atoi(argv[3]);
  const int i_smpl = atoi(argv[4]);
  const double tarea = 0.01*(double)i_tarea;

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

  
  container con(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,n_x,n_y,n_z,
		true,true,true,8);
    
  // Randomly add particles into the container
  setInitialConfiguration(con, vp, bd, mt);

  con.print_custom("%i %v","packing.custom3");

  int endtime = 20;
  int msdstarttime = endtime/20;
  double hsttics = (double)(particles*(endtime-msdstarttime));
  hsttics = 1.0/hsttics;
  for (int time = 0; time != endtime; ++time) {

    container eon(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,n_x,n_y,n_z,
		  true,true,true,8);

    //(random+selfpropel)motion
    for (int i=0; i<particles; ++i) {
      //determine direction for selfpropel
      double dh = Pi*rnd(mt);
      vp[i].theta_ = dh;
      //while (theta[i]<0.0) theta[i] += Pi;
      //while (theta[i]>=Pi) theta[i] -= Pi;
      dh = 2.0*Pi*rnd(mt);
      vp[i].phi_ = dh;
      //if (phi[i]<0.0) phi[i] += 2.0*Pi;
      //else if (phi[i]>=2.0*Pi) phi[i] -= 2.0*Pi;
      //motion
      vp[i].x_ += vp[i].r_*sin(vp[i].theta_)*cos(vp[i].phi_);
      vp[i].y_ += vp[i].r_*sin(vp[i].theta_)*sin(vp[i].phi_);
      vp[i].z_ += vp[i].r_*cos(vp[i].theta_);
      
      eon.put(i,vp[i].x_,vp[i].y_,vp[i].z_);
    }
    
    double tot_eng = retTotalEnergy(eon, tarea);

    double ox, oy, oz, vx, vy, vz;
    for (int i = 0;i<particles;i++) {

      ox = vp[i].x_;
      oy = vp[i].y_;
      oz = vp[i].z_;
      vx = vp[i].x_+dpos;
      vy = vp[i].y_+dpos;
      vz = vp[i].z_+dpos;

      container xon(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,n_x,n_y,n_z,
		    true,true,true,8);
      for (int j = 0;j<particles;j++) {
	if (i!=j) {
	  xon.put(j,vp[j].x_,vp[j].y_,vp[j].z_);
	}
      }
      xon.put(i,vx,oy,oz);
      double tot_eng_x = retTotalEnergy(xon, tarea);
      vp[i].xn_ = -0.01*(tot_eng_x-tot_eng)/dpos;

      
      container yon(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,n_x,n_y,n_z,
		    true,true,true,8);
      for (int j = 0;j<particles;j++) {
	if (i!=j) {
	  yon.put(j,vp[j].x_,vp[j].y_,vp[j].z_);
	}
      }      
      yon.put(i,ox,vy,oz);
      double tot_eng_y = retTotalEnergy(yon, tarea);
      vp[i].yn_ = -0.01*(tot_eng_y-tot_eng)/dpos;

      
      container zon(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,n_x,n_y,n_z,
		    true,true,true,8);
      for (int j = 0;j<particles;j++) {
	if (i!=j) {
	  zon.put(j,vp[j].x_,vp[j].y_,vp[j].z_);
	}
      }
      zon.put(i,ox,oy,vz);
      double tot_eng_z = retTotalEnergy(zon, tarea);
      vp[i].zn_ = -0.01*(tot_eng_z-tot_eng)/dpos;

    }

    //renew position
    container don(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,n_x,n_y,n_z,
		  true,true,true,8);
    renewPosition(don, vp, bd);


    if (time > msdstarttime) {
      msd<<time-msdstarttime<<" ";
      double meandiff = 1.0;
      for ( int i=0; i<particles; ++i) {
	meandiff *= pow(squaredDistanceForMSD(vp[i].x_,vp[i].y_,vp[i].z_,vp[i].x0_,vp[i].y0_,vp[i].z0_,bd->x_axe_leng_,bd->y_axe_leng_,bd->z_axe_leng_),1.0/(double)particles);
        //meandiff = min(meandiff, squaredDistanceForMSD(vp[i].x_,vp[i].y_,vp[i].z_,x0[i],y0[i],z0[i],bd->x_axe_leng_,bd->y_axe_leng_,bd->z_axe_leng_));
      }
      //meandiff /= (double)particles;
      msd<<meandiff<<endl;

      c_loop_all cmh(don);
      voronoicell_neighbor ch;
      if(cmh.start()) do if(don.compute_cell(ch,cmh)) {
	    hstdata[ch.number_of_faces()] += hsttics;
	  }while (cmh.inc());
    }
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
      fprintf(gp, "sp [%f:%f][%f:%f][%f:%f]\"%s\" u 2:3:4 w p ps 2 pt 7 lc rgbcolor \"dark-green\" ti \"vpos\", \"%s\" u 1:2:3 w l lw 0.2 lc rgbcolor \"gray50\" ti \"edge\" \n", bd->xmin_, bd->xmax_, bd->ymin_, bd->ymax_, bd->zmin_, bd->zmax_, str_p.c_str(), str_v.c_str());
      fprintf(gp, "set output\n");
      pclose(gp);
      
    }
    //for (auto & h : hstdata) cout<<h<<" ";
    //cout<<endl;
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
