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
const double x_min=-3.0,x_max=3.0;
const double y_min=-3.0,y_max=3.0;
const double z_min=-3.0,z_max=3.0;

const double x_axe_leng = x_max-x_min;
const double y_axe_leng = y_max-y_min;
const double z_axe_leng = z_max-z_min;

const double cvol=x_axe_leng*y_axe_leng*z_axe_leng;
 
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
  //internal variables
  vector<double> x(particles),y(particles),z(particles),theta(particles),phi(particles),r(particles,0.01);
  vector<double> x0(particles),y0(particles),z0(particles);
  vector<double> xold(particles),yold(particles),zold(particles);
  vector<double> xnew(particles),ynew(particles),znew(particles);
 
  // Create a container with the geometry given above, and make it
  // non-periodic in each of the three coordinates. Allocate space for
  // eight particles within each computational block
  container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
		true,true,true,8);

    
  // Randomly add particles into the container
  for(int i = 0;i<particles;i++) {
    theta[i] = rnd(mt)*Pi;
    phi[i] = rnd(mt)*2.0*Pi;
    
    x[i]=x_min+rnd(mt)*(x_max-x_min);
    y[i]=y_min+rnd(mt)*(y_max-y_min);
    z[i]=z_min+rnd(mt)*(z_max-z_min);
    con.put(i,x[i],y[i],z[i]);
  }

  con.print_custom("%i %v","packing.custom3");

  int endtime = 20000;
  int msdstarttime = endtime/20;
  double hsttics = (double)(particles*(endtime-msdstarttime));
  hsttics = 1.0/hsttics;
  for (int time = 0; time != endtime; ++time) {

    container eon(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
		  true,true,true,8);

    //(random+selfpropel)motion
    for (int i=0; i<particles; ++i) {
      //determine direction for selfpropel
      double dh = Pi*rnd(mt);
      theta[i] = dh;
      //while (theta[i]<0.0) theta[i] += Pi;
      //while (theta[i]>=Pi) theta[i] -= Pi;
      dh = 2.0*Pi*rnd(mt);
      phi[i] = dh;
      //if (phi[i]<0.0) phi[i] += 2.0*Pi;
      //else if (phi[i]>=2.0*Pi) phi[i] -= 2.0*Pi;
      //motion
      x[i] += r[i]*sin(theta[i])*cos(phi[i]);
      y[i] += r[i]*sin(theta[i])*sin(phi[i]);
      z[i] += r[i]*cos(theta[i]);
      
      eon.put(i,x[i],y[i],z[i]);
    }
    
    double tot_eng = 0.0;
    c_loop_all cl(eon);
    voronoicell_neighbor c;
    if(cl.start()) do if(eon.compute_cell(c,cl)) {
	  double ve = (c.volume()-1.0);
	  double ae = (c.surface_area()-tarea);
	  //cout<<"ae"<<ae<<endl;
	  tot_eng += ve*ve+ae*ae;
	}while (cl.inc());
    //cout<<tot_eng<<endl;

    double ox, oy, oz, vx, vy, vz;
    for (int i = 0;i<particles;i++) {

      ox = x[i];
      oy = y[i];
      oz = z[i];
      vx = x[i]+dpos;
      vy = y[i]+dpos;
      vz = z[i]+dpos;

      container xon(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
		    true,true,true,8);
      for (int j = 0;j<particles;j++) {
	if (i!=j) {
	  xon.put(j,x[j],y[j],z[j]);
	}
      }

      xon.put(i,vx,oy,oz);
      double tot_eng_x = 0.0;
      c_loop_all cmx(xon);
      voronoicell_neighbor cx;
      if(cmx.start()) do if(xon.compute_cell(cx,cmx)) {
	    double ve = (cx.volume()-1.0);
	    double ae = (cx.surface_area()-tarea);
	    tot_eng_x += ve*ve+ae*ae;
	  }while (cmx.inc());
      xnew[i] = -0.01*(tot_eng_x-tot_eng)/dpos;

      container yon(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
		    true,true,true,8);
      for (int j = 0;j<particles;j++) {
	if (i!=j) {
	  yon.put(j,x[j],y[j],z[j]);
	}
      }
      
      yon.put(i,ox,vy,oz);
      double tot_eng_y = 0.0;
      c_loop_all cmy(yon);
      voronoicell_neighbor cy;
      if(cmy.start()) do if(yon.compute_cell(cy,cmy)) {
	    double ve = (cy.volume()-1.0);
	    double ae = (cy.surface_area()-tarea);
	    tot_eng_y += ve*ve+ae*ae;
	  }while (cmy.inc());
      ynew[i] = -0.01*(tot_eng_y-tot_eng)/dpos;

      container zon(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
		    true,true,true,8);
      for (int j = 0;j<particles;j++) {
	if (i!=j) {
	  zon.put(j,x[j],y[j],z[j]);
	}
      }

      zon.put(i,ox,oy,vz);
      double tot_eng_z = 0.0;
      c_loop_all cmz(zon);
      voronoicell_neighbor cz;
      if(cmz.start()) do if(zon.compute_cell(cz,cmz)) {
	    double ve = (cz.volume()-1.0);
	    double ae = (cz.surface_area()-tarea);
	    tot_eng_z += ve*ve+ae*ae;
	  }while (cmz.inc());
      znew[i] = -0.01*(tot_eng_z-tot_eng)/dpos;
    }

    //renew position
    container don(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
		  true,true,true,8);
    double cx, cy, cz;
    cx = 0.0; cy = 0.0; cz = 0.0;
    for (int i = 0; i<particles; i++) {
      xold[i] = x[i]; yold[i] = y[i]; zold[i] = z[i];
      x[i] += xnew[i]; y[i] += ynew[i]; z[i] += znew[i];
      cx += x[i]-xold[i]; cy += y[i]-yold[i]; cz += z[i]-zold[i];
    }
    cx /= (double)particles; cy /= (double)particles; cz /= (double)particles;
    for (int i = 0; i<particles; i++) {
      x[i] -= cx; y[i] -= cy; z[i] -= cz;

      while (x[i] < x_min) x[i] += x_axe_leng;
      while (x[i] >= x_max) x[i] -= x_axe_leng;
      while (y[i] < y_min) y[i] += y_axe_leng;
      while (y[i] >= y_max) y[i] -= y_axe_leng;
      while (z[i] < z_min) z[i] += z_axe_leng;
      while (z[i] >= z_max) z[i] -= z_axe_leng;

      don.put(i,x[i],y[i],z[i]);
    }

    if (time > msdstarttime) {
      msd<<time-msdstarttime<<" ";
      double meandiff = 1.0;
      for ( int i=0; i<particles; ++i) {
	meandiff *= pow(squaredDistanceForMSD(x[i],y[i],z[i],x0[i],y0[i],z0[i],x_axe_leng,y_axe_leng,z_axe_leng),1.0/(double)particles);
        //meandiff = min(meandiff, squaredDistanceForMSD(x[i],y[i],z[i],x0[i],y0[i],z0[i],x_axe_leng,y_axe_leng,z_axe_leng));
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
	x0[i] = x[i];
	y0[i] = y[i];
	z0[i] = z[i];
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
      fprintf(gp, "sp [%f:%f][%f:%f][%f:%f]\"%s\" u 2:3:4 w p ps 2 pt 7 lc rgbcolor \"dark-green\" ti \"vpos\", \"%s\" u 1:2:3 w l lw 0.2 lc rgbcolor \"gray50\" ti \"edge\" \n", x_min, x_max, y_min, y_max, z_min, z_max, str_p.c_str(), str_v.c_str());
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
  return 1;
 }
