/*
argv[1] directory name packing data files !!must be ended with "/" !!
argv[2] directory name packing image files !!must be ended with "/" !!
argv[3] target value of area; timed 0.01           eg. 300-1000
argv[4] sample id; also used for random seed       eg. 0-10
argv[5] (int)time step                             eg. 4,8,10
argv[6] diffusion constant                         eg. 1
argv[7] velocity                                   eg. 1-10
argv[8] endtime                                    eg. 200(debug),10000(real)
*/

#include "./src/directorymake.h"
#include "./src/spaceutility.h"
#include "./src/voro_func.h"

#include <vector>
#include <iostream>
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
  const int kint_tarea = atoi(argv[3]);
  const int kint_smpl = atoi(argv[4]);
  const double ktarea = 0.01*(double)kint_tarea;
  const int ktimestep = atoi(argv[5]);
  const double ktimetics = 1.0/(double)ktimestep;
  const int kint_diffconst = atoi(argv[6]);
  const double kdiffconst = 0.1*(double)(kint_diffconst);
  const int kint_unitvel = atoi(argv[7]);
  const double kunitvel = 0.001*(double)kint_unitvel;
  const int kendtime = atoi(argv[8]);

  //data packing directories
  directoryMake(argv[1]);
  directoryMake(argv[2]);
  string datdirectoryname = argv[1];
  string imagedirectoryname = argv[2];

  //filename settings
  const string kstr_tarea = "ar"+to_string(kint_tarea);
  const string kstr_smpl = "sm"+to_string(kint_smpl);
  const string kstr_diff = "df"+to_string(kint_diffconst);
  const string kstr_vel = "vl"+to_string(kint_unitvel);
  const string kstr_endtime = "et"+to_string(kendtime);
  const string kstr_timestep = "ts"+to_string(ktimestep);
  const string kstr_param = kstr_tarea+kstr_smpl+kstr_diff+kstr_vel+kstr_endtime+kstr_timestep;

  //open msd dat file
  const string kstr_msd = datdirectoryname+"MSD"+kstr_param+".dat";
  ofstream msd;
  msd.open(kstr_msd);

  //open histgram dat file
  const string kstr_hst = datdirectoryname+"HST"+kstr_param+".dat";
  ofstream hst;
  hst.open(kstr_hst);
  vector<double> hstdata(40,0.0);
  
  //initialize random function
  mt.seed(kint_smpl);
  uniform_real_distribution<double> rnd(0.0,1.0);

  //spacetics
  const double kdpos = 0.01;

  //boundary
  Boundary* bd = new Boundary(3.0);
  //internal variables
  vector<VoronoiPoint> vp;
  for( int i = 0; i != particles; ++i) {
    vp.push_back(VoronoiPoint());
    vp[i].r_ = kunitvel;
  }

  
  container con(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
		true,true,true,8);
    
  // Randomly add particles into the container
  setInitialConfiguration(con, vp, bd, mt);

  con.print_custom("%i %v","packing.custom3");

  const int kmsdstarttime = kendtime/20;
  const double khsttics = 1.0/( (double)(particles*(kendtime-kmsdstarttime)) ); 

  for (int time = 0; time != kendtime; ++time) {

     for( int ts = 0; ts != ktimestep; ++ts) {
        cout<<time<<" "<<ts<<endl;
        container eon(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
		  true,true,true,8);
        execLangevinStep(eon, vp, bd, mt, ktarea, kdpos, ktimetics, kdiffconst);
     }

     container don(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
                   true,true,true,8);
     setContainer(don, vp);

     restoreHSTData(don,hstdata,khsttics);

     writeMSDData(time,kmsdstarttime,vp,bd,msd);

     if (time%5 == 0) {
        string stime(to_string(time/5));
        string kstr_p = datdirectoryname+stime+kstr_param+"point.dat";
        string kstr_v = datdirectoryname+stime+kstr_param+"edges.dat";
        ofstream pf;
        pf.open(kstr_p.c_str());
        for (const auto& v : vp) pf<<v.x_<<" "<<v.y_<<" "<<v.z_<<" "<<v.dirx_<<" "<<v.diry_<<" "<<v.dirz_<<endl;
        pf.close();
        //don.draw_particles(kstr_p.c_str());
        don.draw_cells_gnuplot(kstr_v.c_str());

        string kstr_pvj = imagedirectoryname+stime+kstr_param+"pointedges.png";
      
        FILE* gp;
        gp = popen("gnuplot -persist","w");
        fprintf(gp, "set term png size 1200, 1200\n");

        fprintf(gp, "set output \"%s\" \n", kstr_pvj.c_str());
        fprintf(gp, "sp[%f:%f][%f:%f][%f:%f]\"%s\" u 1:2:3:(0.3*$4):(0.3*$5):(0.3*$6) w vector notitle, \"%s\" u 1:2:3 w l lw 0.2 lc rgbcolor \"gray50\" notitle \n", bd->xmin_, bd->xmax_, bd->ymin_, bd->ymax_, bd->zmin_, bd->zmax_, kstr_p.c_str(), kstr_v.c_str());
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
  return EXIT_SUCCESS;
 }
