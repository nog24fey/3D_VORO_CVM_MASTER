/*
argv[1] directory name packing data files !!must be ended with "/" !!
argv[2] directory name packing image files !!must be ended with "/" !!
argv[3] target value of area; timed 0.01           eg. 3000-10000
argv[4] sample id; also used for random seed       eg. 0-10
argv[5] force shorten MD steps                     eg. 100-1000
argv[6] initial relaxation MD steps                eg. 10000(real)
argv[7] file numbering for read file               eg, 20 , 50, 2000
*/

#include "./src/directorymake.h"
#include "./src/spaceutility.h"
#include "./src/lvoro_func.h"

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
  const double ktarea = 0.001*(double)kint_tarea;
  const int kfsMDsteps = atoi(argv[5]);
  const int kiniMDsteps = atoi(argv[6]);
  const int kreadtime = atoi(argv[7]);
  
  //data packing directories
  directoryMake(argv[1]);
  directoryMake(argv[2]);
  string datdirectoryname = argv[1];
  string imagedirectoryname = argv[2];

  //filename settings
  const string kstr_tarea = "ar"+to_string(kint_tarea);
  const string kstr_smpl = "sm"+to_string(kint_smpl);
  const string kstr_iniMDsteps = "is"+to_string(kiniMDsteps);
  const string kstr_fsMDsteps = "fs"+to_string(kfsMDsteps);
  const string kstr_readtime = "in"+to_string(kreadtime);
  const string kstr_param = kstr_tarea+kstr_smpl+kstr_iniMDsteps+kstr_fsMDsteps+kstr_readtime;

  //open histgram dat file
  //const string kstr_hst = datdirectoryname+"MCHST"+kstr_param+".dat";
  //ofstream hst;
  //hst.open(kstr_hst);
  //open totalenergy dat file
  const string kstr_eng = datdirectoryname+"MCEng"+kstr_param+".dat";
  ofstream eng;
  eng.open(kstr_eng);
  
  //initialize random function
  mt.seed(kint_smpl);
  uniform_real_distribution<double> rnd(0.0,1.0);

  //spacetics
  const double kdpos = 0.01;

  //boundary
  Boundary* bd = new Boundary(3.0);
  //internal variables
  vector<LVoronoiPoint> vp;
  for( int i = 0; i != particles; ++i) {
    vp.push_back(LVoronoiPoint());
  }

  container con(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
		true,true,true,8);    
  // Randomly add particles into the container
  readInitialConfiguration(con, vp, bd, mt);
  double totalenergy = retTotalEnergy(don, ktarea);
  eng<<0<<" "<<totalenergy<<endl;
  for (int time = 0; time != kiniMDsteps; ++time) {

    container eon(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
                  true,true,true,8);
    execMonteCarloStep(eon, vp, bd, mt, ktarea, kdpos);

    container don(bd->xmin_,bd->xmax_,bd->ymin_,bd->ymax_,bd->zmin_,bd->zmax_,bd->nx_,bd->ny_,bd->nz_,
                  true,true,true,8);
    setContainer(don, vp);

    totalnergy = retTotalEnergy(don, ktarea);
    eng<<time+1<<" "<<totalenergy<<endl;
    //writeHSTData(time, 5, don,hst);
    writeSnapShotFile(time+1, 5, don, vp, bd, datdirectoryname, imagedirectoryname, kstr_param);
  }

  //hst.close();
  eng.close();

  delete bd;
  return EXIT_SUCCESS;
 }
