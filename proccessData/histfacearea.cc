/*
argv[1] input directory name includeing input data files !!must be ended with "/" !!
argv[2] output directory name !!must be ended with "/" !!
data file prameters as follows.
argv[3] target value of area; timed 0.01           eg. 300-1000
argv[4] sample id; also used for random seed       eg. 0-10
argv[5] (int)time step                             eg. 4,8,10
argv[6] diffusion constant                         eg. 1
argv[7] velocity                                   eg. 1-10
argv[8] endtime                                    eg. 200(debug),10000(real)
*/

#include "../src/directorymake.h"

#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

const int particles=6*6*6;

using namespace std;

int getBin(const int numBinAreas, const double area) {
  if ( area > 2.0) return numBinAreas-1;
  int intarea = (int)( (floor)(area*(double)numBinAreas) );
  /* range 0-199*/
  return intarea/2;
}

int main(int argc, char **argv) {
  //system parameters
  const int kint_tarea = atoi(argv[3]);
  const int kint_smpl = atoi(argv[4]);
  const int ktimestep = atoi(argv[5]);
  const int kint_diffconst = atoi(argv[6]);

  const int kint_unitvel = atoi(argv[7]);
  const int kendtime = atoi(argv[8]);

  //data including directories
  directoryMake(argv[2]);
  string inputdirectoryname = argv[1];
  string outputdirectoryname = argv[2];

  //filename settings
  const string kstr_tarea = "ar"+to_string(kint_tarea);
  const string kstr_smpl = "sm"+to_string(kint_smpl);
  const string kstr_diff = "df"+to_string(kint_diffconst);
  const string kstr_vel = "vl"+to_string(kint_unitvel);
  const string kstr_endtime = "et"+to_string(kendtime);
  const string kstr_timestep = "ts"+to_string(ktimestep);
  const string kstr_param = kstr_tarea+kstr_smpl+kstr_diff+kstr_vel+kstr_endtime+kstr_timestep;

  //open input file
  const string kreadfilename = inputdirectoryname+"HST"+kstr_param+".dat";
  ifstream hst;
  hst.open(kreadfilename);
  //open output file
  string writefilename = outputdirectoryname+"HSTFaceArea"+kstr_param+".dat";
  ofstream facearea;
  facearea.open(writefilename);
  writefilename = outputdirectoryname+"HSTFace13Area"+kstr_param+".dat";
  ofstream face13area;
  face13area.open(writefilename);
  writefilename = outputdirectoryname+"HSTFace14Area"+kstr_param+".dat";
  ofstream face14area;
  face14area.open(writefilename);
  writefilename = outputdirectoryname+"HSTFace15Area"+kstr_param+".dat";
  ofstream face15area;
  face15area.open(writefilename);
  writefilename = outputdirectoryname+"HSTFace16Area"+kstr_param+".dat";
  ofstream face16area;
  face16area.open(writefilename);
  writefilename = outputdirectoryname+"HSTFace17Area"+kstr_param+".dat";
  ofstream face17area;
  face17area.open(writefilename);
  writefilename = outputdirectoryname+"HSTFace18Area"+kstr_param+".dat";
  ofstream face18area;
  face18area.open(writefilename);

  writefilename = outputdirectoryname+"HSTFaces"+kstr_param+".dat";
  ofstream hface;
  hface.open(writefilename);
  


  int igbg;
  double dgbg;
  const int maxNumFaces = 30;
  const int numBinArea = 100; // [0, 0.02], ... ,[1.98, 2.00]
  vector< vector<int> > heat(maxNumFaces, vector<int>( numBinArea, 0));
  vector<int> histface( maxNumFaces, 0);
  int face;
  double area;
  int scalefactor = 0;
  int targetBin;
  int linecount = 0;
  while (!hst.eof()) {
    ++linecount;
    hst>>dgbg;
    hst>>face;
    scalefactor += face;
    ++histface[face];
    hst>>igbg;
    hst>>igbg;
    for ( int ai = 0; ai != face; ++ai) {
      hst>>area;
      targetBin = getBin(numBinArea, area);
      ++heat[face][targetBin];
    }
  }

  for ( int fi = 0; fi != maxNumFaces; ++fi) hface<<fi<<" "<<(double)histface[fi]/(double)linecount<<" "<<kint_tarea<<endl;
  for ( int fi = 0; fi != maxNumFaces; ++fi) {
    for ( int bi = 0; bi != numBinArea; ++bi) {
      facearea<<fi<<" "<<(double)(2*bi)/(double)numBinArea<<" "<<heat[fi][bi]/(double)scalefactor<<endl;
      if ( fi == 13) face13area<<fi<<" "<<(double)(2*bi)/(double)numBinArea<<" "<<(double)heat[fi][bi]/(double)scalefactor<<endl;
      if ( fi == 14) face14area<<fi<<" "<<(double)(2*bi)/(double)numBinArea<<" "<<(double)heat[fi][bi]/(double)scalefactor<<endl;
      if ( fi == 15) face15area<<fi<<" "<<(double)(2*bi)/(double)numBinArea<<" "<<(double)heat[fi][bi]/(double)scalefactor<<endl;
      if ( fi == 16) face16area<<fi<<" "<<(double)(2*bi)/(double)numBinArea<<" "<<(double)heat[fi][bi]/(double)scalefactor<<endl;      
      if ( fi == 17) face17area<<fi<<" "<<(double)(2*bi)/(double)numBinArea<<" "<<(double)heat[fi][bi]/(double)scalefactor<<endl;
      if ( fi == 18) face18area<<fi<<" "<<(double)(2*bi)/(double)numBinArea<<" "<<(double)heat[fi][bi]/(double)scalefactor<<endl;
    }
    facearea<<endl;
  }

  hst.close();
  facearea.close();
  face13area.close();
  face14area.close();
  face15area.close();
  face16area.close();
  face17area.close();
  face18area.close();
  hface.close();
    
  return EXIT_SUCCESS;
}
