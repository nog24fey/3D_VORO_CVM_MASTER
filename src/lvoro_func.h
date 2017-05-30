#ifndef VORO_FUNC_H_
#define VORO_FUNC_H_ 1

#include "../voro++/voro++.hh"
#include "./lvoronoipoint.h"
#include "./boundary.h"

#include <omp.h>
#include <random>
#include <fstream>

#define Pi 3.14159265358979323846
using std::mt19937;
using std::uniform_real_distribution;
using std::normal_distribution;

using namespace voro;

void setPeriodicity(double &x, double &y, double &z, const Boundary* bdr) {
  while (x < bdr->xmin_) x += bdr->x_axe_leng_;
  while (x >= bdr->xmax_) x -= bdr->x_axe_leng_;
  while (y < bdr->ymin_) y += bdr->y_axe_leng_;
  while (y >= bdr->ymax_) y -= bdr->y_axe_leng_;
  while (z < bdr->zmin_) z += bdr->z_axe_leng_;
  while (z >= bdr->zmax_) z -= bdr->z_axe_leng_;
}


void setContainer(container& con, vector<LVoronoiPoint>& vps) {
   int index = 0;
   for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v){
      con.put(index,(*v).x_,(*v).y_,(*v).z_);
      ++index;
   }
}

void setInitialConfiguration(container& con, vector<LVoronoiPoint>& vps, const Boundary* bdr, mt19937& mt) {

  uniform_real_distribution<double> randm(-1.0,1.0);
  
  int index = 0;
  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {

    (*v).x_ = 0.5*randm(mt)*(bdr->xmax_-bdr->xmin_);
    (*v).y_ = 0.5*randm(mt)*(bdr->ymax_-bdr->ymin_);
    (*v).z_ = 0.5*randm(mt)*(bdr->zmax_-bdr->zmin_);
    con.put(index,(*v).x_,(*v).y_,(*v).z_);
    ++index;
  }
}

void readInitialConfiguration(container& con, vector<LVoronoiPoint>& vps, const Boundary* bdr, mt19937& mt, ifstream ifs) {

  int index = 0;
  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    ifs>>(*v).x_; ifs>>(*v).y_; ifs>>(*v).z_;
    con.put(index,(*v).x_,(*v).y_,(*v).z_);
    ++index;
  }
}

double retTotalEnergy(container& con, const double value_tarea, const double value_eratio=1.0) {

  double totalenergy = 0.0;

  c_loop_all cm(con);
  voronoicell_neighbor cl;
  if(cm.start()) do if(con.compute_cell(cl,cm)) {
        double ve = (cl.volume()-1.0);
        double ae = (cl.surface_area()-value_tarea);
        totalenergy += ve*ve+ae*ae/value_eratio;
      } while (cm.inc());

  return totalenergy;
}

void renewPositions(vector<LVoronoiPoint>& vps, const Boundary* bdr) {

  double cx = 0.0; double cy = 0.0; double cz = 0.0;

  const double kdvsize = (double)(vps.size());

  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    const double px = (*v).xn2_;
    const double py = (*v).yn2_;
    const double pz = (*v).zn2_;
    (*v).x_ += px; (*v).y_ += py; (*v).z_ += pz; 
    cx += px; cy += py; cz += pz;
  }
 
  cx /= kdvsize; cy /= kdvsize; cz /= kdvsize;

  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    (*v).x_ -= cx; (*v).y_ -= cy; (*v).z_ -= cz;
    setPeriodicity((*v).x_,(*v).y_,(*v).z_,bdr);
  }

}

void execMonteCarloStep(container& base_con, vector<LVoronoiPoint>& vps, const Boundary* bdr, mt19937& mt, const double value_tarea, const double value_eratio, const double dps) {
  //prepare constants for speed
  const double xmin = bdr->xmin_; const double xmax = bdr->xmax_;
  const double ymin = bdr->ymin_; const double ymax = bdr->ymax_;
  const double zmin = bdr->zmin_; const double zmax = bdr->zmax_;
  const int nx = bdr->nx_; const int ny = bdr->ny_; const int nz = bdr->nz_;

  double alpha;
  if ( value_eratio >= 10.0) {
    alpha = 0.1;
  } else if ( value_eratio >= 1.0) {
    alpha = 0.01;
  } else if ( value_eratio >= 0.1) {
    alpha = 0.001;
  } else {
    alpha = 0.0001;
  }

  normal_distribution<double> nml(0.0,1.0);

  int index = 0;
  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    base_con.put(index,(*v).x_,(*v).y_,(*v).z_);
    ++index;
  }
  const double totalenergy_base = retTotalEnergy(base_con, value_tarea, value_eratio);

  int index_self = 0;
  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) { 

    (*v).xn1_ = (*v).x_+dps; (*v).yn1_ = (*v).y_+dps; (*v).zn1_ = (*v).z_+dps;
    setPeriodicity((*v).xn1_, (*v).yn1_, (*v).zn1_, bdr);

    container xon(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,true,true,true,8);
    int index_x = 0;
    for (std::vector<LVoronoiPoint>::iterator s = vps.begin(); s != vps.end(); ++s) {
      if ((*v)!=(*s)) {
        xon.put(index_x,(*s).x_,(*s).y_,(*s).z_);
      }
      ++index_x;
    }
    xon.put(index_self,(*v).xn1_,(*v).y_,(*v).z_);
    const double tot_eng_x = retTotalEnergy(xon, value_tarea, value_eratio);
    (*v).xn2_ = -alpha*(tot_eng_x-totalenergy_base)/dps ;

    container yon(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,true,true,true,8);
    int index_y = 0;
    for (std::vector<LVoronoiPoint>::iterator s = vps.begin(); s != vps.end(); ++s) {
      if ((*v)!=(*s)) {
        yon.put(index_y,(*s).x_,(*s).y_,(*s).z_);
      }
      ++index_y;
    }      
    yon.put(index_self,(*v).x_,(*v).yn1_,(*v).z_);
    const double tot_eng_y = retTotalEnergy(yon, value_tarea, value_eratio);
    (*v).yn2_ = -alpha*(tot_eng_y-totalenergy_base)/dps;

    container zon(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,true,true,true,8);
    int index_z = 0;
    for (std::vector<LVoronoiPoint>::iterator s = vps.begin(); s != vps.end(); ++s) {
      if ((*v)  !=(*s)) {
        zon.put(index_z,(*s).x_,(*s).y_,(*s).z_);
      }
      ++index_z;
    }
    zon.put(index_self,(*v).x_,(*v).y_,(*v).zn1_);
    const double tot_eng_z = retTotalEnergy(zon, value_tarea, value_eratio);
    (*v).zn2_ = -alpha*(tot_eng_z-totalenergy_base)/dps;

    ++index_self;
  }

  renewPositions(vps, bdr);

}

void writeHSTData(const int time, const int period, container& con, ofstream& ofsHST) {
  if (time%period != 0) return;
  c_loop_all cmh(con);
  voronoicell_neighbor ch;
  if (cmh.start()) do if (con.compute_cell(ch,cmh)) {
        ofsHST<<ch.volume()<<" "<<ch.number_of_faces()<<" "
              <<ch.number_of_edges()<<" "<<ch.p<<" ";
        vector<double> areas;
        ch.face_areas(areas);
        for (auto & ar : areas) ofsHST<<ar<<" ";
        ofsHST<<endl;
      } while (cmh.inc());
  ofsHST<<endl;
}

void writeSnapShotFile(const int time, const int period, container& cont, const vector<LVoronoiPoint>& vps, const Boundary* bdr, const string datdirectoryname, const string imgdirectoryname, const string parasetnames) {
  if (time%period != 0) return;

  const string stime(to_string(time/5));
  const string kstr_p = datdirectoryname+stime+parasetnames+"point.dat";
  const string kstr_v = datdirectoryname+stime+parasetnames+"edges.dat";

  ofstream pf;
  pf.open(kstr_p.c_str());
  for (const auto& v : vps) pf<<v.x_<<" "<<v.y_<<" "<<v.z_<<endl;
  pf.close();
  cont.draw_cells_gnuplot(kstr_v.c_str());

  const string kstr_pvj = imgdirectoryname+stime+parasetnames+"pointedges.png";

  FILE* gp;

  gp = popen("gnuplot -persist","w");
  fprintf(gp, "set term png size 1200, 1200\n");

  fprintf(gp, "set output \"%s\" \n", kstr_pvj.c_str());
  fprintf(gp, "sp[%f:%f][%f:%f][%f:%f]\"%s\" u 1:2:3 w p ps 2 pt 7 lc rgb \"purple\" notitle, \"%s\" u 1:2:3 w l lw 1 lc rgbcolor \"gray70\" notitle \n", bdr->xmin_, bdr->xmax_, bdr->ymin_, bdr->ymax_, bdr->zmin_, bdr->zmax_, kstr_p.c_str(), kstr_v.c_str());
  fprintf(gp, "set output\n");
  pclose(gp);

}

#endif
