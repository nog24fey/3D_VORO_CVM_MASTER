#ifndef VORO_FUNC_H_
#define VORO_FUNC_H_ 1

#include "../voro++/voro++.hh"
#include "./voronoipoint.h"
#include "./boundary.h"

#include <random>
#define Pi 3.14159265358979323846
using std::mt19937;
using std::uniform_real_distribution;

using namespace voro;

void setInitialConfiguration(container& con, vector<VoronoiPoint>& vps, Boundary* bdr, mt19937& mt) {
  uniform_real_distribution<double> randm(0.0,1.0);
  
  int index = 0;
  for (auto & v: vps) {
    v.theta_ = randm(mt)*Pi;
    v.phi_ = randm(mt)*2.0*Pi;

    v.x_ = bdr->xmin_+randm(mt)*(bdr->xmax_-bdr->xmin_);
    v.y_ = bdr->ymin_+randm(mt)*(bdr->ymax_-bdr->ymin_);
    v.z_ = bdr->zmin_+randm(mt)*(bdr->zmax_-bdr->zmin_);
    con.put(index,v.x_,v.y_,v.z_);
    ++index;
  }
}

double retTotalEnergy(container& con, double value_tarea) {

  double totalenergy = 0.0;

  c_loop_all cm(con);
  voronoicell_neighbor cl;
  if(cm.start()) do if(con.compute_cell(cl,cm)) {
        double ve = (cl.volume()-1.0);
        double ae = (cl.surface_area()-value_tarea);
        totalenergy += ve*ve+ae*ae;
      } while (cm.inc());

  return totalenergy;
}

void execLangevinStep(container& base_con, vector<VoronoiPoint>& vps, Boundary* bdr, mt19937& mt) {
}

void renewPosition(container& con, vector<VoronoiPoint>& vps, Boundary* bdr) {
  double cx = 0.0; double cy = 0.0; double cz = 0.0;
  
  for (auto & v : vps) {
    v.xo_ = v.x_; v.yo_ = v.y_; v.zo_ = v.z_;
    v.x_ += v.xn_; v.y_ += v.yn_; v.z_ += v.zn_;
    cx += v.xn_; cy += v.yn_; cz += v.zn_;
  }
  double vsize = (double)(vps.size());
  cx /= vsize; cy /= vsize; cz /= vsize;
  int index = 0;
  for (auto & v : vps) {
    v.x_ -= cx; v.y_ -= cy; v.z_ -=cz;
    
    while (v.x_ < bdr->xmin_) v.x_ += bdr->x_axe_leng_;
    while (v.x_ >= bdr->xmax_) v.x_ -= bdr->x_axe_leng_;
    while (v.y_ < bdr->ymin_) v.y_ += bdr->y_axe_leng_;
    while (v.y_ >= bdr->ymax_) v.y_ -= bdr->y_axe_leng_;
    while (v.z_ < bdr->zmin_) v.z_ += bdr->z_axe_leng_;
    while (v.z_ >= bdr->zmax_) v.z_ -= bdr->z_axe_leng_;

    con.put(index, v.x_, v.y_, v.z_);
    ++index;
  }
}
#endif
