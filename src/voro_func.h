#ifndef VORO_FUNC_H_
#define VORO_FUNC_H_ 1

#include "../voro++/voro++.hh"
#include "./voronoipoint.h"
#include "./boundary.h"

#include <random>
#define Pi 3.14159265358979323846
using std::mt19937;
using std::uniform_real_distribution;
using std::normal_distribution;

using namespace voro;

void scaleUnitVector(VoronoiPoint& vps) {
  scaleUnitVector(vps.dirx_, vps.diry_, vps.dirz_);
}

void setInitialConfiguration(container& con, vector<VoronoiPoint>& vps, Boundary* bdr, mt19937& mt) {
  normal_distribution<double> randm(0.0,1.0);
  
  int index = 0;
  for (auto & v: vps) {
    v.dirx_ = randm(mt); v.diry_ = randm(mt); v.dirz_ = randm(mt);
    scaleUnitVector(v);

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

void execLangevinStep(container& base_con, vector<VoronoiPoint>& vps, Boundary* bdr, mt19937& mt, double value_tarea, double dps) {
  //prepare constants for speed
  double xmin = bdr->xmin_; double xmax = bdr->xmax_;
  double ymin = bdr->ymin_; double ymax = bdr->ymax_;
  double zmin = bdr->zmin_; double zmax = bdr->zmax_;
  int nx = bdr->nx_; int ny = bdr->ny_; int nz = bdr->nz_;
        
  //uniform_real_distribution<double> randm(0.0,1.0);
  normal_distribution<double> nml(0.0,1.0);
  double D = 1;//あとでパラメータにする.

  int index = 0;
  //(random+selfpropel)motion
  for (auto & v : vps) {
    double rscale = sqrt(2*D/3);
    double rx = nml(mt)*rscale; double ry = nml(mt)*rscale; double rz = nml(mt)*rscale;
    v.dirx_ += ry*v.dirz_ - rz*v.diry_;
    v.diry_ += rz*v.dirx_ - rx*v.dirz_;
    v.dirz_ += rx*v.diry_ - ry*v.dirx_;
    scaleUnitVector(v);
    
    v.x_ += v.r_*v.dirx_;
    v.y_ += v.r_*v.diry_;
    v.z_ += v.r_*v.dirz_;
    
    base_con.put(index,v.x_,v.y_,v.z_);
    ++index;
  }
    
  double totalenergy_base = retTotalEnergy(base_con, value_tarea);

  double ox, oy, oz, vx, vy, vz;
  int index_self = 0;
  for (auto & v : vps) {
    ox = v.x_;
    oy = v.y_;
    oz = v.z_;
    vx = v.x_+dps;
    vy = v.y_+dps;
    vz = v.z_+dps;

      container xon(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,true,true,true,8);
      int index_x = 0;
      for (auto & s : vps) {
	if (v!=s) {
	  xon.put(index_x,s.x_,s.y_,s.z_);
	}
        ++index_x;
      }
      xon.put(index_self,vx,oy,oz);
      double tot_eng_x = retTotalEnergy(xon, value_tarea);
      v.xn_ = -0.01*(tot_eng_x-totalenergy_base)/dps;// + v.r_*v.dirx_;;

      
      container yon(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,true,true,true,8);
      int index_y = 0;
      for (auto & s : vps) {
	if (v!=s) {
	  yon.put(index_y,s.x_,s.y_,s.z_);
	}
        ++index_y;
      }      
      yon.put(index_self,ox,vy,oz);
      double tot_eng_y = retTotalEnergy(yon, value_tarea);
      v.yn_ = -0.01*(tot_eng_y-totalenergy_base)/dps;// + v.r_*v.diry_;;

      
      container zon(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,true,true,true,8);
      int index_z = 0;
      for (auto & s : vps) {
	if (v!=s) {
	  zon.put(index_z,s.x_,s.y_,s.z_);
	}
        ++index_z;
      }
      zon.put(index_self,ox,oy,vz);
      double tot_eng_z = retTotalEnergy(zon, value_tarea);
      v.zn_ = -0.01*(tot_eng_z-totalenergy_base)/dps;// + v.r_*v.dirz_;;

      ++index_self;
  }

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
