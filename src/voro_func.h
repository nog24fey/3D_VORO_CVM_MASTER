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

void setPeriodicity(double &x, double &y, double &z, Boundary* bdr) {
  while (x < bdr->xmin_) x += bdr->x_axe_leng_;
  while (x >= bdr->xmax_) x -= bdr->x_axe_leng_;
  while (y < bdr->ymin_) y += bdr->y_axe_leng_;
  while (y >= bdr->ymax_) y -= bdr->y_axe_leng_;
  while (z < bdr->zmin_) z += bdr->z_axe_leng_;
  while (z >= bdr->zmax_) z -= bdr->z_axe_leng_;
}

void scaleUnitVector(VoronoiPoint& vps) {
  scaleUnitVector(vps.dirx_, vps.diry_, vps.dirz_);
}

void setContainer(container& con, vector<VoronoiPoint>& vps) {
   int index = 0;
   for (auto & v: vps){
      con.put(index,v.x_,v.y_,v.z_);
      ++index;
   }
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

void renewPositions(vector<VoronoiPoint>& vps) {
  //renewPosition
  double cx = 0.0; double cy = 0.0; double cz = 0.0;

  for (auto & v : vps) {
    v.xo_ = v.x_; v.yo_ = v.y_; v.zo_ = v.z_;
    double px = 0.5*(v.xn1_+v.xn2_);
    doublw py = 0.5*(v.yn1_+v.yn2_);
    double pz = 0.5*(v.zn1_+v.zn2_);
    v.x_ += px; v.y_ += py; v.z_ += pz; 
    cx += px; cy += py; cz += pz;
  }
  double vsize = (double)(vps.size());
  cx /= vsize; cy /= vsize; cz /= vsize;
  int index = 0;
  for (auto & v : vps) {
    v.x_ -= cx; v.y_ -= cy; v.z_ -=cz;
    setPeriodicity(v.x_,v.y_,v.z_,bdr);
    ++index;
  }
}

void execLangevinStep(container& base_con, vector<VoronoiPoint>& vps, Boundary* bdr, mt19937& mt, double value_tarea, double dps, double delta_t) {
  //prepare constants for speed
  double xmin = bdr->xmin_; double xmax = bdr->xmax_;
  double ymin = bdr->ymin_; double ymax = bdr->ymax_;
  double zmin = bdr->zmin_; double zmax = bdr->zmax_;
  int nx = bdr->nx_; int ny = bdr->ny_; int nz = bdr->nz_;

  double rscale = sqrt(2.0*D*delta_t/3.0);
  
  //uniform_real_distribution<double> randm(0.0,1.0);
  normal_distribution<double> nml(0.0,1.0);
  double D = 1;//あとでパラメータにする.

  int index = 0;
  //(random+selfpropel)motion
  for (auto & v : vps) {

    v.nsx_ = nml(mt)*rscale; v.nsy_ = nml(mt)*rscale; v.nsz_ = nml(mt)*rscale;
    v.dirx_ += v.nsy_*v.dirz_ - v.nsz_*v.diry_;
    v.diry_ += v.nsz_*v.dirx_ - v.nsx_*v.dirz_;
    v.dirz_ += v.nsx_*v.diry_ - v.nsy_*v.dirx_;
    scaleUnitVector(v);
    
    base_con.put(index,v.x_,v.y_,v.z_);
    ++index;
  }
    
  double totalenergy_base = retTotalEnergy(base_con, value_tarea);

  double ox, oy, oz, vx, vy, vz;
  for (int rktime = 0; rktime !=2; ++rktime) { //calc F1 in rktime == 0 and F2 in rktime == 1;

    int index_self = 0;

    for (auto & v : vps) {

      if (rktime == 0) {
        ox = v.x_; oy = v.y_; oz = v.z_;
      }else {
        ox = v.x_+delta_t*(v.xn1_); oy = v.y_+delta_t*(v.yn1_); oz = v.z_+delta_t*(v.zn1_);
        setPeriodicity(ox, oy, oz, bdr);
        
      }
      vx = ox+dps; vy = oy+dps; vz = oz+dps;
      setPeriodicity(vx, vy, vz, bdr);

      const double kspl_x = v.r_*v.dirx_;
      const double kspl_y = v.r_*v.diry_;
      const double kspl_z = v.r_*v.dirz_;      
      
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
      if (rktime == 0) {
        v.xn1_ = -0.01*(tot_eng_x-totalenergy_base)/dps + kspl_x;
      } else {
        v.xn2_ = -0.01*(tot_eng_x-totalenergy_base)/dps + kspl_x;
      }
      
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
      if (rktime == 0) {
        v.yn1_ = -0.01*(tot_eng_y-totalenergy_base)/dps + kspl_y;
      } else {
        v.yn2_ = -0.01*(tot_eng_y-totalenergy_base)/dps + kspl_y;
      }
      
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
      if (rktime == 0) {
        v.zn1_ = -0.01*(tot_eng_z-totalenergy_base)/dps + kspl_z;
      } else {
        v.zn1_ = -0.01*(tot_eng_z-totalenergy_base)/dps + kspl_z;
      }
      
      ++index_self;
    }
  }

  renewPositions(vps);

}
  
#endif
