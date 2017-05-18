#ifndef VORO_FUNC_H_
#define VORO_FUNC_H_ 1

#include "../voro++/voro++.hh"
#include "./voronoipoint.h"
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

void scaleUnitVector(VoronoiPoint& vps) {
  scaleUnitVector(vps.dirx_, vps.diry_, vps.dirz_);
}

void setContainer(container& con, vector<VoronoiPoint>& vps) {
   int index = 0;
   for (std::vector<VoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v){
      con.put(index,(*v).x_,(*v).y_,(*v).z_);
      ++index;
   }
}
void setInitialConfiguration(container& con, vector<VoronoiPoint>& vps, const Boundary* bdr, mt19937& mt) {

  uniform_real_distribution<double> randm(-1.0,1.0);
  
  int index = 0;
  for (std::vector<VoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    (*v).dirx_ = randm(mt); (*v).diry_ = randm(mt); (*v).dirz_ = randm(mt);
    scaleUnitVector((*v));

    (*v).x_ = 0.5*randm(mt)*(bdr->xmax_-bdr->xmin_);
    (*v).y_ = 0.5*randm(mt)*(bdr->ymax_-bdr->ymin_);
    (*v).z_ = 0.5*randm(mt)*(bdr->zmax_-bdr->zmin_);
    con.put(index,(*v).x_,(*v).y_,(*v).z_);
    ++index;
  }
}

double retTotalEnergy(container& con, const double value_tarea) {

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

void renewPositions(vector<VoronoiPoint>& vps, const Boundary* bdr, const double delta_t) {

  double cx = 0.0; double cy = 0.0; double cz = 0.0;
  const double dt = delta_t;

  const double kdvsize = (double)(vps.size());

  for (std::vector<VoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    (*v).xo_ = (*v).x_; (*v).yo_ = (*v).y_; (*v).zo_ = (*v).z_;
    const double px = 0.5*dt*((*v).xn1_+(*v).xn2_);
    const double py = 0.5*dt*((*v).yn1_+(*v).yn2_);
    const double pz = 0.5*dt*((*v).zn1_+(*v).zn2_);
    (*v).x_ += px; (*v).y_ += py; (*v).z_ += pz; 
    cx += px; cy += py; cz += pz;
  }
 
  cx /= kdvsize; cy /= kdvsize; cz /= kdvsize;

  for (std::vector<VoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    (*v).x_ -= cx; (*v).y_ -= cy; (*v).z_ -= cz;
    setPeriodicity((*v).x_,(*v).y_,(*v).z_,bdr);
  }

}

void execLangevinStep(container& base_con, vector<VoronoiPoint>& vps, const Boundary* bdr, mt19937& mt, const double value_tarea, const double dps, const double delta_t, const double diffusion_constant) {
  //prepare constants for speed
  const double xmin = bdr->xmin_; const double xmax = bdr->xmax_;
  const double ymin = bdr->ymin_; const double ymax = bdr->ymax_;
  const double zmin = bdr->zmin_; const double zmax = bdr->zmax_;
  const int nx = bdr->nx_; const int ny = bdr->ny_; const int nz = bdr->nz_;

  const double D = diffusion_constant;
  const double alpha = 0.01;
  const double rscale = sqrt(2.0*D*delta_t/3.0);
  
  normal_distribution<double> nml(0.0,1.0);

  int index = 0;

  for (std::vector<VoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    (*v).nsx_ = nml(mt)*rscale; (*v).nsy_ = nml(mt)*rscale; (*v).nsz_ = nml(mt)*rscale;
    (*v).dirx_ += (*v).nsy_*(*v).dirz_ - (*v).nsz_*(*v).diry_;
    (*v).diry_ += (*v).nsz_*(*v).dirx_ - (*v).nsx_*(*v).dirz_;
    (*v).dirz_ += (*v).nsx_*(*v).diry_ - (*v).nsy_*(*v).dirx_;
    scaleUnitVector((*v));
    
    base_con.put(index,(*v).x_,(*v).y_,(*v).z_);
    ++index;
  }
    
  const double totalenergy_base = retTotalEnergy(base_con, value_tarea);

  for (int rktime = 0; rktime !=2; ++rktime) { //calc F1 in rktime == 0 and F2 in rktime == 1;

    int index_self = 0;

    for (std::vector<VoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) { 

      if (rktime == 0) {
        (*v).tpox_ = (*v).x_; (*v).tpoy_ = (*v).y_; (*v).tpoz_ = (*v).z_;
      }else {
        (*v).tpox_ = (*v).x_+delta_t*((*v).xn1_); (*v).tpoy_ = (*v).y_+delta_t*((*v).yn1_); (*v).tpoz_ = (*v).z_+delta_t*((*v).zn1_);
        setPeriodicity((*v).tpox_, (*v).tpoy_, (*v).tpoz_, bdr);   
      }
      (*v).tpnx_ = (*v).tpox_+dps; (*v).tpny_ = (*v).tpoy_+dps; (*v).tpnz_ = (*v).tpoz_+dps;
      setPeriodicity((*v).tpnx_, (*v).tpny_, (*v).tpnz_, bdr);

      const double kspl_x = (*v).r_*(*v).dirx_;
      const double kspl_y = (*v).r_*(*v).diry_;
      const double kspl_z = (*v).r_*(*v).dirz_;      

      container xon(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,true,true,true,8);
      int index_x = 0;
      for (std::vector<VoronoiPoint>::iterator s = vps.begin(); s != vps.end(); ++s) {
        if ((*v)!=(*s)) {
          xon.put(index_x,(*s).x_,(*s).y_,(*s).z_);
        }
        ++index_x;
      }
      xon.put(index_self,(*v).tpnx_,(*v).tpoy_,(*v).tpoz_);
      const double tot_eng_x = retTotalEnergy(xon, value_tarea);
      if (rktime == 0) {
        (*v).xn1_ = -alpha*(tot_eng_x-totalenergy_base)/dps + kspl_x;
      } else {
        (*v).xn2_ = -alpha*(tot_eng_x-totalenergy_base)/dps + kspl_x;
      }

      container yon(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,true,true,true,8);
      int index_y = 0;
      for (std::vector<VoronoiPoint>::iterator s = vps.begin(); s != vps.end(); ++s) {
        if ((*v)!=(*s)) {
          yon.put(index_y,(*s).x_,(*s).y_,(*s).z_);
        }
        ++index_y;
      }      
      yon.put(index_self,(*v).tpox_,(*v).tpny_,(*v).tpoz_);
      const double tot_eng_y = retTotalEnergy(yon, value_tarea);
      if (rktime==0) {
        (*v).yn1_ = -alpha*(tot_eng_y-totalenergy_base)/dps + kspl_y;
      } else {
        (*v).yn2_ = -alpha*(tot_eng_y-totalenergy_base)/dps + kspl_y;
      }

      container zon(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,true,true,true,8);
      int index_z = 0;
      //int vjz = 0;
      for (std::vector<VoronoiPoint>::iterator s = vps.begin(); s != vps.end(); ++s) {
        if ((*v)!=(*s)) {
          zon.put(index_z,(*s).x_,(*s).y_,(*s).z_);
        }
        ++index_z;
      }
      zon.put(index_self,(*v).tpox_,(*v).tpoy_,(*v).tpnz_);
      const double tot_eng_z = retTotalEnergy(zon, value_tarea);
      if (rktime == 0) {
        (*v).zn1_ = -alpha*(tot_eng_z-totalenergy_base)/dps + kspl_z;
      } else {
        (*v).zn2_ = -alpha*(tot_eng_z-totalenergy_base)/dps + kspl_z;
      }
      
      ++index_self;
    }
  }

  renewPositions(vps, bdr, delta_t);

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

void writeMSDData(const int time, const int starttime, vector<VoronoiPoint>& vps, const Boundary* bdr, ofstream& ofsMSD, bool finishflag) {
  if (time > starttime) {

    ofsMSD<<time-starttime<<" ";
    double meandiff = 1.0;
    //double meandiff = 0.0;
    double size = (double)(vps.size());
    for (const auto& v : vps) {
      meandiff *=pow(squaredDistanceForMSD(v.x_,v.y_,v.z_,v.x0_,v.y0_,v.z0_,bdr->x_axe_leng_,bdr->y_axe_leng_,bdr->z_axe_leng_),1.0/size);
      //meandif += squaredDistanceForMSD(v.x_,v.y_,v.x0_,v.y0_,bdr->x_axe_leng_,bdr->y_axe_leng_);
    }
    //meandiff /= size;
    ofsMSD<<meandiff<<endl;
    if ( meandiff > bdr->x_axe_leng_/4.0 ) finishflag = true;
    
  } else if (time == starttime) {
    
    for (auto& v : vps) {
      v.x0_ = v.x_; v.y0_ = v.y_;; v.z0_ = v.z_;
    }
    
  } else {
    return;
  }
  
}

void writeSnapShotFile(const int time, const int period, container& cont, const vector<VoronoiPoint>& vps, const Boundary* bdr, const string datdirectoryname, const string imgdirectoryname, const string parasetnames) {
  if (time%period != 0) return;

  const string stime(to_string(time/5));
  const string kstr_p = datdirectoryname+stime+parasetnames+"point.dat";
  const string kstr_v = datdirectoryname+stime+parasetnames+"edges.dat";

  ofstream pf;
  pf.open(kstr_p.c_str());
  for (const auto& v : vps) pf<<v.x_<<" "<<v.y_<<" "<<v.z_<<" "<<v.dirx_<<" "<<v.diry_<<" "<<v.dirz_<<endl;
  pf.close();
  cont.draw_cells_gnuplot(kstr_v.c_str());

  const string kstr_pvj = imgdirectoryname+stime+parasetnames+"pointedges.png";

  FILE* gp;

  gp = popen("gnuplot -persist","w");
  fprintf(gp, "set term png size 1200, 1200\n");

  fprintf(gp, "set output \"%s\" \n", kstr_pvj.c_str());
  fprintf(gp, "sp[%f:%f][%f:%f][%f:%f]\"%s\" u 1:2:3:(0.3*$4):(0.3*$5):(0.3*$6) w vector lw 2 lc rgb \"purple\" notitle, \"%s\" u 1:2:3 w l lw 1 lc rgbcolor \"gray70\" notitle \n", bdr->xmin_, bdr->xmax_, bdr->ymin_, bdr->ymax_, bdr->zmin_, bdr->zmax_, kstr_p.c_str(), kstr_v.c_str());
  fprintf(gp, "set output\n");
  pclose(gp);

}

#endif
