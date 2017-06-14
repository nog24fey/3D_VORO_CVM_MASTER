#ifndef LVORO_FUNC_2_H_
#define LVORO_FUNC_2_H_ 1

#include "./lvoro_func.h"

#include <iostream>
#include <iomanip>
using std::string;
using std::cout;

using namespace voro;

void readData(string s, vector<LVoronoiPoint>& vps, container& con, const int numparticles) {
  ifstream ifs;
  ifs.open(s.c_str());

  for (int i = 0; i != numparticles; ++i) {
    double x, y, z;
    ifs >> x; ifs >> y; ifs >> z;
    vps[i].x_ = x;
    vps[i].y_ = y;
    vps[i].z_ = z;
    con.put(i, x, y, z);
  }

  ifs.close();

}

void linkVPtoContainer( vector<LVoronoiPoint>& vps, container& con) {
  int index = 0;
  for ( auto &v : vps) {
    con.put(index, v.x_, v.y_, v.z_);
    ++index;
  }
}

void copyVpstoVps( vector<LVoronoiPoint>& vpsorigin, vector<LVoronoiPoint>& vps) {
  for ( const auto & vo : vpsorigin) {
    vps.push_back(vo);
  }
}

/*
 *list = i*numparticles + j
 */
void makeTransitionPairs( container& con, vector<int>& list, const int numparticles) {
  const int listsize = 100;

  c_loop_all cmh(con);
  voronoicell_neighbor ch;

  bool breaktag = false;

  if (cmh.start()) do if (con.compute_cell(ch, cmh)) {
        int id = cmh.pid();
        vector<int> neighbors;
        ch.neighbors(neighbors);

        for (auto ni : neighbors) {
          if ( ni > id) list.push_back(id*numparticles+ni);
          if ( list.size() >= listsize ) {
            breaktag = true;
            break;
          }
        }

        if (breaktag) break;
      } while (cmh.inc());

  assert( list.size() == listsize);
}

bool checkPairness( container& con, const int pair1, const int pair2) {
  assert( pair1 < pair2);

  c_loop_all cmh(con);
  voronoicell_neighbor ch;

  bool breaktag = false;

  if (cmh.start()) do if (con.compute_cell(ch, cmh)) {
        int id = cmh.pid();
        if ( id != pair1 ) continue;

        vector<int> neighbors;
        ch.neighbors(neighbors);
        for (auto ni : neighbors) {
          if ( ni == pair2 ) return true;
        }

        return false;

      } while (cmh.inc());
}

/*
 * ret = (x,y,z,rr) rr is distance of the P_pair1 and P_pair2
 */
void getDirectionVectorBetweenPairs( const int pair1, const int pair2, vector<LVoronoiPoint>& vps, const Boundary* bdr, vector<double>& ret) {
   assert( pair1 < pair2);
   LVoronoiPoint _V1 = vps[pair1];
   LVoronoiPoint _V2 = vps[pair2];

   const double _xax = bdr->x_axe_leng_;
   const double _yax = bdr->y_axe_leng_;
   const double _zax = bdr->z_axe_leng_;
   getDirectionVector(_V1.x_, _V1.y_, _V1.z_, _V2.x_, _V2.y_, _V2.z_, _xax, _yax, _zax, ret);
}

void awayPairsSlowly( const int pair1, const int pair2, vector<LVoronoiPoint>&vps, const Boundary* bdr, const double factor) {
  /* factor(eg. 0.1) make  pair distance 1.1 times longer
   *
   */
  vector<double> ret;
  getDirectionVectorBetweenPairs(pair1, pair2, vps, bdr, ret);
  vps[pair1].x_ -= 0.5*factor*ret[0];
  vps[pair1].y_ -= 0.5*factor*ret[1];
  vps[pair1].z_ -= 0.5*factor*ret[2];
  vps[pair2].x_ += 0.5*factor*ret[0];
  vps[pair2].y_ += 0.5*factor*ret[1];
  vps[pair2].z_ += 0.5*factor*ret[2];
}

/*
 * ret = area, true -> have a contact face.
 */
bool getAreaBetweenPairs(container& con, const int pair1, const int pair2, double& ret) {

   assert( pair1 < pair2);

   vector<double> area1, area2;

   c_loop_all cmh(con);
   voronoicell_neighbor ch;
   int breaktag = 0;
   if (cmh.start()) do if (con.compute_cell(ch, cmh)) {
            int id = cmh.pid();
            if ( (id != pair1) && (id != pair2)  ) continue;
            if ( id == pair1) {
               ch.face_areas(area1);
               ++breaktag;
            } else {
               assert( id == pair2);
               ch.face_areas(area2);
               ++breaktag;
            }

            if (breaktag == 2) break;

         } while (cmh.inc());
   assert( area1.size() != 0);
   assert( area2.size() != 0);

   for ( const auto a1 : area1) {
      for ( const auto a2 : area2) {
        if (  fabs(a1 - a2) < 1e-7 ) {
            ret = a1;
            return true;
         }
      }
   }
   return false;
}

void renewPositionsExceptPairs(const int pair1, const int pair2, vector<LVoronoiPoint>& vps, const Boundary* bdr) {

  double cx = 0.0; double cy = 0.0; double cz = 0.0;

  const double kdvsize = (double)(vps.size());

  int index_self = 0;
  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    if ( index_self == pair1 || index_self == pair2) {
       ++index_self;
       continue;
    }
    const double px = (*v).xn2_;
    const double py = (*v).yn2_;
    const double pz = (*v).zn2_;
    (*v).x_ += px; (*v).y_ += py; (*v).z_ += pz; 
    cx += px; cy += py; cz += pz;
    ++index_self;
  }

  cx /= kdvsize-2.0; cy /= kdvsize-2.0; cz /= kdvsize-2.0;
  index_self = 0;
  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    if ( index_self == pair1 || index_self == pair2) {
       ++index_self;
       continue;
    }
    (*v).x_ -= cx; (*v).y_ -= cy; (*v).z_ -= cz;
    setPeriodicity((*v).x_,(*v).y_,(*v).z_,bdr);
    ++index_self;
  }

}

void execMonteCarloStepExceptPairs(const int pair1, const int pair2, container& base_con, vector<LVoronoiPoint>& vps, const Boundary* bdr, mt19937& mt, const double value_tarea, const double value_eratio, const double dps) {
  assert(pair1<pair2);
  //prepare constants for speed
  const double xmin = bdr->xmin_; const double xmax = bdr->xmax_;
  const double ymin = bdr->ymin_; const double ymax = bdr->ymax_;
  const double zmin = bdr->zmin_; const double zmax = bdr->zmax_;
  const int nx = bdr->nx_; const int ny = bdr->ny_; const int nz = bdr->nz_;

  const double alpha = getAlpha(value_eratio);

  normal_distribution<double> nml(0.0,1.0);

  int index = 0;
  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) {
    base_con.put(index,(*v).x_,(*v).y_,(*v).z_);
    ++index;
  }
  const double totalenergy_base = retTotalEnergy(base_con, value_tarea, value_eratio);

  int index_self = 0;
  for (std::vector<LVoronoiPoint>::iterator v = vps.begin(); v != vps.end(); ++v) { 
    if (index_self == pair1 || index_self ==pair2) {
       ++index_self;
       continue;
    }
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

  renewPositionsExceptPairs(pair1, pair2, vps, bdr);

}

#endif
