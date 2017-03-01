#ifndef SPACEUTILITY_H_
#define SPACEUTILITY_H_ 1

#include <cmath>

using std::min;

double squaredDistanceForMSD(double x, double y, double z, double x0, double y0, double z0, double x_axe_leng, double y_axe_leng, double z_axe_leng){
  //double r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
  //if (r2<0.25*x_axe_leng*x_axe_leng) return sqrt(r2);

  double fx = min(fabs(x-x0),fabs(x-x0-x_axe_leng));
  fx = min(fx,fabs(x-x0+x_axe_leng));
  double fy = min(fabs(y-y0),fabs(y-y0-y_axe_leng));
  fy = min(fy,fabs(y-y0+y_axe_leng));
  double fz = min(fabs(z-z0),fabs(z-z0-z_axe_leng));
  fx = min(fz,fabs(z-z0+z_axe_leng));
  return fx*fx+fy*fy+fz*fz;
}

#endif
