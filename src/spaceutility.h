#ifndef SPACEUTILITY_H_
#define SPACEUTILITY_H_ 1

#include <cmath>
#include <algorithm>
#include <cassert>

using namespace std;

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

double getDirectionElement(double x1, double x2, double x_axe_leng) {
  double x;
   if ( fabs( x2 - x1) < fabs( x2 - x1 - x_axe_leng) ) {
      if ( fabs( x2 -x1) < fabs( x2 - x1 + x_axe_leng) ) {
         x = x2 - x1;
      } else {
         x = x2 - x1 + x_axe_leng;
      }
   } else {
      if ( fabs( x2 - x1 - x_axe_leng) < fabs( x2 - x1 + x_axe_leng) ){
        x = x2 - x1 - x_axe_leng;
      } else {
        x = x2 - x1 +x_axe_leng;
      }
   }
   assert( fabs(x) < x_axe_leng);
   return x;
}

/*
 * ret = (x, y, z, rr)
 */
void getDirectionVector(double x1, double y1, double z1, double x2, double y2, double z2, double x_axe_leng, double y_axe_leng, double z_axe_leng, vector<double> &ret) {

   double x, y, z, rr;

   x = getDirectionElement(x1,x2,x_axe_leng);
   y = getDirectionElement(y1,y2,y_axe_leng);
   z = getDirectionElement(z1,z2,z_axe_leng);
   rr = sqrt(x*x+y*y+z*z);
   ret = {x,y,z,rr};
   assert( ret.size() == 4);
}

void scaleUnitVector(double& xx, double& yy, double& zz) {
  double scale = sqrt( xx*xx + yy*yy + zz*zz );
  xx /= scale; yy /= scale; zz /= scale;
}

#endif
