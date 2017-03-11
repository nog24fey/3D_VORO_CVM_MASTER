#ifndef BOUNDARY_
#define BOUNDARY_ 1

class Boundary {
 public:
  Boundary(double a) : xmin_(-a),ymin_(-a),zmin_(-a),
                       xmax_(a), ymax_(a),zmax_(a),
                       x_axe_leng_(2.0*a),y_axe_leng_(2.0*a),z_axe_leng_(2.0*a),
                       cvol_(8.0*a*a*a) {} 
  virtual ~Boundary() {};
  double xmin_,ymin_,zmin_,xmax_,ymax_,zmax_;
  double x_axe_leng_, y_axe_leng_, z_axe_leng_;
  double cvol_;                       
};

class PosSet {
 public:
  PosSet(double x, double y, double z, double dx) : ox(x), oy(y), oz(z),
                                                    vx(x+dx), vy(y+dx), vz(z+dx) {}
  virtual ~PosSet() {}
  double ox, oy, oz, vx, vy, vz;
};
#endif
