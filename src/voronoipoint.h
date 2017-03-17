#ifndef VORONOIPOINT_
#define VORONOIPOINT_ 1

using namespace std;

class VoronoiPoint {
 public:
  VoronoiPoint();
  VoronoiPoint(const VoronoiPoint & v);
  VoronoiPoint(VoronoiPoint && v) noexcept;
  bool operator==(const VoronoiPoint & v) {
    return (x_ == v.x_) && (y_ == v.y_) && (z_ == v.z_);
  }
  bool operator!=(const VoronoiPoint & v) {
    return (x_ != v.x_) || (y_ != v.y_) || (z_ != v.z_);
  }

  virtual ~VoronoiPoint();
  
  double x_, y_, z_;
  double xn1_, yn1_, zn1_;
  double xn2_, yn2_, zn2_;
  double xo_, yo_, zo_;
  double dirx_, diry_, dirz_, r_;
  double x0_, y0_, z0_;
  double nsx_, nsy_, nsz_;// to restore noise temporally

};

VoronoiPoint::VoronoiPoint() : r_(0.01) {}
VoronoiPoint::VoronoiPoint(const VoronoiPoint &v) {
  x_ = v.x_; y_ = v.y_; z_ = v.z_;
  xn1_ = v.xn1_; yn1_ = v.yn1_; zn1_ = v.zn1_;
  xn2_ = v.xn2_; yn2_ = v.yn2_; zn2_ = v.zn2_;
  xo_ = v.xo_; yo_ = v.yo_; zo_ = v.zo_;
  dirx_ = v.dirx_; diry_ = v.diry_; dirz_ = v.dirz_; r_ = v.r_;
  x0_ = v.x0_; y0_ = v.y0_; z0_ = v.z0_;
  nsx_ = v.nsx_; nsy_ = v.nsy_; nsz_ = v.nsz_;
}

VoronoiPoint::VoronoiPoint(VoronoiPoint && v) noexcept {
  x_ = v.x_; y_ = v.y_; z_ = v.z_;
  xn1_ = v.xn1_; yn1_ = v.yn1_; zn1_ = v.zn1_;
  xn2_ = v.xn2_; yn2_ = v.yn2_; zn2_ = v.zn2_;
  xo_ = v.xo_; yo_ = v.yo_; zo_ = v.zo_;
  dirx_ = v.dirx_; diry_ = v.diry_; dirz_ = v.dirz_; r_ = v.r_;
  x0_ = v.x0_; y0_ = v.y0_; z0_ = v.z0_;
  nsx_ = v.nsx_; nsy_ = v.nsy_; nsz_ = v.nsz_;
}

VoronoiPoint::~VoronoiPoint() {}
#endif
