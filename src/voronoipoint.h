#ifndef VORONOIPOINT_
#define VORONOIPOINT_ 1

using namespace std;

class VoronoiPoint {
 public:
  VoronoiPoint();
  VoronoiPoint(const VoronoiPoint & v);
  VoronoiPoint(VoronoiPoint && v) noexcept;

  virtual ~VoronoiPoint();
  
  double x_, y_, z_;
  double xn_, yn_, zn_;
  double xo_, yo_, zo_;
  double theta_, phi_, r_;
  double x0_, y0_, z0_;

};

VoronoiPoint::VoronoiPoint() : r_(0.01) {}
VoronoiPoint::VoronoiPoint(const VoronoiPoint &v) {
  x_ = v.x_; y_ = v.y_; z_ = v.z_;
  xn_ = v.xn_; yn_ = v.yn_; zn_ = v.zn_;
  xo_ = v.xo_; yo_ = v.yo_; zo_ = v.zo_;
  theta_ = v.theta_; phi_ = v.phi_; r_ = v.r_;
  x0_ = v.x0_; y0_ = v.y0_; z0_ = v.z0_;
}

VoronoiPoint::VoronoiPoint(VoronoiPoint && v) noexcept {
  x_ = v.x_; y_ = v.y_; z_ = v.z_;
  xn_ = v.xn_; yn_ = v.yn_; zn_ = v.zn_;
  xo_ = v.xo_; yo_ = v.yo_; zo_ = v.zo_;
  theta_ = v.theta_; phi_ = v.phi_; r_ = v.r_;
  x0_ = v.x0_; y0_ = v.y0_; z0_ = v.z0_;
}

VoronoiPoint::~VoronoiPoint() {}
#endif
