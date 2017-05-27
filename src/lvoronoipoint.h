#ifndef LVORONOIPOINT_
#define LVORONOIPOINT_ 1

using namespace std;

class LVoronoiPoint {
public:
   LVoronoiPoint();
   LVoronoiPoint(const LVoronoiPoint & v);
   LVoronoiPoint(LVoronoiPoint && v) noexcept;
   bool operator==(const LVoronoiPoint & v) {
      return (x_ == v.x_) && (y_ == v.y_) && (z_ == v.z_);
   }
   bool operator!=(const LVoronoiPoint & v) {
      return (x_ != v.x_) || (y_ != v.y_) || (z_ != v.z_);
   }

   virtual ~LVoronoiPoint();
  
   double x_, y_, z_;
   double xn1_, yn1_, zn1_;
   double xn2_, yn2_, zn2_;
   double xo_, yo_, zo_;

};

LVoronoiPoint::LVoronoiPoint() {}
LVoronoiPoint::LVoronoiPoint(const LVoronoiPoint &v) {
   x_ = v.x_; y_ = v.y_; z_ = v.z_;
   xn1_ = v.xn1_; yn1_ = v.yn1_; zn1_ = v.zn1_;
   xn2_ = v.xn2_; yn2_ = v.yn2_; zn2_ = v.zn2_;
   xo_ = v.xo_; yo_ = v.yo_; zo_ = v.zo_;
}

LVoronoiPoint::LVoronoiPoint(LVoronoiPoint && v) noexcept {
   x_ = v.x_; y_ = v.y_; z_ = v.z_;
   xn1_ = v.xn1_; yn1_ = v.yn1_; zn1_ = v.zn1_;
   xn2_ = v.xn2_; yn2_ = v.yn2_; zn2_ = v.zn2_;
   xo_ = v.xo_; yo_ = v.yo_; zo_ = v.zo_;
}

LVoronoiPoint::~LVoronoiPoint() {}
#endif
