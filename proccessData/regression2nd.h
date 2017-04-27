/*
  y = a*x*x + b*x + c
*/
#ifndef REGRESSION2nd_H_
#define REGRESSION2nd_H_
class Regression2nd {
 public:
  Regression2nd(): a_(0.0), b_(0.0), c_(0.0), n_(0), Sx2_(0.0), Sx_(0.0), Sxy_(0.0), Sy_(0.0), Sx3_(0.0), Sx2y_(0.0), Sx4_(0.0) {};
  virtual ~Regression2nd() {};
  void pushData(double x, double y);
  void getResult(double &aa, double &bb, double &cc);
 private:
  double a_;
  double b_;
  double c_;

  int n_;
  double Sx2_;
  double Sx_;
  double Sxy_;
  double Sy_;
  double Sx3_;
  double Sx2y_;
  double Sx4_;
};

void Regression2nd::pushData(double x, double y) {
  ++n_;

  const double _xx = x*x;
  Sx2_ += _xx;
  Sx_ += x;
  Sxy_ += x*y;
  Sy_ += y;
  Sx3_ += _xx*x;
  Sx2y_ += _xx*y;
  Sx4_ += _xx*_xx;
}


void Regression2nd::getResult(double &aa, double &bb, double &cc) {
  const double _nn = (double)n_;
  const double _SSxx = Sx2_- Sx_*Sx_/_nn;
  const double _SSxy = Sxy_ - Sx_*Sy_/_nn;
  const double _SSxx2 = Sx3_- Sx_*Sx2_/_nn;
  const double _SSx2y = Sx2y_ - Sx2_*Sy_/_nn;
  const double _SSx2x2 = Sx4_ - Sx2_*Sx2_/_nn;
  const double _D = _SSxx*_SSx2x2 - _SSxx2*_SSxx2;

  a_ = (_SSx2y*_SSxx - _SSxy*_SSxx2) / _D;
  b_ = (_SSxy*_SSx2x2 - _SSx2y*_SSxx2) / _D;
  c_ = Sy_/_nn - b_*Sx_/_nn - a_*Sx2_/_nn;
  
  aa = a_; bb = b_; cc = c_;
}

#endif
