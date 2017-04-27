/*
  y = a*x + b
*/
#ifndef REGRESSION1st_H_
#define REGRESSION1st_H_
class Regression1st {
 public:
  Regression1st(): a_(0.0), b_(0.0), n_(0), Sx2_(0.0), Sx_(0.0), Sxy_(0.0), Sy_(0.0) {};
  virtual ~Regression1st() {};
  void pushData(double x, double y);
  void getResult(double &aa, double &bb);
 private:
  double a_;
  double b_;

  int n_;
  double Sx2_;
  double Sx_;
  double Sxy_;
  double Sy_;
};

void Regression1st::pushData(double x, double y) {
  ++n_;

  const double _xx = x*x;
  Sx2_ += _xx;
  Sx_ += x;
  Sxy_ += x*y;
  Sy_ += y;
}


void Regression1st::getResult(double &aa, double &bb) {
  const double _nn = (double)n_;

  const double _D = Sx2_*_nn - Sx_*Sx_;

  a_ = (Sxy_*_nn - Sx_*Sy_) / _D;
  b_ = (Sx2_*Sxy_ - Sxy_*Sx_) / _D;
  
  aa = a_; bb = b_;
}

#endif
