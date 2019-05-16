#ifndef _IMPLIEDVOL_H
#define _IMPLIEDVOL_H

#include <vector>
#include <utility>

using namespace std;

// CubicSpline interpolated smile, extrapolate flat
class Smile
{
 public:
  // constructor
  Smile( const vector< pair<double,double> >& _marks );
  
  double Vol(double strike);

 private:
  // strike to implied vol marks
  vector< pair<double, double> > marks;
  vector<double> y2; // second derivatives
  
};

class ImpliedVol 
{
 public:
  ImpliedVol( const vector< pair<double, Smile> >& );

  // linear interpolation in variance, along the strike line
  double Vol(double t, double k);

  double dVoldK(double t, double k) {return (Vol(t, k+0.01) - Vol(t, k-0.01)) / 0.02; }
  double dVoldT(double t, double k) {return (Vol(t+0.005, k) - Vol(t, k)) / 0.005; }
  double dVol2dK2(double t, double k) {return (Vol(t, k+0.01) + Vol(t, k-0.01) -  2*Vol(t, k) ) / 0.01 / 0.01; }
  
 private:
  vector< pair<double, Smile> > pillarSmiles;
};

#endif
