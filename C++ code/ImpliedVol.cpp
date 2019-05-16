#include "ImpliedVol.h"
#include <iostream>
#include <cmath>

Smile::Smile( const vector< pair<double, double> >& _marks )
  : marks(_marks)
{
  int n = marks.size();
  // end y' are zero, flat extrapolation
  double yp1 = 0;
  double ypn = 0;
  y2.resize(n);
  vector<double> u(n-1);

  y2[0] = -0.5;
  u[0]=(3.0/(marks[1].first-marks[0].first)) *
    ((marks[1].second-marks[0].second) / (marks[1].first-marks[0].first) - yp1);

  for(int i = 1; i < n-1; i++) {
    double sig=(marks[i].first-marks[i-1].first)/(marks[i+1].first-marks[i-1].first);
    double p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(marks[i+1].second-marks[i].second)/(marks[i+1].first-marks[i].first)
      - (marks[i].second-marks[i-1].second)/(marks[i].first-marks[i-1].first);
    u[i]=(6.0*u[i]/(marks[i+1].first-marks[i-1].first)-sig*u[i-1])/p;
  }

  double qn=0.5;
  double un=(3.0/(marks[n-1].first-marks[n-2].first)) *
    (ypn-(marks[n-1].second-marks[n-2].second)/(marks[n-1].first-marks[n-2].first));

  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);

  
  //std::cout << "y2[" << n-1 << "] = " << y2[n-1] << std::endl;
  for (int i=n-2;i>=0;i--) {
    y2[i]=y2[i]*y2[i+1]+u[i];
    //std::cout << "y2[" << i << "] = " << y2[i] << std::endl;
  }
  
}

double Smile::Vol(double strike)
{
  int i;
  // we use trivial search, but can consider binary search for better performance
  for (i = 0; i < marks.size(); i++ ) 
    if (strike < marks[i].first )
      break; // i stores the index of the right end of the bracket

  // extrapolation
  if (i == 0)
    return marks[i].second;
  if (i == marks.size() )
    return marks[i-1].second;

  // interpolate
  double h = marks[i].first - marks[i-1].first;
  double a = (marks[i].first - strike) / h;
  double b = 1 - a;
  double c = (a*a*a - a) * h * h / 6.0;
  double d = (b*b*b - b) * h * h / 6.0;
  return a*marks[i-1].second + b*marks[i].second + c*y2[i-1] + d*y2[i];
}

ImpliedVol::ImpliedVol( const vector< pair<double, Smile> >& _pillarSmiles )
  : pillarSmiles(_pillarSmiles)
{
}

double ImpliedVol::Vol(double t, double k)
{
  // we use trivial search, but can consider binary search for better performance
  int i;
  for (i = 0; i < pillarSmiles.size(); i++ ) 
    if (t < pillarSmiles[i].first )
      break; // i stores the index of the right end of the bracket

  if (i == 0)
    return pillarSmiles[0].second.Vol(k);
  if (i == pillarSmiles.size())
    return pillarSmiles[i-1].second.Vol(k);
  double t1 = pillarSmiles[i-1].first;
  double t2 = pillarSmiles[i].first;
  double a = (t2 - t) / (t2 - t1);
  double b = 1 - a;
  double v1 = pillarSmiles[i-1].second.Vol(k);
  double v2 = pillarSmiles[i].second.Vol(k);  
  return std::sqrt( (a*v1*v1*t1 + b*v2*v2*t2)/t );
}
