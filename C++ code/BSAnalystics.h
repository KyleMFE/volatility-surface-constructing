#ifndef _BS_ANALYTICS
#define _BS_ANALYTICS

#include <cmath>

enum OptionType { Call, Put };
double cnorm(double x)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x) / sqrt(2.0);
	double t = 1.0 / (1.0 + p * x);
	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x * x);
	return 0.5*(1.0 + sign * y);
}

double bsPricer(OptionType optType, double K, double T, double S_0, double sigma, double rate)
{
	double sigmaSqrtT = sigma * std::sqrt(T);
	double d1 = (std::log(S_0 / K) + rate * T) / sigmaSqrtT + 0.5 * sigmaSqrtT;
	double d2 = d1 - sigmaSqrtT;

	double V_0;
	switch (optType)
	{
	case Call:
		V_0 = S_0 * cnorm(d1) - K * exp(-rate * T) * cnorm(d2);
		break;
	case Put:
		V_0 = K * exp(-rate * T) * cnorm(-d2) - S_0 * cnorm(-d1);
		break;
	default:
		throw "unsupported optionType";
	}
	return V_0;
}

#endif