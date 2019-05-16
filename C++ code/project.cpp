#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <functional>
#include <cassert>
#include<map>
#include"Date.h"
#include"BSAnalystics.h"
#include"ImpliedVol.h"

using namespace std;

const double exchange_rate = 3880.0; // exchange rate of bitcon to US dollar;

// function object with template
template<typename T>
bool RootBracketingT(T f, double &a, double &b)
{
	const int NTRY = 50;
	const double FACTOR = 1.6;
	if (a >= b) throw("wrong input a and b in RootBracketing");
	double f1 = f(a);
	double f2 = f(b);
	for (int j = 0;j < NTRY;j++) {
		if (f1*f2 < 0.0) return true;
		if (std::abs(f1) < std::abs(f2))
			f1 = f(a += FACTOR * (a - b));
		else
			f2 = f(b += FACTOR * (b - a));
	}
	return false;
}

// root finding method: brent method
double rfbrent(std::function<double(double)> f, double a, double b, double tol)
{
	const int ITMAX = 100;
	const double EPS = std::numeric_limits<double>::epsilon();

	double c = b, d, e, fa = f(a), fb = f(b); //,fc,p,q,r,s,tol1,xm;
	assert(fa * fb <= 0);
	double fc = fb, p, q, r, s, tol1, xm;
	for (int iter = 0;iter < ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c = a;
			fc = fa;
			e = d = b - a;
		}
		if (std::abs(fc) < std::abs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0*EPS*std::abs(b) + 0.5*tol;
		xm = 0.5*(c - b);
		if (std::abs(xm) <= tol1 || fb == 0.0) return b;
		if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
			s = fb / fa;
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0 - s;
			}
			else {
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
				q = (q - 1.0)*(r - 1.0)*(s - 1.0);
			}
			if (p > 0.0) q = -q;
			p = std::abs(p);
			double min1 = 3.0*xm*q - std::abs(tol1*q);
			double min2 = std::abs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e = d;
				d = p / q;
			}
			else {
				d = xm;
				e = d;
			}
		}
		else {
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (std::abs(d) > tol1)
			b += d;
		else
			b += xm > 0 ? tol1 * xm : tol1 * (-xm);
		fb = f(b);
		//std::cout << "(a, b) = (" << a << ", " << b << ")" << std::endl;
	}
	throw("Maximum number of iterations exceeded in rfbrent");
}

class VolToBSPrice
{
public:
	VolToBSPrice(OptionType _optType, double _spot, double _rate, double K, double _T, double _price)
		: optType(_optType),spot(_spot), rate(_rate), strike(K), T(_T), price(_price) {}
	// vol to price function
	double operator()(double vol)
	{
		return bsPricer(optType, strike, T, spot, vol, rate) - price;
	}
private:
	OptionType optType;
	double spot;
	double rate;
	double strike;
	double T;
	double price;
};

double rfbisect(std::function<double(double)> f, double a, double b, double tol)
{
	assert(a < b && f(a) * f(b) < 0);
	double c;
	while ((b - a) / 2 > tol) {
		c = (a + b) / 2.0;
		std::cout << "(a, b) = (" << a << ", " << b << ")" << std::endl;
		if (std::abs(f(c)) < tol)
			return c;
		else {
			if (f(a)*f(c) < 0)
				b = c;
			else
				a = c;
		}
	}
	return c;
}

double getSpotPriceFromCSV(string str)
{
	ifstream inFile(str);
	string lineStr;
	vector<vector<string>> strArray;
	while (getline(inFile, lineStr))
	{
		// save to two dimension data
		stringstream ss(lineStr);
		string str;
		vector<string> lineArray;
		// seperate by ,
		while (getline(ss, str, ','))
			lineArray.push_back(str);
		strArray.push_back(lineArray);
	}
	inFile.close();
	vector<double> fwd;
	for (size_t i = 1;i < strArray.size();i++)
	{
		fwd.push_back(stod(strArray[i][5]));
	}// loading and preprocess data
	int size = fwd.size();
	double runningSum = 0.0;
	for (int i = 0; i < size;i++)
		runningSum += fwd[i];
	double spotPrice = runningSum / fwd.size(); // compute spotPrice: average of forward price
	return spotPrice;
}

vector < pair<double, double> > dataPreprocessing(string str, double T)
{
	ifstream inFile(str); 
	string lineStr;
	vector<vector<string>> strArray;
	while (getline(inFile, lineStr))
	{
		// save to two dimension data
		stringstream ss(lineStr);
		string str;
		vector<string> lineArray;
		// seperate by ,
		while (getline(ss, str, ','))
			lineArray.push_back(str);
		strArray.push_back(lineArray);
	}
	inFile.close();
	vector<OptionType> optType;
	vector<double> strike;
	vector<double> bid;
	vector<double>ask;
	vector<double> r;
	vector<double> fwd;
	map<string, OptionType> nodemap;
	nodemap["call"] = Call;
	nodemap["put"] = Put;
	for (size_t i = 1;i < strArray.size();i++)
	{
		optType.push_back(nodemap[strArray[i][0]]); // convert string to OptionType
		strike.push_back(stod(strArray[i][1]));
		bid.push_back(stod(strArray[i][2]));
		ask.push_back(stod(strArray[i][3]));
		r.push_back(stod(strArray[i][4]));
		fwd.push_back(stod(strArray[i][5]));
	}// loading and preprocess data


	int size = optType.size();
	cout << "size :" << size << endl; // print size

	
	for (int j = 0; j < size;j++)
		cout << optType[j] << " " << strike[j] << " " << bid[j] << " " << ask[j] << " " << r[j] << " " << fwd[j] << endl; // print data
	
	double runningSum = 0.0;
	for (int i = 0; i < size;i++)
		runningSum += fwd[i];
	double spotPrice = runningSum / fwd.size(); // compute spotPrice: average of forward price
	cout << "spot_price: " << spotPrice << endl;

	function<double(VolToBSPrice&, double)> fun;
	fun = &VolToBSPrice::operator();
	vector < pair<double, double> > mark;

	for (int i = 0; i < size;i++)
	{

		// To compute the volatility we use the average of the implied volatility of price ask and bid
		double a = 0.0;
		double b = 1.0;
		VolToBSPrice f1(optType[i], spotPrice, r[i], strike[i], T, ask[i] * exchange_rate);
		RootBracketingT(f1, a, b);
		double vol1 = rfbrent(f1, a, b, 1e-6);


		double c = 0.0;
		double d = 1.0;
		VolToBSPrice f2(optType[i], spotPrice, r[i], strike[i], T, bid[i] * exchange_rate);
		RootBracketingT(f2, c, d);
		double vol2 = rfbrent(f2, c, d, 1e-6);
		double vol = (vol1 + vol2) / 2.0;


		mark.push_back(pair<double, double>(strike[i], vol));
	}
	/*
	for (unsigned i = 0; i < size; i++)
		cout << mark[i].first << "\t" << mark[i].second << endl;
	*/
	// in mark vector, some strike are the same but have different volatility, so we need to set the strike with the average volatility

	vector < pair<double, double> > mark2;
	mark2.push_back(pair<double, double>(mark[0].first, mark[0].second));
	for (int i = 1; i < size; i++)
	{
		if (mark[i].first == mark[i - 1].first)
		{
			mark2[mark2.size() - 1].second = (mark[i].second + mark[i - 1].second) / 2.0;
		}
		else
		{
			mark2.push_back(pair<double, double>(mark[i].first, mark[i].second));
		}
	}
	
	for (unsigned i = 0; i < mark2.size(); i++)
		cout << mark2[i].first << " " << mark2[i].second << endl;
	
	return mark2;
}

int findTimeRange(double T, vector<double> TSet) {
	int size = TSet.size();
	if (T<TSet[0] || T>TSet[size-1] )
	{
		throw "invalid Time";
	}
	else
	{
		for (int i = 1;i < size;i++)
		{
			if (T >= TSet[i-1] && T < TSet[i])
			{
				return i - 1 ;
			}
		}
	}
}

int findStrikeRange(double S)
{
	assert(S >= 1500 && S <= 13000);
	int k = (S - 1500) / 25;
	return k;
}

double linearInterpolation(double x0, double y0, double x1, double y1, double x)
{
	assert (x0 <= x1 && x >= x0 && x <= x1);
	return y0 * (x1 - x) / (x1 - x0) + y1 * (x - x0) / (x1 - x0);
}

double volByTimeStrikeInterpolation(double T, double Strike, const vector<double>& T_set, const vector<vector<double>>& smilevol)
{
	int T_index = findTimeRange(T, T_set);
	int S_index = findStrikeRange(Strike);
	int pos = T_index * 461 + S_index;
	cout << "======================================================" << endl;
	cout << "time length:  " << T << "    index of time Interval:  " << T_index << endl;
	cout << "strike you want is :   " << Strike << endl;
	cout << "index of strike interval: " << findStrikeRange(Strike) << endl;
	cout << "======================================================" << endl;
	cout << " Interpolation range:" << endl;
	cout << smilevol[pos][0] << " " << smilevol[pos][1] << " " << smilevol[pos][2] << endl;
	cout << smilevol[pos + 1][0] << " " << smilevol[pos + 1][1] << " " << smilevol[pos + 1][2] << endl;
	double vol1 = linearInterpolation(smilevol[pos][1], smilevol[pos][2], smilevol[pos + 1][1], smilevol[pos + 1][2], Strike);
	cout << " volatility is : " << vol1 << endl;
	cout << "======================================================" << endl;

	int pos2 = (T_index + 1) * 461 + S_index;
	cout << "======================================================" << endl;
	cout << " Interpolation range:" << endl;
	cout << smilevol[pos2][0] << " " << smilevol[pos2][1] << " " << smilevol[pos2][2] << endl;
	cout << smilevol[pos2 + 1][0] << " " << smilevol[pos2 + 1][1] << " " << smilevol[pos2 + 1][2] << endl;
	double vol2 = linearInterpolation(smilevol[pos2][1], smilevol[pos2][2], smilevol[pos2 + 1][1], smilevol[pos2 + 1][2], Strike);
	cout << " volatility is : " << vol2 << endl;
	cout << "======================================================" << endl;
	// interpolate between vol1, vol2
	double vol = linearInterpolation(smilevol[pos][0], vol1, smilevol[pos2][0], vol2, T);
	return vol;
}

double spotPriceInterpolation(double T, const vector<double>& T_set, const vector<double>& Spotprice_set)
{
	int size = T_set.size();
	int pos;
	for (int i = 1; i < size;i++)
	{
		if (T >= T_set[i - 1] && T < T_set[i])
		{
			pos = i - 1;
			break;
		}
	}

	cout << "======================================================" << endl;
	cout << "Interval of spotPrice" << endl;
	cout << T_set[pos] << "  " << Spotprice_set[pos] << endl;
	cout << T_set[pos+1] << "  " << Spotprice_set[pos+1] << endl;
	cout << "======================================================" << endl;
	double spotPrice = linearInterpolation(T_set[pos], Spotprice_set[pos], T_set[pos + 1], Spotprice_set[pos + 1], T);
	return spotPrice;
}

void saveSmilesMarks(Smile &sm1, Smile &sm2, Smile &sm3, Smile &sm4, Smile &sm5, Smile &sm6, vector < pair<double, double> > &marks1, vector < pair<double, double> > &marks2, vector < pair<double, double> > &marks3, vector < pair<double, double> > &marks4, vector < pair<double, double> > &marks5, vector < pair<double, double> > &marks6)
{
	ofstream fout1("smile1.txt");
	for (int i = 0; i < 50; i++) {
		double k = 3200 + 25 * i;
		fout1 << k << "\t " << sm1.Vol(k) << endl;
	}
	ofstream fout11("marks1.txt");
	for (int i = 0; i < marks1.size(); i++) {
		fout11 << marks1[i].first << "\t" << marks1[i].second << endl;
	}
	ofstream fout2("smile2.txt");
	for (int i = 0; i < 61; i++) {
		double k = 3000 + 25 * i;
		fout2 << k << "\t " << sm2.Vol(k) << endl;
	}
	ofstream fout22("marks2.txt");
	for (int i = 0; i < marks2.size(); i++) {
		fout22 << marks2[i].first << "\t" << marks2[i].second << endl;
	}
	ofstream fout3("smile3.txt");
	for (int i = 0; i < 321; i++) {
		double k = 2000 + 25 * i;
		fout3 << k << "\t " << sm3.Vol(k) << endl;
	}
	ofstream fout33("marks3.txt");
	for (int i = 0; i < marks3.size(); i++) {
		fout33 << marks3[i].first << "\t" << marks3[i].second << endl;
	}
	ofstream fout4("smile4.txt");
	for (int i = 0; i < 141; i++) {
		double k = 2000 + 25 * i;
		fout4 << k << "\t " << sm4.Vol(k) << endl;
	}
	ofstream fout44("marks4.txt");
	for (int i = 0; i < marks4.size(); i++) {
		fout44 << marks4[i].first << "\t" << marks4[i].second << endl;
	}
	ofstream fout5("smile5.txt");
	for (int i = 0; i < 461; i++) {
		double k = 1500 + 25 * i;
		fout5 << k << "\t " << sm5.Vol(k) << endl;
	}
	ofstream fout55("marks5.txt");
	for (int i = 0; i < marks5.size(); i++) {
		fout55 << marks5[i].first << "\t" << marks5[i].second << endl;
	}
	ofstream fout6("smile6.txt");
	for (int i = 0; i < 421; i++) {
		double k = 1500 + 25 * i;
		fout6 << k << "\t " << sm6.Vol(k) << endl;
	}
	ofstream fout66("marks6.txt");
	for (int i = 0; i < marks6.size(); i++) {
		fout66 << marks6[i].first << "\t" << marks6[i].second << endl;
	}
}
int main()
{
	string str1 = "2019-03-08 08_00_00 GMT.csv";
	string str2 = "2019-03-15 08_00_00 GMT.csv";
	string str3 = "2019-03-29 08_00_00 GMT.csv";
	string str4 = "2019-04-26 08_00_00 GMT.csv";
	string str5 = "2019-06-28 08_00_00 GMT.csv";
	string str6 = "2019-09-27 08_00_00 GMT.csv";

	Date beginDate(2019, 3, 7);
	Date expiry1(2019, 3, 8);
	Date expiry2(2019, 3, 15);
	Date expiry3(2019, 3, 29);
	Date expiry4(2019, 4, 26);
	Date expiry5(2019, 6, 28);
	Date expiry6(2019, 9, 27);

	double T1 = expiry1 - beginDate;
	double T2 = expiry2 - beginDate;
	double T3 = expiry3 - beginDate;
	double T4 = expiry4 - beginDate;
	double T5 = expiry5 - beginDate;
	double T6 = expiry6 - beginDate;


	vector < pair<double, double> > marks1 = dataPreprocessing(str1, T1);
	vector < pair<double, double> > marks2 = dataPreprocessing(str2, T2);
	vector < pair<double, double> > marks3 = dataPreprocessing(str3, T3);
	vector < pair<double, double> > marks4 = dataPreprocessing(str4, T4);
	vector < pair<double, double> > marks5 = dataPreprocessing(str5, T5);
	vector < pair<double, double> > marks6 = dataPreprocessing(str6, T6);


// build a smile for each expiry
	Smile sm1(marks1);
	Smile sm2(marks2);
	Smile sm3(marks3);
	Smile sm4(marks4);
	Smile sm5(marks5);
	Smile sm6(marks6);
	saveSmilesMarks(sm1, sm2, sm3, sm4, sm5, sm6, marks1, marks2, marks3, marks4, marks5, marks6);
	
	// construct a volatility surface out of the smiles
	vector< pair<double, Smile> > pillarSmiles;
	pillarSmiles.push_back(pair<double, Smile>(T1, sm1));
	pillarSmiles.push_back(pair<double, Smile>(T2, sm2));
	pillarSmiles.push_back(pair<double, Smile>(T3, sm3));
	pillarSmiles.push_back(pair<double, Smile>(T4, sm4));
	pillarSmiles.push_back(pair<double, Smile>(T5, sm5));
	pillarSmiles.push_back(pair<double, Smile>(T6, sm6));

	ImpliedVol iv(pillarSmiles);

	vector<vector<double>> smilevol; // to save all the data of impliedVol.txt
	vector<double> maturitySet; // to save the time date of impliedVol.txt
	ofstream fout("impliedvol.txt");
	for (int i = 0; i < 13; i++) {
		double t = T1 + 0.05 * i; // range from one day to seven months
		for (int j = 0; j < 461; j++) {
			vector<double> row;
			double k = 1500 + 25 * j; // range of 1500, 13000 , interval:25
			fout << t << "\t" << k << "\t " << iv.Vol(t, k) << endl;
			row.push_back(t);
			row.push_back(k);
			row.push_back(iv.Vol(t, k));
			smilevol.push_back(row);
		}
		maturitySet.push_back(t);
	}

	//cout << "Size of implied volatility:" << smilevol.size() <<" " << smilevol[0].size() << endl;
	
	int _opt;
	double strike;
	Date date;
	cout << "======================================================" << endl;
	cout << "Please choose the option type(0,1)"  << endl;
	cin >> _opt;
	OptionType optType = OptionType(_opt);
	cout << "enter the strike of the option :" << endl;
	cin >> strike;
	cout << "enter the expiry date of the option:(2019 3 8 - 2019 9 27) " << endl;
	cin	 >> date;
	
	vector<double> timeSet;
	timeSet.push_back(T1);
	timeSet.push_back(T2);
	timeSet.push_back(T3);
	timeSet.push_back(T4);
	timeSet.push_back(T5);
	timeSet.push_back(T6);

	vector<double> spotPriceSet;
	spotPriceSet.push_back(getSpotPriceFromCSV(str1));
	spotPriceSet.push_back(getSpotPriceFromCSV(str2));
	spotPriceSet.push_back(getSpotPriceFromCSV(str3));
	spotPriceSet.push_back(getSpotPriceFromCSV(str4));
	spotPriceSet.push_back(getSpotPriceFromCSV(str5));
	spotPriceSet.push_back(getSpotPriceFromCSV(str6));

	// find mid, ask, bid price in US dollar
	double r = 0.0;
	double deltaT = date - beginDate;
	double vol = volByTimeStrikeInterpolation(deltaT, strike, maturitySet, smilevol);
	cout << "Implied volatility at time   "<< date << "at strike   "<< strike<<"   is : " << vol << endl;
	double spotPrice = spotPriceInterpolation(deltaT, timeSet, spotPriceSet);
	cout << "Spot of that option is :  " << spotPrice << endl;
	cout << "======================================================" << endl;
	double midUS = bsPricer(optType, strike, deltaT, spotPrice, vol, r);
	double askUS = bsPricer(optType, strike, deltaT, spotPrice, vol*1.005, r);
	double bidUS = bsPricer(optType, strike, deltaT, spotPrice, vol*0.995, r);
	cout << "mid price of the option is:  BTC  " << midUS/exchange_rate << endl;
	cout << "bid price of the option is:  BTC  " << bidUS/exchange_rate << endl;
	cout << "ask price of the option is:  BTC  " << askUS/exchange_rate << endl;

	// Compute Delta in BTC, Vega and Theta in USD
	// Use finite difference method
	double delta = (bsPricer(optType, strike, deltaT, spotPrice + 1, vol, r) - bsPricer(optType, strike, deltaT, spotPrice - 1, vol, r)) / 2;
	double vega =  (bsPricer(optType, strike, deltaT, spotPrice, vol+0.01, r) - bsPricer(optType, strike, deltaT, spotPrice, vol-0.01, r)) / 0.02;
	double theta = (bsPricer(optType, strike, deltaT-0.005, spotPrice, vol, r) - bsPricer(optType, strike, deltaT+0.005, spotPrice, vol, r)) / 0.01;
	cout << "======================================================" << endl;
	cout << "Option Delta in BTC:   " << delta << endl;
	cout << "Option Vega in USD:   " << vega << endl;
	cout << "Option Theta in USD:   " << theta	<< endl;

}
