//
//  optionBS.h
//  Option Pricing By BS-formula
//
// Jiexin Duan
// May/02/2016

#ifndef optionBS_h
#define optionBS_h

#include <iostream>
#include <string>
#include <cmath>
#include "stock.h"

using namespace std;

//cdf of standard normal distribution 
double sndcdf(double x){
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
    
    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*3.1415926,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - sndcdf(-x);
    }
}


//Pricing European Call & Put by BS-formula
class option :public stock{
private:
    int currenttime;
    double currentprice;
    double r;
    int T;
    double K;
    string type;
    stock Stock;
    int No;
    
public:
    
    option(int _No, string _type,stock  _Stock, int _currenttime, double _r, int _T, double _K){
        
        type = _type;
        Stock = _Stock;
		currenttime = Stock.GetTime();
        r = _r;
        T = _T;
        K = _K;
        No = _No;
        
    }
    
	double GetPrice(){
		return currentprice;
	}

	int GetTime(){
		return currenttime;
	}
    
    int GetNo(){
        return No;
    }
    

	void Update(){
        if (currenttime <= T) {
			
			double S = Stock.GetPrice();
            double sigma = Stock.GetSigma();
			int t = Stock.GetTime();

			currenttime = t;
            double d1 = 1/(sigma*sqrt(T - t)) * (log(S/K) + (r+sigma*sigma/2) * (T - t));
            double d2 = d1 - sigma*sqrt(T - t);
            
			//Option Pricing by BS-Formula
            if (type == "European_call") {
                currentprice = sndcdf(d1)*S - sndcdf(d2)*K*exp(-r*(T - t));
            }
            if (type == "European_put") {
                currentprice = sndcdf(-d2)*K*exp(-r*(T - t)) - sndcdf(-d1)*S;
            }
			Stock.Update();
        }

    }

};


#endif /* optionBS_h */
