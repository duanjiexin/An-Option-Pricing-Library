//
//  stock.h
//  Stock Portofolio
//
//  Jiexin Duan
//  May/02/2016


#ifndef stock_h
#define stock_h

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <cstring>
#include <random>

using namespace std;


class stock{
    
private:
    string symbol;
    string name;
    int No;  //the number for stocks i.e No1, No2, ...
    int currenttime;
    double currentprice;
    double mu;
    double sigma;
    int currentstate;
    
public:
    
    stock(){
        
        name = symbol = "N/A";
        No = 0;
        currenttime = 0;
        currentprice = 1;
        mu = 0;
        sigma = 1;
        currentstate = 1;
    }
    
    
    stock(int _No, int _currenttime, double _currentprice, double _mu, double _sigma, int _currentstate){
        
        name = symbol = "N/A";
        No = _No;
        currenttime = _currenttime;
        currentprice = _currentprice;
        mu = _mu;
        sigma = _sigma;
        currentstate = _currentstate;
    }
  
    stock(string _name, string _symbol,int _No, int _currenttime, double _currentprice, double _mu, double _sigma, int _currentstate){
        
       // name = _name;
        symbol = _symbol;
        No = _No;
        currenttime = _currenttime;
        currentprice = _currentprice;
        mu = _mu;
        sigma = _sigma;
        currentstate = _currentstate;
    }
        
    int GetTime(){
        return currenttime;
    }
    
    double GetPrice(){
        return currentprice;
    }
    
	double GetMu(){
		return mu;
	}

    double GetSigma(){
        return sigma;
    }
    
    int GetNo(){
        return No;
	}

	int GetCurrentState(){
		return currentstate;
	}
    
    void Update(){   //update  the stock price by Geometric Brownian Motion
		
        default_random_engine generator (currentstate);
        normal_distribution<double> distribution (0.0,1.0);
        
        
        currentprice *= exp(mu - 0.5*sigma*sigma +sigma*distribution(generator));

		currenttime += 1;
		currentstate += 1;
    }
};


#endif /* stock_h */
