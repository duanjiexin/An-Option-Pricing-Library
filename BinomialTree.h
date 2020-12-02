//
//  BinomialTree.h
//  BinomialTree Opition Pricing
//
//  Jiexin Duan   
//  May/02/2016

#ifndef BinomialTree_h
#define BinomialTree_h

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <cstring>
#include <random>
#include <vector>
#include <algorithm>

using namespace std;

#define MAX_Period  200   //define MAX_Period in the tree

class SingleAssetBinomialOptionPricing{
private:
	string type;
	double St, K, u, d, sigma, qu, qd, r, R, t, T , barrier;
	int nPeriod;

public:
	SingleAssetBinomialOptionPricing(double _St, double _K, double _sigma, double _r, double _t, double _T, int _nPeriod, double _barrier){
		

		St = _St;
		K = _K;
		sigma = _sigma;
		r = _r;   //riskless interest rate
		t = _t;
		T = _T;
		nPeriod = _nPeriod;
		barrier = _barrier; //if barrier =0, it isn't a barrier option, just a common option

		//calculate parameters by CRR model
		R = exp(r*(T-t)/nPeriod);  //return rate for a period in tree

		//get u, d by CRR model
		u = exp(sigma*sqrt((T-t)/nPeriod));
		d = 1 / u;
		qu = (R-d)/(u-d);
		qd = 1-qu;
		}
	
	//Calculate European Call Price by Binomial tree
	double GetEuropeanCallPrice(){
		
		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double EuropeanCallPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			EuropeanCallPrice[i][nPeriod] = max((s[i][nPeriod] - K), 0.0);
		}

		//Calculate option price backwardly
		for (int j = nPeriod-1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				EuropeanCallPrice[i][j] = 1 / R *(qu * EuropeanCallPrice[i][j+1] + qd * EuropeanCallPrice[i + 1][j + 1]);
			}
		}
		return EuropeanCallPrice[0][0];
	}
	
	//Calculate European Put Price by Binomial tree
	double GetEuropeanPutPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double EuropeanPutPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			EuropeanPutPrice[i][nPeriod] = max((K - s[i][nPeriod]), 0.0);
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				EuropeanPutPrice[i][j] = 1 / R *(qu * EuropeanPutPrice[i][j + 1] + qd * EuropeanPutPrice[i + 1][j + 1]);
			}
		}
		return EuropeanPutPrice[0][0];
	}

	//Calculate European General Price by Binomial tree
	double GetEuropeanGeneralPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double EuropeanGeneralPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			EuropeanGeneralPrice[i][nPeriod] = s[i][nPeriod] * s[i][nPeriod]; //eg \Phi =ST^2
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				EuropeanGeneralPrice[i][j] = 1 / R *(qu * EuropeanGeneralPrice[i][j + 1] + qd * EuropeanGeneralPrice[i + 1][j + 1]);
			}
		}
		return EuropeanGeneralPrice[0][0];
	}


	//Calculate American Call Price by Binomial tree
	double GetAmericanCallPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double AmericanCallPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			AmericanCallPrice[i][nPeriod] = max((s[i][nPeriod] - K), 0.0);
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				AmericanCallPrice[i][j] = max(1 / R *(qu * AmericanCallPrice[i][j + 1] + qd * AmericanCallPrice[i + 1][j + 1]), s[i][j] - K);
			}
		}
		return AmericanCallPrice[0][0];
	}


	//Calculate American Put Price by Binomial tree
	double GetAmericanPutPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double AmericanPutPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			AmericanPutPrice[i][nPeriod] = max((K -s[i][nPeriod]), 0.0);
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				AmericanPutPrice[i][j] = max(1 / R *(qu * AmericanPutPrice[i][j + 1] + qd * AmericanPutPrice[i + 1][j + 1]), K - s[i][j] );
			}
		}
		return AmericanPutPrice[0][0];
	}


	//Calculate American General Option Price by Binomial tree
	double GetAmericanGeneralPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double AmericanGeneralPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			AmericanGeneralPrice[i][nPeriod] = s[i][nPeriod] * s[i][nPeriod];
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				AmericanGeneralPrice[i][j] = max(1 / R *(qu * AmericanGeneralPrice[i][j + 1] + qd * AmericanGeneralPrice[i + 1][j + 1]), s[i][j] * s[i][j]);
			}
		}
		return AmericanGeneralPrice[0][0];
	}


	//Calculate  down-and-out call option by Binomial tree
	double GetDownAndOutCallPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double DownAndOutCallPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			if (s[i][nPeriod] < barrier)
				DownAndOutCallPrice[i][nPeriod] = 0;
			else
				DownAndOutCallPrice[i][nPeriod] = max((s[i][nPeriod] - K), 0.0);
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				if (s[i][j] < barrier)
					DownAndOutCallPrice[i][j] = 0;
				else
					DownAndOutCallPrice[i][j] = 1 / R *(qu * DownAndOutCallPrice[i][j + 1] + qd * DownAndOutCallPrice[i + 1][j + 1]);
			}
		}
		return DownAndOutCallPrice[0][0];
	}


	//Calculate  down-and-out Put option by Binomial tree
	double GetDownAndOutPutPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double DownAndOutPutPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			if (s[i][nPeriod] < barrier)
				DownAndOutPutPrice[i][nPeriod] = 0;
			else
				DownAndOutPutPrice[i][nPeriod] = max((K - s[i][nPeriod]), 0.0);
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				if (s[i][j] < barrier)
					DownAndOutPutPrice[i][j] = 0;
				else
					DownAndOutPutPrice[i][j] = 1 / R *(qu * DownAndOutPutPrice[i][j + 1] + qd * DownAndOutPutPrice[i + 1][j + 1]);
			}
		}
		return DownAndOutPutPrice[0][0];
	}


	//Calculate  up-and-out call option by Binomial tree
	double GetUpAndOutCallPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double UpAndOutCallPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			if (s[i][nPeriod] > barrier)
				UpAndOutCallPrice[i][nPeriod] = 0;
			else
				UpAndOutCallPrice[i][nPeriod] = max((s[i][nPeriod] - K), 0.0);
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				if (s[i][j] > barrier)
					UpAndOutCallPrice[i][j] = 0;
				else
					UpAndOutCallPrice[i][j] = 1 / R *(qu * UpAndOutCallPrice[i][j + 1] + qd * UpAndOutCallPrice[i + 1][j + 1]);
			}
		}
		return UpAndOutCallPrice[0][0];
	}


	//Calculate  up-and-out Put option by Binomial tree
	double GetUpAndOutPutPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double UpAndOutPutPrice[MAX_Period + 1][MAX_Period + 1] = { 0 };
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			if (s[i][nPeriod] > barrier)
				UpAndOutPutPrice[i][nPeriod] = 0;
			else
				UpAndOutPutPrice[i][nPeriod] = max((K - s[i][nPeriod]), 0.0);
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				if (s[i][j] > barrier)
					UpAndOutPutPrice[i][j] = 0;
				else
					UpAndOutPutPrice[i][j] = 1 / R *(qu * UpAndOutPutPrice[i][j + 1] + qd * UpAndOutPutPrice[i + 1][j + 1]);
			}
		}
		return UpAndOutPutPrice[0][0];
	}

	//Calculate  down-and-in call option by Binomial tree
	double GetDownAndInCallPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}
		
		double DownAndInCallPrice_act[MAX_Period + 1][MAX_Period + 1] = { 0 }; //define price when active at time j
		double DownAndInCallPrice_not[MAX_Period + 1][MAX_Period + 1] = { 0 }; //define price when not_active at time j
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			DownAndInCallPrice_act[i][nPeriod] = max((s[i][nPeriod] - K), 0.0);
			if (s[i][nPeriod] < barrier){
				DownAndInCallPrice_not[i][nPeriod] = max((s[i][nPeriod] - K), 0.0);
			}
			else{
				DownAndInCallPrice_not[i][nPeriod] = 0;
			}
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				DownAndInCallPrice_act[i][j] = 1 / R *(qu * DownAndInCallPrice_act[i][j + 1] + qd * DownAndInCallPrice_act[i + 1][j + 1]);
				
				if (s[i][j] < barrier){
					DownAndInCallPrice_not[i][j] = DownAndInCallPrice_act[i][j];
				}
				else{
					DownAndInCallPrice_not[i][j] = 1 / R *(qu * DownAndInCallPrice_not[i][j + 1] + qd * DownAndInCallPrice_not[i + 1][j + 1]);
				}
			}
		}
		if(s[0][0] < barrier)
			return DownAndInCallPrice_act[0][0];
		else
			return DownAndInCallPrice_not[0][0];
	}


	//Calculate  down-and-in put option by Binomial tree
	double GetDownAndInPutPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double DownAndInPutPrice_act[MAX_Period + 1][MAX_Period + 1] = { 0 }; //define price when active at time j
		double DownAndInPutPrice_not[MAX_Period + 1][MAX_Period + 1] = { 0 }; //define price when not_active at time j
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			DownAndInPutPrice_act[i][nPeriod] = max((K - s[i][nPeriod]), 0.0);
			if (s[i][nPeriod] < barrier){
				DownAndInPutPrice_not[i][nPeriod] = max((K - s[i][nPeriod]), 0.0);
			}
			else{
				DownAndInPutPrice_not[i][nPeriod] = 0;
			}
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				DownAndInPutPrice_act[i][j] = 1 / R *(qu * DownAndInPutPrice_act[i][j + 1] + qd * DownAndInPutPrice_act[i + 1][j + 1]);

				if (s[i][j] < barrier){
					DownAndInPutPrice_not[i][j] = DownAndInPutPrice_act[i][j];
				}
				else{
					DownAndInPutPrice_not[i][j] = 1 / R *(qu * DownAndInPutPrice_not[i][j + 1] + qd * DownAndInPutPrice_not[i + 1][j + 1]);
				}
			}
		}
		if (s[0][0] < barrier)
			return DownAndInPutPrice_act[0][0];
		else
			return DownAndInPutPrice_not[0][0];
	}


	//Calculate  up-and-in call option by Binomial tree
	double GetUpAndInCallPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double UpAndInCallPrice_act[MAX_Period + 1][MAX_Period + 1] = { 0 }; //define price when active at time j
		double UpAndInCallPrice_not[MAX_Period + 1][MAX_Period + 1] = { 0 }; //define price when not_active at time j
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			UpAndInCallPrice_act[i][nPeriod] = max((s[i][nPeriod] - K), 0.0);
			if (s[i][nPeriod] > barrier){
				UpAndInCallPrice_not[i][nPeriod] = max((s[i][nPeriod] - K), 0.0);
			}
			else{
				UpAndInCallPrice_not[i][nPeriod] = 0;
			}
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				UpAndInCallPrice_act[i][j] = 1 / R *(qu * UpAndInCallPrice_act[i][j + 1] + qd * UpAndInCallPrice_act[i + 1][j + 1]);

				if (s[i][j] > barrier){
					UpAndInCallPrice_not[i][j] = UpAndInCallPrice_act[i][j];
				}
				else{
					UpAndInCallPrice_not[i][j] = 1 / R *(qu * UpAndInCallPrice_not[i][j + 1] + qd * UpAndInCallPrice_not[i + 1][j + 1]);
				}
			}
		}
		if (s[0][0] > barrier)
			return UpAndInCallPrice_act[0][0];
		else
			return UpAndInCallPrice_not[0][0];
	}


	//Calculate  up-and-in put option by Binomial tree
	double GetUpAndInPutPrice(){

		//generate binomial tree of stock price;
		double s[MAX_Period + 1][MAX_Period + 1] = { 0 };
		for (int j = 0; j <= nPeriod; j++){   			 // j for jth time inteval
			for (int i = 0; i <= j; i++) {  // i for ith point in jth time
				s[i][j] = St * pow(u, (j - i)) * pow(d, i);
			}
		}

		double UpAndInPutPrice_act[MAX_Period + 1][MAX_Period + 1] = { 0 }; //define price when active at time j
		double UpAndInPutPrice_not[MAX_Period + 1][MAX_Period + 1] = { 0 }; //define price when not_active at time j
		//Get the option price of T
		for (int i = 0; i <= nPeriod; i++){
			UpAndInPutPrice_act[i][nPeriod] = max((K -s[i][nPeriod]), 0.0);
			if (s[i][nPeriod] > barrier){
				UpAndInPutPrice_not[i][nPeriod] = max((K -s[i][nPeriod]), 0.0);
			}
			else{
				UpAndInPutPrice_not[i][nPeriod] = 0;
			}
		}

		//Calculate option price backwardly
		for (int j = nPeriod - 1; j >= 0; j--){   			 // i for ith time inteval
			for (int i = 0; i <= j; i++) {  // j for jth point in ith time
				UpAndInPutPrice_act[i][j] = 1 / R *(qu * UpAndInPutPrice_act[i][j + 1] + qd * UpAndInPutPrice_act[i + 1][j + 1]);

				if (s[i][j] > barrier){
					UpAndInPutPrice_not[i][j] = UpAndInPutPrice_act[i][j];
				}
				else{
					UpAndInPutPrice_not[i][j] = 1 / R *(qu * UpAndInPutPrice_not[i][j + 1] + qd * UpAndInPutPrice_not[i + 1][j + 1]);
				}
			}
		}
		if (s[0][0] > barrier)
			return UpAndInPutPrice_act[0][0];
		else
			return UpAndInPutPrice_not[0][0];
	}

};

#endif /* BinomialTree_h */