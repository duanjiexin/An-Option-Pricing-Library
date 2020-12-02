//
//  MonteCarlo.h
//  Monte Carlo Opition Pricing
//
//  Jiexin Duan
//  May/02/2016


#ifndef MonteCarlo_h
#define MonteCarlo_h

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <cstring>
#include <random>
#include <vector>
#include <algorithm>

using namespace std;

#define MAX_Periods 10000
#define MAX_simuNo  10000

class SingleAssetMonteCarloOptionPricing{
private:
	string type;
	double St, K, u, d, mu, sigma, qu, qd, r, R, t, T, barrier, dt;
	int nPeriod, simuNo;

public:
	SingleAssetMonteCarloOptionPricing(double _St, double _K, double _sigma, double _r, double _t, double _T, int _nPeriod, double _barrier, int _simuNo){


		St = _St;
		K = _K;
		sigma = _sigma;
		r = _r;   //riskless interest rate
		t = _t;
		T = _T;
		nPeriod = _nPeriod;
		barrier = _barrier; //if barrier =0, it isn't a barrier option, just a common option
		simuNo = _simuNo;  //Simulation times
		dt = (T - t) / nPeriod; //time interval in simulation
	}

	//Calculate European Call Price by MonteCarlo
	double GetMonteCarloEuropeanCallPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector

		for (int i=0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			StockPricePath[0] = St;
			double Payoff = 0;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T-t)/nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				//StockPricePath[j] = StockPricePath[j - 1] * exp((r - 0.5 * sigma * sigma) * dt + sigma * dWt);
			}
			Payoff = max(StockPricePath[nPeriod] - K, 0.0);
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T-t)) * PayoffSum / simuNo;
		//MonteCarloPrice = 1 / pow(1+r*dt,nPeriod) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}

	//Calculate European Put Price by MonteCarlo
	double GetMonteCarloEuropeanPutPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector

		for (int i = 0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			StockPricePath[0] = St;
			double Payoff = 0;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T - t) / nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				//StockPricePath[j] = StockPricePath[j - 1] * exp((r - 0.5 * sigma * sigma) * dt + sigma * dWt);
			}
			Payoff = max(K - StockPricePath[nPeriod], 0.0);
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T - t)) * PayoffSum / simuNo;
		//MonteCarloPrice = 1 / pow(1+r*dt,nPeriod) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}


	//Calculate Binary/Digit Call Price by MonteCarlo
	double GetMonteCarloBinaryCallPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector

		for (int i = 0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			double Payoff = 0;
			StockPricePath[0] = St;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T - t) / nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				//StockPricePath[j] = StockPricePath[j - 1] * exp((r - 0.5 * sigma * sigma) * dt + sigma * dWt);
			}
			if (StockPricePath[nPeriod] - K > 0.0)
				Payoff = 1;
			else
				Payoff = 0;
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T - t)) * PayoffSum / simuNo;
		//MonteCarloPrice = 1 / pow(1+r*dt,nPeriod) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}


	//Calculate Binary/Digit Put Price by MonteCarlo
	double GetMonteCarloBinaryPutPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector


		for (int i = 0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			double Payoff = 0;
			StockPricePath[0] = St;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T - t) / nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				//StockPricePath[j] = StockPricePath[j - 1] * exp((r - 0.5 * sigma * sigma) * dt + sigma * dWt);
			}
			if (StockPricePath[nPeriod] - K < 0.0)
				Payoff = 1;
			else
				Payoff = 0;
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T - t)) * PayoffSum / simuNo;
		//MonteCarloPrice = 1 / pow(1+r*dt,nPeriod) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}


	//Calculate Lookback Call Price by MonteCarlo
	double GetMonteCarloLookbackCallPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector

		for (int i = 0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			double Payoff = 0;
			StockPricePath[0] = St;
			double Smax = St;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T - t) / nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				if (StockPricePath[j] > Smax)
					Smax = StockPricePath[j];
			}
			Payoff = max(Smax-K, 0.0);
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T - t)) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}

	//Calculate Lookback Put Price by MonteCarlo
	double GetMonteCarloLookbackPutPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector

		for (int i = 0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			double Payoff = 0;
			StockPricePath[0] = St;
			double Smin = St;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T - t) / nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				if (StockPricePath[j] < Smin)
					Smin = StockPricePath[j];
			}
			Payoff = max(K - Smin, 0.0);
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T - t)) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}


	//Calculate Arithmetic Asian Call Price by MonteCarlo
	double GetMonteCarloArithmeticAsianCallPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector

		for (int i = 0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			double Payoff = 0;
			double SAverage = 0;
			double SSum = St;
			StockPricePath[0] = St;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T - t) / nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				SSum += StockPricePath[j];  //calculate sum of price of this price path
			}
			SAverage = SSum/(nPeriod + 1);
			Payoff = max(SAverage - K, 0.0);
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T - t)) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}


	//Calculate Arithmetic Asian Put Price by MonteCarlo
	double GetMonteCarloArithmeticAsianPutPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector

		for (int i = 0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			double Payoff = 0;
			double SAverage = 0;
			double SSum = St;
			StockPricePath[0] = St;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T - t) / nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				SSum += StockPricePath[j];  //calculate sum of price of this price path
			}
			SAverage = SSum / (nPeriod + 1);
			Payoff = max(K - SAverage, 0.0);
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T - t)) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}


	//Calculate Geometric Asian Call Price by MonteCarlo
	double GetMonteCarloGeometricAsianCallPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector

		for (int i = 0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			double Payoff = 0;
			double SAverage = log(St) / (nPeriod + 1);  //use log method to calculate geometric mean
			//double SProduct = St;
			StockPricePath[0] = St;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T - t) / nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				SAverage += log(StockPricePath[j]) / (nPeriod + 1);  //calculate sum of price of this price path
			}
			Payoff = max(exp(SAverage) - K, 0.0);
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T - t)) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}


	//Calculate Geometric Asian Put Price by MonteCarlo
	double GetMonteCarloGeometricAsianPutPrice(){
		double MonteCarloPrice = 0;
		double PayoffSum = 0; //used to calculate average of PayoffVector

		for (int i = 0; i < simuNo; i++){   //simulation several StockPrice paths
			double StockPricePath[MAX_Periods + 1];
			double Payoff = 0;
			double SAverage = log(St) / (nPeriod + 1);  //use log method to calculate geometric mean
			//double SProduct = St;
			StockPricePath[0] = St;
			for (int j = 1; j <= nPeriod; j++){    //generate each stock price path
				int seed = i * nPeriod + j;
				default_random_engine generator(seed);
				normal_distribution<double> distribution(0.0, 1.0);
				double dWt = sqrt((T - t) / nPeriod) * distribution(generator);  //define Brownian motion
				StockPricePath[j] = StockPricePath[j - 1] + StockPricePath[j - 1] * r  *dt + StockPricePath[j - 1] * sigma * dWt;
				SAverage += log(StockPricePath[j]) / (nPeriod + 1);  //calculate sum of price of this price path
			}
			Payoff = max(K - exp(SAverage), 0.0);
			PayoffSum += Payoff;
		}

		MonteCarloPrice = exp(-r*(T - t)) * PayoffSum / simuNo;
		return MonteCarloPrice;
	}
};


#endif /* MonteCarlo_h */