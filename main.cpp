//
//  main.cpp
//  Final Project: Option Pricing
//
//  Jiexin Duan
//  May/02/2016

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <cstring>
#include <random>
#include <vector>
#include <algorithm>
#include "stock.h"
#include "optionBS.h"
#include "portfolio.h"
#include "BinomialTree.h"
#include "MonteCarlo.h"
#include "numericalPDE.h"

using namespace std;


int main() {



	vector<stock> StockPool;   //define the portfolio
	int seed = 1;
    stock A1(1,0,100,0.1,0.01,seed);
    StockPool.push_back(A1);
    stock A2(2,0,100,0.2,0.02,seed);
    StockPool.push_back(A2);
    stock A3(3,0,100,0.3,0.03,seed);
    StockPool.push_back(A3);
    stock A4(4,0,100,0.4,0.04,seed);
    StockPool.push_back(A4);
    stock A5(5,0,100,0.5,0.05,seed);
    StockPool.push_back(A5);

    
    for (int i = 0; i<5; ++i) {
    
        StockPool[i].Update();
        cout << "The Price of " << i+1 << "th stock is " << StockPool[i].GetPrice() << endl;
    }

	//  Option Pricing By BS-formula
	option B1(1, "European_call", A1, 0, 0.01, 10, 60);
	option B2(1, "European_put", A1, 0, 0.01, 10, 220);
	
	cout << "The price of A1 at time " << A1.GetTime() << " is " << A1.GetPrice() << endl;
	B1.Update();
	B2.Update();
	cout << "The price of B1 at time " << B1.GetTime() << " is " << B1.GetPrice() << endl;
	cout << "The price of B2 at time " << B2.GetTime() << " is " << B2.GetPrice() << endl;
	A1.Update();
	B1.Update();
	//cout << A1.GetCurrentState() << endl;
	B2.Update();
	cout << "The price of A1 at time " << A1.GetTime() << " is " << A1.GetPrice() << endl;
	cout << "The price of B1 at time " << B1.GetTime() << " is " << B1.GetPrice() << endl;
	cout << "The price of B2 at time " << B2.GetTime() << " is " << B2.GetPrice() << endl;
	A1.Update();
	B1.Update();
	B2.Update();
	cout << "The price of A1 at time " << A1.GetTime() << " is " << A1.GetPrice() << endl;
	cout << "The price of B1 at time " << B1.GetTime() << " is " << B1.GetPrice() << endl;
	cout << "The price of B2 at time " << B2.GetTime() << " is " << B2.GetPrice() << endl;
	A1.Update();
	B1.Update();
	B2.Update();
	cout << "The price of A1 at time " << A1.GetTime() << " is " << A1.GetPrice() << endl;
	cout << "The price of B1 at time " << B1.GetTime() << " is " << B1.GetPrice() << endl;
	cout << "The price of B2 at time " << B2.GetTime() << " is " << B2.GetPrice() << endl;
		
	//Portfolio Pricing    
	portfolio P;

	P.InsertOption(B1, 1);	//Add an option
	P.InsertOption(B2, 2);
	P.InsertStock(A1, 3);	//Add a stock
	P.InsertStock(A2, 4);

	cout << "The Value of Portolio is " << P.GetPortofolioValue() << endl;
 
	//Binomial Tree Pricing for Euoropean and American option
	SingleAssetBinomialOptionPricing O1(100, 105, 0.02, 0.02, 0, 3, 200, 0);
	cout << "The Price of European Call is " << O1.GetEuropeanCallPrice() << endl;
	cout << "The Price of European Put is " << O1.GetEuropeanPutPrice() << endl;
	cout << "The Price of European General (ST^2) is " << O1.GetEuropeanGeneralPrice() << endl; //eg \Phi =ST^2
	cout << endl;

	cout << "The Price of American Call is " << O1.GetAmericanCallPrice() << endl;
	cout << "The Price of American Put is " << O1.GetAmericanPutPrice() << endl;
	cout << "The Price of American General (ST^2) is " << O1.GetAmericanGeneralPrice() << endl; //eg \Phi =ST^2
	cout << endl;

	//Binomial Tree Pricing for Barrier option
	SingleAssetBinomialOptionPricing O2(100, 105, 0.03, 0.03, 0, 1, 200, 107);
	cout << "The Price of down-and-out call option is " << O2.GetDownAndOutCallPrice() << endl;
	cout << "The Price of down-and-out put option is " << O2.GetDownAndOutPutPrice() << endl;
	cout << "The Price of up-and-out call option is " << O2.GetUpAndOutCallPrice() << endl;
	cout << "The Price of up-and-out put option is " << O2.GetUpAndOutPutPrice() << endl;
	cout << "The Price of down-and-in call option is " << O2.GetDownAndInCallPrice() << endl;
	cout << "The Price of down-and-in put option is " << O2.GetDownAndInPutPrice() << endl;
	cout << "The Price of up-and-in call option is " << O2.GetUpAndInCallPrice() << endl;
	cout << "The Price of up-and-in put option is " << O2.GetUpAndInPutPrice() << endl;
	cout << endl;

	SingleAssetBinomialOptionPricing O3(100, 102, 0.03, 0.03, 0, 1, 200, 99);
	cout << "The Price of down-and-out call option is " << O3.GetDownAndOutCallPrice() << endl;
	cout << "The Price of down-and-out put option is " << O3.GetDownAndOutPutPrice() << endl;
	cout << "The Price of up-and-out call option is " << O3.GetUpAndOutCallPrice() << endl;
	cout << "The Price of up-and-out put option is " << O3.GetUpAndOutPutPrice() << endl;
	cout << "The Price of down-and-in call option is " << O3.GetDownAndInCallPrice() << endl;
	cout << "The Price of down-and-in put option is " << O3.GetDownAndInPutPrice() << endl;
	cout << "The Price of up-and-in call option is " << O3.GetUpAndInCallPrice() << endl;
	cout << "The Price of up-and-in put option is " << O3.GetUpAndInPutPrice() << endl;
	cout << endl;

	//Monte Carlo Pricing for 0ption
	SingleAssetMonteCarloOptionPricing M1(100, 100, 0.03, 0.02, 0, 1, 500, 0, 3000);
	cout << "The Price of European Call option is " << M1.GetMonteCarloEuropeanCallPrice() << endl;
	cout << "The Price of European Put  option is " << M1.GetMonteCarloEuropeanPutPrice() << endl;
	cout << "The Price of Binary Call option is " << M1.GetMonteCarloBinaryCallPrice() << endl;
	cout << "The Price of Binary Put  option is " << M1.GetMonteCarloBinaryPutPrice() << endl;
	cout << "The Price of LookBack Call option is " << M1.GetMonteCarloLookbackCallPrice() << endl;
	cout << "The Price of LookBack Put option is " << M1.GetMonteCarloLookbackPutPrice() << endl;
	
	SingleAssetMonteCarloOptionPricing M2(100, 100, 0.03, 0.02, 0, 1, 365, 0, 3000); //daily price average
	cout << "The Price of Arithmetic Asian Call option is " << M2.GetMonteCarloArithmeticAsianCallPrice() << endl;
	cout << "The Price of Arithmetic Asian Put  option is " << M2.GetMonteCarloArithmeticAsianPutPrice() << endl;
	cout << "The Price of Geometric Asian Call option is " << M2.GetMonteCarloGeometricAsianCallPrice() << endl;
	cout << "The Price of Geometric Asian Put option is " << M2.GetMonteCarloGeometricAsianPutPrice() << endl;
	cout << endl;
	
	//numericalPDE for Heat Diffusion Equation
	HeatDiffusionEquation U1(0, 10, 0.5, 200, 200, 0.25);   //rho =alpha(T/M)/((b-a)/N)^2 = 0.25(0.5/200)/((10-0)/200)^2 =0.25
	U1.ExplicitMethod();
	HeatDiffusionEquation U2(0, 10, 0.5, 200, 200, 0.45);   //rho =alpha(T/M)/((b-a)/N)^2 = 0.25(0.5/200)/((10-0)/200)^2 =0.25
	U2.ExplicitMethod();
	HeatDiffusionEquation U3(0, 10, 0.5, 200, 200, 0.49);   //rho =alpha(T/M)/((b-a)/N)^2 = 0.25(0.5/200)/((10-0)/200)^2 =0.25
	U3.ExplicitMethod();
	HeatDiffusionEquation U4(0, 10, 0.5, 200, 200, 0.51);   //rho =alpha(T/M)/((b-a)/N)^2 = 0.25(0.5/200)/((10-0)/200)^2 =0.25
	U4.ExplicitMethod();
	HeatDiffusionEquation U5(0, 10, 0.5, 200, 200, 0.55);   //rho =alpha(T/M)/((b-a)/N)^2 = 0.25(0.5/200)/((10-0)/200)^2 =0.25
	U5.ExplicitMethod();
	HeatDiffusionEquation U6(0, 10, 0.5, 200, 200, 0.75);   //rho =alpha(T/M)/((b-a)/N)^2 = 0.25(0.5/200)/((10-0)/200)^2 =0.25
	U6.ExplicitMethod();
	

	//numericalPDE for original BSPDE
	OriginalBSPDE Call1(0.0, 200.0, 1.0, 200.0, 200.0, 0.05, 0.2, 90.0);
	Call1.ExplicitMethod();  //unstable when sigma=0.2
	OriginalBSPDE Call2(0.0,  200.0,  1.0,  200.0,  200.0,  0.05,  0.02,  90.0);
	Call2.ExplicitMethod();  //stable when sigma=0.02

	return 0;

}