//
//  numericalPDE.h
//  numericalPDE Opition Pricing
//
//  Jiexin Duan
//  May/02/2016

#ifndef numericalPDE_h
#define numericalPDE_h

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <cstring>
#include <random>
#include <vector>
#include <algorithm>


using namespace std;

#define MAX_M  200    //numbers of intervals of time t
#define MAX_N  400     //numbers of intervals of space x
#define Pi    3.1415926  //define Pi


//class of numericalPDE to solve Heat Diffusion Equation in one dimension
class HeatDiffusionEquation{
private:
	double a, b, T, alpha, dt, dx, rho;
	int M, N;

public:
	HeatDiffusionEquation(double _a, double _b, double _T, int _N, int _M, double _alpha){

		a = _a;  //lower bound of space x
		b = _b;  //upper bound of space x
		T = _T;   //upper bound of time T
		M = _M;   //numbers of intervals of time t
		N = _N;   //numbers of intervals of space x
		alpha = _alpha;   //parameter alpha
		dt = T/M;     //length of time step
		dx = (b-a)/N; //length of space step
		rho = alpha*dt/(dx*dx);  //Simulation times

	}

	
	//Implementation of Explicit Method
	void ExplicitMethod(){
		double u[MAX_M + 1][MAX_N + 1];   //define u(t,x)  matrix  
		double t[MAX_M + 1];       //define grid sequnce of t, t=0,1,...,M
		double x[MAX_N + 1];       //define grid sequnce of x, t=0,1,...,N

		//Initialize u(t,x) at time 0: u(0,x)=sin(Pi*x) 
		for (int j = 0; j<= N ; j++){
			x[j] = j*dx + a;    //x_j =j/N *(b-a)+a
			u[0][j] = sin(Pi*x[j]);  //u(0,x)=sin(Pi*x) 
		}

		//Boundary Value of u(t,x) at a,b : u(t,a)=0.0, u(t,b)=0.0 
		for (int k = 0; k <= M; k++){
			t[k] = k*dt ;    //t_k =k/M 
			u[k][0] = 0.0;    //u(t,a)=0.0
			u[k][N] = 0.0;	  //u(t,b)=0.0 
		}

		//update value of u(t,x) step by step
		for (int k = 0; k <= M-1; k++){   //Time loop for t
			for (int j = 1; j <= N-1; j++){  //Space loop  for x
				u[k + 1][j] = rho * u[k][j + 1] + (1 - 2 * rho)*u[k][j] + rho * u[k][j - 1];
			}
		}

		//output result of u(T,x)
		for (int j = 0; j <= N; j++){
			cout << "x: " << x[j] << "  u(T,x): " << u[M][j] << endl;   //  output x , u(T,x) sequence
		}
		//output value of rho
		cout << "rho=" << rho << endl;
	}
	
};


//class of numericalPDE to price European Call Option by Original BSPDE
class OriginalBSPDE{
private:
	double a, b, T, r, sigma, dt, ds ,K;
	int M, N;

public:
	OriginalBSPDE(double _a, double _b, double _T, int _N, int _M, double _r, double _sigma, double _K){

		a = _a;  //lower bound of stock
		b = _b;  //upper bound of stock
		T = _T;   //upper bound of time T
		M = _M;   //numbers of intervals of time t
		N = _N;   //numbers of intervals of space x
		r = _r;   //riskless interest rate
		sigma = _sigma;  //volatility of stock
		dt = T / M;     //length of time step
		ds = (b - a) / N; //length of space step
		K =_K;  //Exercise/Strike price of option

	}


	//Implementation of Explicit Method
	void ExplicitMethod(){
		double F[MAX_M + 1][MAX_N + 1];   //define F(t,s)  matrix  
		double t[MAX_M + 1];       //define grid sequnce of t, t=0,1,...,M
		double s[MAX_N + 1];       //define grid sequnce of s, t=0,1,...,N

		//Initialize F(t,x) at time T (i.e k=0) : F(T,s)=max(ST-K,0) =max(s[0]-K,0) 
		for (int j = 0; j <= N; j++){
			s[j] = j*ds + a;    //s_j =j/N *(b-a)+a
			F[0][j] = max(s[j]-K, 0.0);  //F(T,s)=max(ST-K,0) =max(s[0]-K,0) 
		}

		//Boundary Value of F(t,s) at a,b : uF(t,a)=0.0, u(t,b)=0.0 
		for (int k = 0; k <= M; k++){
			t[k] = T - k*dt;    //t_k =T -k dt
			F[k][0] = 0.0;    //F(t,0)=0.0
			F[k][N] = s[N] - K*exp(-r*(T-t[k])); //F(t,Smax)= Smax - Kexp(-r(T-t))  if s is very large
		}

		//update value of F(t,s) step by step
		for (int k = 0; k <= M - 1; k++){   //Time loop for t
			for (int j = 1; j <= N - 1; j++){  //Space loop  for x
				F[k + 1][j] = 0.5* (sigma*sigma*j*j - r*j) *dt* F[k][j - 1] + (1 - (sigma*sigma*j*j + r)*dt)*F[k][j] + 0.5* (sigma*sigma*j*j + r*j) *dt* F[k][j + 1];
			}
		}

		//output result of u(T,x)
		for (int j = 0; j <= N; j++){
			cout << "s: " << s[j] << "  F(0,s): " << F[M][j] << endl;   //  output s , F(T,s) sequence
		}
	}

};


#endif /* numericalPDE_h */