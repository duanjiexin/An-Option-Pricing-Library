//
//  portfolio.h
//  Stock Portofolio
//
//  Jiexin Duan 
//  May/02/2016

#ifndef portfolio_h
#define portfolio_h

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "stock.h"
#include "optionBS.h"

using namespace std;

class portfolio{
    
private:
    vector<stock> Stocklist;
    vector<option> Optionlist;
    vector<int> Stock_position;
    vector<int> Option_position;
    
public:
    portfolio(){
        
        vector<option> OL;
        Optionlist = OL;
        vector<int> OP;
        Option_position = OP;
        
        vector<stock> SL;
        Stocklist = SL;
        vector<int> SP;
        Stock_position = SP;
    }
    
    void InsertStock(stock _stock, int _position){
        
        int exist = 0;
        
        for (int i = 0; i<Stocklist.size(); ++i) {
            
            if (Stocklist[i].GetNo() == _stock.GetNo()) {
                Stock_position[i] += _position;
                exist = 1;
            }
            
        }
        if (!exist){
            Stocklist.push_back(_stock);
            Stock_position.push_back(_position);
            
        }
    }
    
    void InsertOption(option _option, int _position){
        
        int exist = 0;
        
        for (int i = 0; i<Optionlist.size(); ++i) {
            
            if (Optionlist[i].GetNo() == _option.GetNo()) {
                Option_position[i] += _position;
                exist = 1;
            }
            
        }
        if (!exist){
            Optionlist.push_back(_option);
            Option_position.push_back(_position);
            
        }
    }
    
    double GetPortofolioValue(){
        
        double value = 0;
        
        for (int i = 0; i<Optionlist.size(); ++i) {
            
            value += Option_position[i]*Optionlist[i].GetPrice();
            
        }
        
        for (int i = 0; i<Stocklist.size(); ++i) {
            
            value += Stock_position[i]*Stocklist[i].GetPrice();
            
        }
        
        return value;
        
    }
    
};


#endif /* portfolio_h */
