//
// Created by polya on 9/12/24.
//

#ifndef POTENTIALFUNCTIONPROTOTYPE_HPP
#define POTENTIALFUNCTIONPROTOTYPE_HPP
#include <cstring> // For memcpy
#include <fstream>
#include <immintrin.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>

class potentialFunction
{
public:
    //base class for potential function
    // virtual double operator()(const double *eta_H,const double *v0,const double *v1, const double *v2)=0;
    virtual void json2Coefs(const std::string &coefsStr)=0;
    virtual  void init()=0;
    virtual ~ potentialFunction() {};
};


std::shared_ptr<potentialFunction>  createPotentialFunction(const std::string& funcName, const std::string &row) ;


#endif //POTENTIALFUNCTIONPROTOTYPE_HPP
