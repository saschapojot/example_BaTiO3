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
#include <thread>
#include <vector>

const auto PI=M_PI;

class potentialFunction
{
public:
    //base class for potential function
    virtual double operator()(const std::shared_ptr<double[]>& eta_H,const std::shared_ptr<double[]>& v0,const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2)=0;
    virtual void json2Coefs(const std::string &coefsStr)=0;
    virtual  void init()=0;
    virtual ~ potentialFunction() {};
};


std::shared_ptr<potentialFunction>  createPotentialFunction(const std::string& funcName, const std::string &row) ;


#endif //POTENTIALFUNCTIONPROTOTYPE_HPP
