//
// Created by polya on 9/12/24.
//

#ifndef MC_READ_LOAD_COMPUTE_HPP
#define MC_READ_LOAD_COMPUTE_HPP
#include "../potentialFunction/potentialFunctionPrototype.hpp"
#include <boost/filesystem.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/object/pickle_support.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <cfenv> // for floating-point exceptions
#include <chrono>
#include <cstdlib>
#include <cxxabi.h>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <thread>
#include <typeinfo>
#include <vector>
namespace fs = boost::filesystem;
const auto PI=M_PI;


namespace py = boost::python;
namespace np = boost::python::numpy;

class mc_computation
{

public:
    mc_computation(const std::string &cppInParamsFileName): e2(std::random_device{}()),distUnif01(0.0, 1.0)
    {
        std::ifstream file(cppInParamsFileName);
        if (!file.is_open()) {
            std::cerr << "Failed to open the file." << std::endl;
            std::exit(20);
        }
        std::string line;
        int paramCounter = 0;

        while (std::getline(file, line))
        {
            // Check if the line is empty
            if (line.empty()) {
                continue; // Skip empty lines
            }
            std::istringstream iss(line);

            //read T
            if (paramCounter == 0) {
                iss >> T;
                if (T <= 0) {
                    std::cerr << "T must be >0" << std::endl;
                    std::exit(1);
                }//end if
                std::cout << "T=" << T << std::endl;
                this->beta = 1 / T*Eh_unit/kB;

                paramCounter++;
                continue;
            }//end reading T

            //read unit cell number
            if (paramCounter==1) {

                iss >> N;

                this->elemNumTot_u=N*N*N;
                this->elemNumTot_v=N*N*N*5;

                v0_init=std::shared_ptr<double[]>(new double[elemNumTot_v], std::default_delete<double[]>());
                v1_init=std::shared_ptr<double[]>(new double[elemNumTot_v], std::default_delete<double[]>());
                v2_init=std::shared_ptr<double[]>(new double[elemNumTot_v], std::default_delete<double[]>());

                eta_H_init=std::shared_ptr<double[]>(new double[6], std::default_delete<double[]>());

                paramCounter++;
                continue;

            }// end reading N

            //read coefficients
            if(paramCounter==2){
                iss>>coefsToPotFunc;
                paramCounter++;
                continue;

            }// end reading coefficients

            //read potential function name
            if(paramCounter==3){
                iss>>potFuncName;
                paramCounter++;
                continue;
            }//end reading potential function name


            //read sweepToWrite
            if(paramCounter==4){
                iss>>sweepToWrite;
                paramCounter++;
                continue;
            }//end reading sweepToWrite


            //read newFlushNum
            if(paramCounter==5){
                iss>>newFlushNum;
                paramCounter++;
                continue;
            }//end reading newFlushNum

            //read flushLastFile
            if(paramCounter==6){
                //if flushLastFileStr is "-1"
                //flushLastFile+1 will be 0
                iss>>flushLastFile;
                paramCounter++;
                continue;
            }//end reading flushLastFile


            //read TDirRoot
            if (paramCounter==7){
                iss>>TDirRoot;
                paramCounter++;
                continue;
            }//end reading TDirRoot

            //read U_dist_dataDir
            if(paramCounter==8){
                iss>>U_dist_dataDir;
                paramCounter++;
                continue;
            }//end reading U_dist_dataDir

            //read h
            if(paramCounter==9){
                iss>>h;
                paramCounter++;
                std::cout << "h=" << h << std::endl;

                continue;
            }// end h


            //read sweep_multiple
            if(paramCounter==10){
                iss>>sweep_multiple;
                paramCounter++;
                std::cout << "sweep_multiple=" << sweep_multiple << std::endl;

                continue;
            }//end sweep_multiple

        }//end while

        this->load_data(flushLastFile);
    }

    /// load data by flushNum
    /// @param flushNum
    void load_data(const int& flushNum);

    void load_pickle_data(const std::string& filename, std::shared_ptr<double[]> data_ptr, std::size_t size);
    template<class T>
    void print_shared_ptr(const std::shared_ptr<T> &ptr,const int& size){
        if (!ptr) {
            std::cout << "Pointer is null." << std::endl;
            return;
        }

        for(int i=0;i<size;i++){
            if(i<size-1){
                std::cout<<ptr[i]<<",";
            }
            else{
                std::cout<<ptr[i]<<std::endl;
            }
        }

    }//end print_shared_ptr
public:
    double T;// temperature
    double beta;
    double h;// step size
    int sweepToWrite;
    int newFlushNum;
    int flushLastFile;
    std::shared_ptr<potentialFunction> potFuncPtr;
    std::string TDirRoot;
    std::string U_dist_dataDir;
    std::shared_ptr<double[]> U_dist_ptr;
    int N;
    int elemNumTot_u;
    int elemNumTot_v;

    const double kB=1.380649e-23;
    const double Eh_unit=4.35974e-18;

    std::string coefsToPotFunc;
    std::string potFuncName;

    int mcNum_1sweep;
    std::ranlux24_base e2;
    std::uniform_real_distribution<> distUnif01;
    std::uniform_int_distribution<int> dist0_2N_minus1;
    int sweep_multiple;

    std::shared_ptr<double[]> v0_init;
    std::shared_ptr<double[]> v1_init;
    std::shared_ptr<double[]> v2_init;
    std::shared_ptr<double[]> eta_H_init;

};

#endif //MC_READ_LOAD_COMPUTE_HPP
