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



namespace py = boost::python;
namespace np = boost::python::numpy;

class mc_computation
{
//set seed 10
private:
    int seed=10;
public:
    mc_computation(const std::string &cppInParamsFileName): e2(seed),distUnif01(0.0, 1.0)
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
                this->h_eta_H=h;
                this->h_v=10*h;
                continue;
            }// end h


            //read sweep_multiple
            if(paramCounter==10){
                iss>>sweep_multiple;
                paramCounter++;
                std::cout << "sweep_multiple=" << sweep_multiple << std::endl;

                continue;
            }//end sweep_multiple

            //read xiVec
            if(paramCounter==11){
                iss>>xiVecStr;
                paramCounter++;
                std::cout << "xiVecStr=" << xiVecStr << std::endl;

                continue;
            }//end xiVec

        }//end while

        // this->load_data(flushLastFile);
        this->potFuncPtr = createPotentialFunction(potFuncName, coefsToPotFunc);
        potFuncPtr->init();

        //allocate memory for data
        try
        {
            this->v0_data_ptr=std::shared_ptr<double[]>(new double[sweepToWrite * elemNumTot_v],
                                                        std::default_delete<double[]>());
            this->v1_data_ptr=std::shared_ptr<double[]>(new double[sweepToWrite * elemNumTot_v],
                                                        std::default_delete<double[]>());

            this->v2_data_ptr=std::shared_ptr<double[]>(new double[sweepToWrite * elemNumTot_v],
                                                        std::default_delete<double[]>());
            this->eta_H_data_ptr=std::shared_ptr<double[]>(new double[sweepToWrite * 6],
                                                        std::default_delete<double[]>());
            this->U_data_ptr=std::shared_ptr<double[]>(new double[sweepToWrite ],
                                                        std::default_delete<double[]>());
        }
        catch (const std::bad_alloc &e) {
            std::cerr << "Memory allocation error: " << e.what() << std::endl;
            std::exit(2);
        } catch (const std::exception &e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            std::exit(2);
        }


        randint_0_N_minus1=std::uniform_int_distribution<int>(0,N-1);
        // std::cout<<"randint_0_N_minus1(e2)="<<randint_0_N_minus1(e2)<<std::endl;
        elemNumTot_eta_H=6;
        randint_0_5=std::uniform_int_distribution<int>(0,5);
        randint_0_4=std::uniform_int_distribution<int>(0,4);
        // std::cout<<"randint_0_5(e2)="<<randint_0_5(e2)<<std::endl;

        std::cout<<"sweepToWrite="<<sweepToWrite<<std::endl;
        // std::cout<<"mcNum_1sweep="<<mcNum_1sweep<<std::endl;
        std::cout<<"newFlushNum="<<newFlushNum<<std::endl;
        std::cout<<"flushLastFile+1="<<flushLastFile+1<<std::endl;
        std::cout<<"TDirRoot="<<TDirRoot<<std::endl;
        std::cout<<"U_dist_dataDir="<<U_dist_dataDir<<std::endl;
        this->out_U_path=this->U_dist_dataDir+"/U/";
        if (!fs::is_directory(out_U_path) || !fs::exists(out_U_path)) {
            fs::create_directories(out_U_path);
        }

        this->out_eta_H_path=this->U_dist_dataDir+"/eta_H/";
        if (!fs::is_directory(out_eta_H_path) || !fs::exists(out_eta_H_path)) {
            fs::create_directories(out_eta_H_path);
        }

        this->out_v0_path=this->U_dist_dataDir+"/v0/";
        if (!fs::is_directory(out_v0_path) || !fs::exists(out_v0_path)) {
            fs::create_directories(out_v0_path);
        }
        this->out_v1_path=this->U_dist_dataDir+"/v1/";
        if (!fs::is_directory(out_v1_path) || !fs::exists(out_v1_path)) {
            fs::create_directories(out_v1_path);
        }
        this->out_v2_path=this->U_dist_dataDir+"/v2/";
        if (!fs::is_directory(out_v2_path) || !fs::exists(out_v2_path)) {
            fs::create_directories(out_v2_path);
        }





        //check potential value
        // double pot= (*potFuncPtr)(eta_H_init,v0_init,v1_init,v2_init);
        //  std::cout<<"pot="<<pot<<std::endl;
    }
public:

    double generate_nearby_normal(const double & x,const double &sigma);

    void proposal(const std::shared_ptr<double[]> & vecCurr,std::shared_ptr<double[]>&vecNext,const int&pos, const int &vecLength,const double& sigma);

    double acceptanceRatio( const double&UCurr, const double& UNext);

    int flattened_ind_for_v(const int& i, const int& j, const int& k, const int& q);
    void execute_mc_one_sweep(std::shared_ptr<double[]>& v0_Curr,std::shared_ptr<double[]>&v1_Curr,std::shared_ptr<double[]>&v2_Curr,
        std::shared_ptr<double[]>& eta_H_Curr,double &UCurr,
        std::shared_ptr<double[]>&v0_Next,std::shared_ptr<double[]>& v1_Next, std::shared_ptr<double[]>&v2_Next,std::shared_ptr<double[]>&eta_H_Next
        , double &U_time, double& proposal_time,double &rand_time,double &acc_reject_time);

    void execute_mc(const std::shared_ptr<double[]>& v0Vec,const std::shared_ptr<double[]>& v1Vec, const std::shared_ptr<double[]>& v2Vec,const std::shared_ptr<double[]>& eta_HVec, const int & flushNum);

    void init_and_run();


    void save_array_to_pickle(const std::shared_ptr<double[]> &ptr,const int& size,const std::string& filename);
    void load_pickle_data(const std::string& filename, std::shared_ptr<double[]>& data_ptr, std::size_t size);

    void initialize_v0_v1_v2_eta_H();

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
    double h_eta_H;//step size for eta_H
    double h_v;//step size for v
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
    int elemNumTot_eta_H;

    const double kB=1.380649e-23;
    const double Eh_unit=4.35974e-18;

    std::string coefsToPotFunc;
    std::string potFuncName;
    std::string xiVecStr;

    // int mcNum_1sweep;
    std::ranlux24_base e2;
    std::uniform_real_distribution<> distUnif01;
    std::uniform_int_distribution<int> randint_0_N_minus1;//to update v's ijk
    std::uniform_int_distribution<int> randint_0_5;//to update eta_H
    std::uniform_int_distribution<int>randint_0_4;//to update v
    int sweep_multiple;

    //initial value
    std::shared_ptr<double[]> v0_init;
    std::shared_ptr<double[]> v1_init;
    std::shared_ptr<double[]> v2_init;
    std::shared_ptr<double[]> eta_H_init;

    //data



    std::shared_ptr<double[]> v0_data_ptr;
    std::shared_ptr<double[]>v1_data_ptr;
    std::shared_ptr<double[]>v2_data_ptr;
    std::shared_ptr<double[]>eta_H_data_ptr;

    std::shared_ptr<double[]>U_data_ptr;

    std::string out_U_path;
    std::string out_eta_H_path;
    std::string out_v0_path;
    std::string out_v1_path;
    std::string out_v2_path;

};

#endif //MC_READ_LOAD_COMPUTE_HPP
