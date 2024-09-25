//
// Created by polya on 9/12/24.
//

#ifndef POTENTIALFUNCTIONPROTOTYPE_HPP
#define POTENTIALFUNCTIONPROTOTYPE_HPP
#include <cstring> // For memcpy
#include <functional> // for std::ref
#include <fstream>
#include <immintrin.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
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

using Combination = std::tuple<int, int, int, int, int, int>;
class V_BaTiO3_parallel: public potentialFunction
{
public:
    V_BaTiO3_parallel(const std::string& coefsStr): potentialFunction()
    {
        this->coefsInStr = coefsStr;
        num_threads_pkl = std::thread::hardware_concurrency();
        // std::cout<<"num_threads_pkl="<<num_threads_pkl<<std::endl;
        results_pkl.resize(num_threads_pkl, 0.0);  // Allocate space for results_pkl
        initialize_combinations_dpl();  // Initialize the combinations_dpl
        num_combinations_pkl = combinations_dpl.size();  // Initialize num_combinations_pkl
        chunk_size_pkl = num_combinations_pkl / num_threads_pkl;  // Initialize chunk_size_pkl once
    }//end constructor

public:
  double operator()(const std::shared_ptr<double[]>& eta_H,const std::shared_ptr<double[]>& v0,const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2)override;
 void json2Coefs(const std::string &coefsStr)override;
  void init()override;
   // ~V_BaTiO3_parallel() override;
    void fill_Q();
    //convert v to u
    void v2u(const std::shared_ptr<double[]>& v, std::shared_ptr<double[]>& u);
    ///self energy
    double E_self(const std::shared_ptr<double[]>& u0, const std::shared_ptr<double[]>& u1,
                  const std::shared_ptr<double[]>& u2);

    // dipole energy
    void initialize_combinations_dpl();
    double E_dpl();
    void compute_chunk_dpl(const std::vector<Combination>& combinations_dpl, double& result, int start, int end, const std::shared_ptr<double[]>& u_left_ptr, const std::shared_ptr<double[]>& u_right_ptr, int Q_elem_ind_part7, int beta);

    //short-range energy
    double E_short();
    //short-range energy, 1NN term
    double E_short_1NN();
    //short-range energy, 2NN term
    double E_short_2NN();
    //short-range energy, 3NN term
    double E_short_3NN();

    // delta function
    double delta(const int& i, const int& j);

    //mod function
    int python_mod(int a, int M);


    //elastic energy

    double E_elas(const std::shared_ptr<double[]>& eta_H, const std::shared_ptr<double[]>& v0,
                  const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2);

    int flattened_ind_for_E_elas(const int& i, const int& j, const int& k, const int& q);

    //elastic-mode interaction
    double E_elas_mode_int(const std::shared_ptr<double[]>& eta_H, const std::shared_ptr<double[]>& v0,
                          const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2);

    int flattened_ind_for_E_elas_mode_int(const int& i, const int& j, const int& k, const int& q);
    // Helper function to clear threads
    void clear_threads();
public:
    std::string coefsInStr;
    double kappa2_val;
    double alpha_val;
    double gamma_val;
    double j1_val;
    double j3_val;
    double j6_val;
    double j2_val;
    double j4_val;
    double j7_val;
    double j5_val;
    double B11_val;
    double B12_val;
    double B44_val;
    double gamma11_val;
    double gamma12_val;
    double gamma44_val;

    int N_power_5;
    int N_power_4;
    int N_power_3;
    int N_power_2;
    int _3_power_2;


    double B100_val, B111_val, B122_val;
    double B200_val, B211_val, B222_val;
    double B300_val, B311_val, B322_val;
    double B412_val, B421_val;
    double B502_val, B520_val;
    double B601_val, B610_val;

    double B1xx_val;
    double B1yy_val;
    double B4yz_val;
    double ZStar_val;

    double epsilon_infty;

    double xi_Ba;
    double xi_Ti;
    double xi_O_parallel;
    double xi_O_perpendicular;
    int N;
    double lambda;
    int elemNumTot_u;

    std::shared_ptr<double[]> u0;
    std::shared_ptr<double[]> u1;
    std::shared_ptr<double[]> u2;
    std::shared_ptr<double[]> Q;
    std::vector<std::shared_ptr<double[]>> ptr2_u0u1u2;

    std::shared_ptr<double[]> R_hat;
    std::shared_ptr<double[]> u_left_ptr;
    std::shared_ptr<double[]> u_right_ptr;
    std::vector<Combination> combinations_dpl;
    // Threads for parallelization
    std::vector<std::thread> threads;

    // Results for each thread
    std::vector<double> results_pkl;


    int num_threads_pkl;
    int num_combinations_pkl;
    int chunk_size_pkl;


};


#endif //POTENTIALFUNCTIONPROTOTYPE_HPP
