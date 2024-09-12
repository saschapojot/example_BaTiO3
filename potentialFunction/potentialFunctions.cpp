//
// Created by polya on 9/12/24.
//

#include "potentialFunctionPrototype.hpp"


class V_BaTiO3: public potentialFunction
{
public:
    V_BaTiO3(const std::string &coefsStr):potentialFunction(){
        this->coefsInStr=coefsStr;
    }
    void json2Coefs(const std::string &coefsStr)override
    {

        std::stringstream iss;
        iss<<coefsStr;
        std::string temp;

        //read kappa2
        if (std::getline(iss, temp, ',')){
            this->kappa2_val=std::stod(temp);
        }

        //read alpha

        if (std::getline(iss, temp, ',')){
            this->alpha_val=std::stod(temp);
        }

        //read gamma

        if (std::getline(iss, temp, ',')){
            this->gamma_val=std::stod(temp);
        }

        //read j1

        if (std::getline(iss, temp, ',')){
            this->j1_val=std::stod(temp);
        }

        //read j2

        if (std::getline(iss, temp, ',')){
            this->j2_val=std::stod(temp);
        }

        //read j3

        if (std::getline(iss, temp, ',')){
            this->j3_val=std::stod(temp);
        }

        //read j4

        if (std::getline(iss, temp, ',')){
            this->j4_val=std::stod(temp);
        }


        //read j5

        if (std::getline(iss, temp, ',')){
            this->j5_val=std::stod(temp);
        }

        //read j6

        if (std::getline(iss, temp, ',')){
            this->j6_val=std::stod(temp);
        }

        //read j7

        if (std::getline(iss, temp, ',')){
            this->j7_val=std::stod(temp);
        }

        //read B11

        if (std::getline(iss, temp, ',')){
            this->B11_val=std::stod(temp);
        }

        //read B12

        if (std::getline(iss, temp, ',')){
            this->B12_val=std::stod(temp);
        }

        //read B44

        if (std::getline(iss, temp, ',')){
            this->B44_val=std::stod(temp);
        }

        //read B1xx

        if (std::getline(iss, temp, ',')){
            this->B1xx_val=std::stod(temp);
        }

        //read B1yy

        if (std::getline(iss, temp, ',')){
            this->B1yy_val=std::stod(temp);
        }

        //read B4yz

        if (std::getline(iss, temp, ',')){
            this->B4yz_val=std::stod(temp);
        }

        //read ZStar

        if (std::getline(iss, temp, ',')){
            this->ZStar_val=std::stod(temp);
        }

        //read epsilon_infty
        if (std::getline(iss, temp, ',')){
            this->epsilon_infty=std::stod(temp);
        }

        //read xi_Ba

        if (std::getline(iss, temp, ',')){
            this->xi_Ba=std::stod(temp);
        }

        //read xi_Ti
        if (std::getline(iss, temp, ',')){
            this->xi_Ti=std::stod(temp);
        }

        //read xi_O_parallel
        if (std::getline(iss, temp, ',')){
            this->xi_O_parallel=std::stod(temp);
        }

        //read xi_O_perpendicular

        if (std::getline(iss, temp, ',')){
            this->xi_O_perpendicular=std::stod(temp);
        }


        //read N

        if (std::getline(iss, temp, ',')){
            this->N=std::stoi(temp);
        }

        this->elemNumTot_u=N*N*N;


    }

    void init() override
    {
        this->json2Coefs(coefsInStr);

        u0=std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());
        u1=std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());
        u2=std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());

        std::cout<<"kappa2="<<kappa2_val<<", alpha="<<alpha_val<<", gamma="<<gamma_val
        <<", j1="<<j1_val<<", j2="<<j2_val<<", j3="<<j3_val
        <<", j4="<<j4_val<<", j5="<<j5_val<<", j6="<<j6_val
        <<", j7="<<j7_val<<", B11="<<B11_val<<", B12="<<B12_val
        <<", B44="<<B44_val<<", B1xx="<<B1xx_val<<", B1yy="<<B1yy_val
        <<", B4yz="<<B4yz_val<<", ZStar="<<ZStar_val<<", epsilon_infty="<<epsilon_infty
        <<", xi_Ba="<<xi_Ba<<", xi_Ti="<<xi_Ti<<", xi_O_parallel="<<xi_O_parallel
        <<", xi_O_perpendicular="<<xi_O_perpendicular
        <<std::endl;
    }

    double operator() (const double *eta_H,const double *v0,const double *v1, const double *v2)override
    {



    return 0;

    }

    double E_self(const double* u0,const double * u1, const double * u2)
    {
    double val=0;



        double val1=0;
        double val2=0;
        double val3=0;
        double tmp;
        double u02tmp;
        double u12tmp;
        double u22tmp;
        for(int j=0;j<elemNumTot_u;j++)

        {
            u02tmp=std::pow(u0[j],2.0);

            u12tmp=std::pow(u1[j],2.0);
            u22tmp=std::pow(u2[j],2.0);

            tmp=u02tmp+u12tmp+u22tmp;
            val1+=tmp;
            val2+=std::pow(tmp,2.0);
            val3+=u02tmp*u12tmp+u12tmp*u22tmp+u22tmp*u02tmp;
        }

        val1*=kappa2_val;
        val2*=alpha_val;
        val3*=gamma_val;

        val=val1+val2+val3;

        return val;





    }
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
    int elemNumTot_u;

    std::shared_ptr<double[]> u0;
    std::shared_ptr<double[]> u1;
    std::shared_ptr<double[]> u2;
};