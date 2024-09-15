//
// Created by polya on 9/12/24.
//

#include <boost/python/proxy.hpp>
#include <sys/stat.h>

#include "potentialFunctionPrototype.hpp"


class V_BaTiO3 : public potentialFunction
{
public:
    V_BaTiO3(const std::string& coefsStr): potentialFunction()
    {
        this->coefsInStr = coefsStr;
    }

    void json2Coefs(const std::string& coefsStr) override
    {
        std::stringstream iss;
        iss << coefsStr;
        std::string temp;
        // std::cout<<"coefsStr is "<<coefsStr<<std::endl;
        //read kappa2
        if (std::getline(iss, temp, ','))
        {
            this->kappa2_val = std::stod(temp);
        }

        //read alpha

        if (std::getline(iss, temp, ','))
        {
            this->alpha_val = std::stod(temp);
        }

        //read gamma

        if (std::getline(iss, temp, ','))
        {
            this->gamma_val = std::stod(temp);
        }

        //read j1

        if (std::getline(iss, temp, ','))
        {
            this->j1_val = std::stod(temp);
        }

        //read j2

        if (std::getline(iss, temp, ','))
        {
            this->j2_val = std::stod(temp);
        }

        //read j3

        if (std::getline(iss, temp, ','))
        {
            this->j3_val = std::stod(temp);
        }

        //read j4

        if (std::getline(iss, temp, ','))
        {
            this->j4_val = std::stod(temp);
        }


        //read j5

        if (std::getline(iss, temp, ','))
        {
            this->j5_val = std::stod(temp);
        }

        //read j6

        if (std::getline(iss, temp, ','))
        {
            this->j6_val = std::stod(temp);
        }

        //read j7

        if (std::getline(iss, temp, ','))
        {
            this->j7_val = std::stod(temp);
        }

        //read B11

        if (std::getline(iss, temp, ','))
        {
            this->B11_val = std::stod(temp);
        }

        //read B12

        if (std::getline(iss, temp, ','))
        {
            this->B12_val = std::stod(temp);
        }

        //read B44

        if (std::getline(iss, temp, ','))
        {
            this->B44_val = std::stod(temp);
        }

        //read B1xx

        if (std::getline(iss, temp, ','))
        {
            this->B1xx_val = std::stod(temp);
        }

        //read B1yy

        if (std::getline(iss, temp, ','))
        {
            this->B1yy_val = std::stod(temp);
        }

        //read B4yz

        if (std::getline(iss, temp, ','))
        {
            this->B4yz_val = std::stod(temp);
        }

        //read ZStar

        if (std::getline(iss, temp, ','))
        {
            this->ZStar_val = std::stod(temp);
        }

        //read epsilon_infty
        if (std::getline(iss, temp, ','))
        {
            this->epsilon_infty = std::stod(temp);
        }

        //read xi_Ba

        if (std::getline(iss, temp, ','))
        {
            this->xi_Ba = std::stod(temp);
        }

        //read xi_Ti
        if (std::getline(iss, temp, ','))
        {
            this->xi_Ti = std::stod(temp);
        }

        //read xi_O_parallel
        if (std::getline(iss, temp, ','))
        {
            this->xi_O_parallel = std::stod(temp);
        }

        //read xi_O_perpendicular

        if (std::getline(iss, temp, ','))
        {
            this->xi_O_perpendicular = std::stod(temp);
        }


        //read N

        if (std::getline(iss, temp, ','))
        {
            this->N = std::stoi(temp);
        }

        this->elemNumTot_u = N * N * N;

        //read lambda
        if (std::getline(iss, temp, ','))
        {
            this->lambda = std::stod(temp);
        }
    } // end json2Coefs

    void init() override
    {
        this->json2Coefs(coefsInStr);

        u0 = std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());
        u1 = std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());
        u2 = std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());

        ptr2_u0u1u2 = std::shared_ptr<std::shared_ptr<double[]>[]>(new std::shared_ptr<double[]>[3]{u0, u1, u2});

        R_hat=std::shared_ptr<double[]>(new double[3], std::default_delete<double[]>());

        int QElemNum = static_cast<int>(std::pow(N, 6) * 9);
        Q = std::shared_ptr<double[]>(new double[QElemNum], std::default_delete<double[]>());


        this->fill_Q();

        std::cout << "kappa2=" << kappa2_val << ", alpha=" << alpha_val << ", gamma=" << gamma_val
            << ", j1=" << j1_val << ", j2=" << j2_val << ", j3=" << j3_val
            << ", j4=" << j4_val << ", j5=" << j5_val << ", j6=" << j6_val
            << ", j7=" << j7_val << ", B11=" << B11_val << ", B12=" << B12_val
            << ", B44=" << B44_val << ", B1xx=" << B1xx_val << ", B1yy=" << B1yy_val
            << ", B4yz=" << B4yz_val << ", ZStar=" << ZStar_val << ", epsilon_infty=" << epsilon_infty
            << ", xi_Ba=" << xi_Ba << ", xi_Ti=" << xi_Ti << ", xi_O_parallel=" << xi_O_parallel
            << ", xi_O_perpendicular=" << xi_O_perpendicular << ", N=" << N << ", lambda=" << lambda
            << std::endl;
    }

    double operator()(const std::shared_ptr<double[]>& eta_H, const std::shared_ptr<double[]>& v0,
                      const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2) override
    {
        //v to u
        this->v2u(v0, u0);
        this->v2u(v1, u1);
        this->v2u(v2, u2);

        double energy_self = this->E_self(u0, u1, u2);
        std::cout << "energy_self=" << energy_self << std::endl;

        double energy_dipole = E_dpl();
        std::cout << "energy_dipole=" << energy_dipole << std::endl;

        double energy_short=E_short();
        std::cout<<"energy_short="<<energy_short<<std::endl;

        return 0;
    }

    void v2u(const std::shared_ptr<double[]>& v, std::shared_ptr<double[]>& u)
    {
        for (int j = 0; j < elemNumTot_u; j++)
        {
            int starting_ind = 5 * j;
            double v_Ba = v[starting_ind + 0];
            double v_Ti = v[starting_ind + 1];
            double v_O1 = v[starting_ind + 2];
            double v_O2 = v[starting_ind + 3];
            double v_O3 = v[starting_ind + 4];

            u[j] = xi_Ba * v_Ba + xi_Ti * v_Ti + xi_O_parallel * v_O1 + xi_O_perpendicular * (v_O2 + v_O3);
        }
    }

    ///self energy
    double E_self(const std::shared_ptr<double[]>& u0, const std::shared_ptr<double[]>& u1,
                  const std::shared_ptr<double[]>& u2)
    {
        double val = 0;


        double val1 = 0;
        double val2 = 0;
        double val3 = 0;
        double tmp;
        double u02tmp;
        double u12tmp;
        double u22tmp;
        for (int j = 0; j < elemNumTot_u; j++)

        {
            u02tmp = std::pow(u0[j], 2.0);

            u12tmp = std::pow(u1[j], 2.0);
            u22tmp = std::pow(u2[j], 2.0);

            tmp = u02tmp + u12tmp + u22tmp;
            val1 += tmp;
            val2 += std::pow(tmp, 2.0);
            val3 += u02tmp * u12tmp + u12tmp * u22tmp + u22tmp * u02tmp;
        }

        val1 *= kappa2_val;
        val2 *= alpha_val;
        val3 *= gamma_val;

        val = val1 + val2 + val3;

        return val;
    } // end E_self

    /// dipole energy
    double E_dpl()
    {
        double val = 0;

        for (int i1 = 0; i1 < N; i1++)
        {
            for (int j1 = 0; j1 < N; j1++)
            {
                for (int k1 = 0; k1 < N; k1++)
                {
                    for (int i2 = 0; i2 < N; i2++)
                    {
                        for (int j2 = 0; j2 < N; j2++)
                        {
                            for (int k2 = 0; k2 < N; k2++)
                            {
                                for (int alpha = 0; alpha < 3; alpha++)
                                {
                                    for (int beta = 0; beta < 3; beta++)
                                    {
                                        int Q_elem_ind = i1 * (N * N * N * N * N * 3 * 3) +
                                            j1 * (N * N * N * N * 3 * 3) +
                                            k1 * (N * N * N * 3 * 3) +
                                            i2 * (N * N * 3 * 3) +
                                            j2 * (N * 3 * 3) +
                                            k2 * (3 * 3) +
                                            alpha * 3 +
                                            beta;
                                        int u_left_ind = i1 * N * N + j1 * N + k1;
                                        int u_right_ind = i2 * N * N + j2 * N + k2;

                                        double Q_elem_val = Q[Q_elem_ind];

                                        std::shared_ptr<double[]> u_left_ptr = ptr2_u0u1u2[alpha];
                                        std::shared_ptr<double[]> u_right_ptr = ptr2_u0u1u2[beta];

                                        double u_left_elem_val = u_left_ptr[u_left_ind];
                                        double u_right_elem_val = u_right_ptr[u_right_ind];
                                        val += Q_elem_val * u_left_elem_val * u_right_elem_val;
                                    } //end beta
                                } // end alpha
                            } // end k2
                        } // end j2
                    } // end i2
                } //end k1
            } //end j1
        } //end i1

        return val;
    } // end E_dpl

    void fill_Q()
    {
        std::shared_ptr<int[]> n0n1n2_vec = std::shared_ptr<int[]>(new int[3], std::default_delete<int[]>());

        for (int i1 = 0; i1 < N; i1++)
        {
            for (int j1 = 0; j1 < N; j1++)
            {
                for (int k1 = 0; k1 < N; k1++)
                {
                    for (int i2 = 0; i2 < N; i2++)
                    {
                        for (int j2 = 0; j2 < N; j2++)
                        {
                            for (int k2 = 0; k2 < N; k2++)
                            {
                                for (int alpha = 0; alpha < 3; alpha++)
                                {
                                    for (int beta = 0; beta < 3; beta++)
                                    {
                                        double elem_val = 0;
                                        int elem_ind = i1 * (N * N * N * N * N * 3 * 3) +
                                            j1 * (N * N * N * N * 3 * 3) +
                                            k1 * (N * N * N * 3 * 3) +
                                            i2 * (N * N * 3 * 3) +
                                            j2 * (N * 3 * 3) +
                                            k2 * (3 * 3) +
                                            alpha * 3 +
                                            beta;

                                        for (int n0 = 0; n0 < N; n0++)
                                        {
                                            for (int n1 = 0; n1 < N; n1++)
                                            {
                                                for (int n2 = 0; n2 < N; n2++)
                                                {
                                                    if (n0 == 0 and n1 == 0 and n2 == 0)
                                                    {
                                                        continue;
                                                    }
                                                    n0n1n2_vec[0] = n0;
                                                    n0n1n2_vec[1] = n1;
                                                    n0n1n2_vec[2] = n2;

                                                    int n_alpha = n0n1n2_vec[alpha];
                                                    int n_beta = n0n1n2_vec[beta];

                                                    double n0_double = static_cast<double>(n0);
                                                    double n1_double = static_cast<double>(n1);
                                                    double n2_double = static_cast<double>(n2);

                                                    double i1_double = static_cast<double>(i1);
                                                    double i2_double = static_cast<double>(i2);

                                                    double j1_double = static_cast<double>(j1);
                                                    double j2_double = static_cast<double>(j2);

                                                    double k1_double = static_cast<double>(k1);
                                                    double k2_double = static_cast<double>(k2);

                                                    double n_alpha_double = static_cast<double>(n_alpha);
                                                    double n_beta_double = static_cast<double>(n_beta);

                                                    double N_double = static_cast<double>(N);

                                                    double denom_tmp = std::pow(n0_double, 2) + std::pow(n1_double, 2) +
                                                        std::pow(n2_double, 2);

                                                    elem_val += n_alpha_double * n_beta_double / denom_tmp
                                                        * std::exp(-std::pow(PI, 2) * denom_tmp / (std::pow(
                                                            N_double * lambda, 2)))
                                                        * std::cos(
                                                            2 * PI / N_double * n0_double * (i2_double - i1_double)
                                                            + 2 * PI / N_double * n1_double * (j2_double - j1_double)
                                                            + 2 * PI / N_double * n2_double * (k2_double - k1_double)
                                                        );
                                                } // end n2
                                            } // end n1
                                        } //end n0
                                        elem_val *= PI;
                                        if (alpha == beta and i1 == i2 and j1 == j2 and k1 == k2)
                                        {
                                            elem_val -= std::pow(lambda, 3) / (3 * std::sqrt(PI));
                                        }

                                        elem_val *= 2 * std::pow(ZStar_val, 2) / epsilon_infty;

                                        Q[elem_ind] = elem_val;
                                    } //end beta
                                } // end for alpha
                            } //end for k2
                        } //end for j2
                    } //end for i2
                } // end for k1
            } //end for j1
        } //end for i1
    } // end fill_Q

    double E_short()
    {
        double energy_short_1NN=E_short_1NN();
        // std::cout<<"energy_short_1NN="<<energy_short_1NN<<std::endl;

        double energy_short_2NN=E_short_2NN();
        // std::cout<<"energy_short_2NN="<<energy_short_2NN<<std::endl;

        double energy_short_3NN=E_short_3NN();
        // std::cout<<"energy_short_3NN="<<energy_short_3NN<<std::endl;

        return energy_short_1NN+energy_short_2NN+energy_short_3NN;

    }

    //short-range energy, 1NN term
    double E_short_1NN()
    {
        double val1 = 0; //term 1NN1
        for (int n0 = 0; n0 < N; n0++)
        {
            for (int n1 = 0; n1 < N; n1++)
            {
                for (int n2 = 0; n2 < N; n2++)
                {
                    for (int alpha = 0; alpha < 3; alpha++)
                    {
                        int ind_u_left = n0 * N * N + n1 * N + n2;
                        int int_u_right = n0 * N * N + n1 * N + (n2 + 1) % N;


                        std::shared_ptr<double[]> u_left_ptr = ptr2_u0u1u2[alpha];
                        std::shared_ptr<double[]> u_right_ptr = ptr2_u0u1u2[alpha];

                        double elem_left = u_left_ptr[ind_u_left];
                        double elem_right = u_right_ptr[int_u_right];

                        val1 += j2_val * elem_left * elem_right;
                    } //end alpha
                } //end n2
            } //end n1
        } //end n0


        double val2 = 0; //term 1NN2
        for (int n0 = 0; n0 < N; n0++)
        {
            for (int n1 = 0; n1 < N; n1++)
            {
                for (int n2 = 0; n2 < N; n2++)
                {
                    for (int alpha = 0; alpha < 3; alpha++)
                    {
                        int ind_u_left = n0 * N * N + n1 * N + n2;
                        int ind_u_right = n0 * N * N + ((n1 + 1) % N) * N + n2;

                        std::shared_ptr<double[]> u_left_ptr = ptr2_u0u1u2[alpha];
                        std::shared_ptr<double[]> u_right_ptr = ptr2_u0u1u2[alpha];

                        double elem_left = u_left_ptr[ind_u_left];
                        double elem_right = u_right_ptr[ind_u_right];
                        val2 += j2_val * elem_left * elem_right;
                    } //end alpha
                } //end n2
            } //end n1
        } //end n0

        double val3 = 0 ;//term 1NN3

        for (int n0 = 0; n0 < N; n0++)
        {
            for (int n1 = 0; n1 < N; n1++)
            {

                for(int n2=0;n2<N;n2++)
                {
                    for(int alpha=0;alpha<3;alpha++)
                    {
                        int ind_u_left = n0 * N * N + n1 * N + n2;
                        int ind_u_right=((n0+1)%N)*N*N+n1*N+n2;

                        std::shared_ptr<double[]> u_left_ptr = ptr2_u0u1u2[alpha];

                        std::shared_ptr<double[]>u_right_ptr=ptr2_u0u1u2[alpha];

                        double elem_left=u_left_ptr[ind_u_left];
                        double elem_right=u_right_ptr[ind_u_right];

                        val3+=j2_val*elem_left*elem_right;



                    }//end alpha
                }//end n2
            } //end n1
        } //end n0

        return val1+val2+val3;
    } //end E_short_1NN

    double E_short_2NN()
    {



    double val1=0;//term 2NN1
        for(int alpha=0;alpha<3;alpha++)
        {
            for(int beta=0;beta<3;beta++)
            {
                for(int n0=0;n0<N;n0++)
                {
                    for(int n1=0;n1<N;n1++)
                    {
                        for(int n2=0;n2<N;n2++)
                        {
                            for(int m1: {python_mod(n1-1,N),python_mod(n1+1,N)})
                            {
                                for(int m2:{python_mod(n2-1,N),python_mod(n2+1,N)})
                                {
                                    std::shared_ptr<double[]> u_left_ptr = ptr2_u0u1u2[alpha];
                                    std::shared_ptr<double[]> u_right_ptr=ptr2_u0u1u2[beta];


                                    int m0=n0;
                                    int left_ind=n0*N*N+n1*N+n2;
                                    int right_ind=m0*N*N+m1*N+m2;

                                    double elem_left=u_left_ptr[left_ind];
                                    double elem_right=u_right_ptr[right_ind];

                                    R_hat[0]=0.0;

                                    double m1_double=static_cast<double>(m1);
                                    double n1_double =static_cast<double>(n1);

                                    double m2_double=static_cast<double>(m2);
                                    double n2_double=static_cast<double>(n2);

                                    R_hat[1]=(m1_double-n1_double)/std::sqrt(2.0);
                                    R_hat[2]=(m2_double-n2_double)/std::sqrt(2.0);

                                    double J=(j4_val+std::sqrt(2.0)*(j3_val-j4_val)*std::abs(R_hat[alpha]))*delta(alpha,beta)
                                            +2.0*j5_val*R_hat[alpha]*R_hat[beta]*(1-delta(alpha,beta));

                                    val1+=J*elem_left*elem_right;

                                }//end m2
                            }//end m1

                        }//end n2
                    }//end n1
                }//end n0
            }//end beta
        }//end alpha

        val1*=0.5;


        double val2=0;//term 2NN2

        for(int alpha=0;alpha<3;alpha++)
        {
            for(int beta=0;beta<3;beta++)
            {
                for(int n0=0;n0<N;n0++)
                {
                    for(int n1=0;n1<N;n1++)
                    {
                        for(int n2=0;n2<N;n2++)
                        {
                            for(int m0 : {python_mod(n0-1,N),python_mod(n0+1,N)})
                            {
                                for(int m2 : {python_mod(n2-1,N),python_mod(n2+1,N)})
                                {
                                    std::shared_ptr<double[]> u_left_ptr = ptr2_u0u1u2[alpha];
                                    std::shared_ptr<double[]> u_right_ptr=ptr2_u0u1u2[beta];

                                    int m1=n1;
                                    int left_ind=n0*N*N+n1*N+n2;
                                    int right_ind=m0*N*N+m1*N+m2;

                                    double elem_left=u_left_ptr[left_ind];
                                    double elem_right=u_right_ptr[right_ind];

                                    double m0_double=static_cast<double>(m0);
                                    double n0_double=static_cast<double>(n0);
                                    double m2_double=static_cast<double>(m2);
                                    double n2_double=static_cast<double>(n2);

                                    R_hat[0]=(m0_double-n0_double)/std::sqrt(2.0);
                                    R_hat[1]=0.0;
                                    R_hat[2]=(m2_double-n2_double)/std::sqrt(2.0);

                                    double J=(j4_val+std::sqrt(2.0)*(j3_val-j4_val)*std::abs(R_hat[alpha]))*delta(alpha,beta)
                                           +2.0*j5_val*R_hat[alpha]*R_hat[beta]*(1-delta(alpha,beta));

                                    val2+=J*elem_left*elem_right;





                                }//end m2
                            }//end m0
                        }//end n2
                    }//end n1
                }//end n0
            }//end beta
        }//end alpha

        val2*=0.5;


        double val3=0;//term 2NN3
        for(int alpha=0;alpha<3;alpha++)
        {
            for(int beta=0;beta<3;beta++)
            {
                for(int n0=0;n0<N;n0++)
                {
                    for(int n1=0;n1<N;n1++)
                    {
                        for(int n2=0;n2<N;n2++)
                        {
                            for(int m0:{python_mod(n0-1,N),python_mod(n0+1,N)})
                            {
                                for(int m1:{python_mod(n1-1,N),python_mod(n1+1,N)})
                                {
                                    std::shared_ptr<double[]> u_left_ptr = ptr2_u0u1u2[alpha];
                                    std::shared_ptr<double[]> u_right_ptr=ptr2_u0u1u2[beta];

                                    int m2=n2;
                                    int left_ind=n0*N*N+n1*N+n2;
                                    int right_ind=m0*N*N+m1*N+m2;

                                    double elem_left=u_left_ptr[left_ind];
                                    double elem_right=u_right_ptr[right_ind];

                                    double m0_double=static_cast<double>(m0);
                                    double n0_double=static_cast<double>(n0);
                                    double m1_double=static_cast<double>(m1);
                                    double n1_double=static_cast<double>(n1);

                                    R_hat[0]=(m0_double-n0_double)/std::sqrt(2.0);
                                    R_hat[1]=(m1_double-n1_double)/std::sqrt(2.0);
                                    R_hat[2]=0.0;

                                    double J=(j4_val+std::sqrt(2.0)*(j3_val-j4_val)*std::abs(R_hat[alpha]))*delta(alpha,beta)
                                           +2.0*j5_val*R_hat[alpha]*R_hat[beta]*(1-delta(alpha,beta));


                                    val3+=J*elem_left*elem_right;
                                }//end m1
                            }//end m0
                        }//end n2
                    }//end n1
                }//end n0
            }//end beta

        }//end alpha

        val3*=0.5;

        return val1+val2+val3;

    }//end E_short_2NN


    double E_short_3NN()
    {
     double val=0;

        for(int alpha=0;alpha<3;alpha++)
        {
            for(int beta=0;beta<3;beta++)
            {
                for(int n0=0;n0<N;n0++)
                {
                    for(int n1=0;n1<N;n1++)
                    {
                        for(int n2=0;n2<N;n2++)
                        {
                            for(int m0 : {python_mod(n0-1,N),python_mod(n0+1,N)})
                            {
                                for(int m1:{python_mod(n1-1,N),python_mod(n1+1,N)})
                                {
                                    for(int m2:{python_mod(n2-1,N),python_mod(n2+1,N)})
                                    {
                                        std::shared_ptr<double[]> u_left_ptr = ptr2_u0u1u2[alpha];
                                        std::shared_ptr<double[]> u_right_ptr=ptr2_u0u1u2[beta];
                                        int left_ind=n0*N*N+n1*N+n2;
                                        int right_ind=m0*N*N+m1*N+m2;

                                        double elem_left=u_left_ptr[left_ind];
                                        double elem_right=u_right_ptr[right_ind];

                                        double m0_double=static_cast<double>(m0);
                                        double n0_double=static_cast<double>(n0);
                                        double m1_double=static_cast<double>(m1);
                                        double n1_double=static_cast<double>(n1);
                                        double m2_double=static_cast<double>(m2);
                                        double n2_double=static_cast<double>(n2);

                                        R_hat[0]=(m0_double-n0_double)/std::sqrt(3.0);
                                        R_hat[1]=(m1_double-n1_double)/std::sqrt(3.0);
                                        R_hat[2]=(m2_double-n2_double)/std::sqrt(3.0);

                                        double J=j6_val*delta(alpha,beta)+3.0*j7_val*R_hat[alpha]*R_hat[beta]*(1-delta(alpha,beta));
                                        val+=J*elem_left*elem_right;

                                    }//end m2
                                }//end m1
                            }//end m0
                        }//end n2
                    }//end n1
                }//end n0
            }//end beta
        }//end alpha
        val*=0.5;

        return val;

    }//end E_short_3NN

    ///delta function
    double delta(const int &i,const int &j)
    {
        if(i==j)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }

    int python_mod(int a, int M) {
        return (a % M + M) % M;
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
    double lambda;
    int elemNumTot_u;

    std::shared_ptr<double[]> u0;
    std::shared_ptr<double[]> u1;
    std::shared_ptr<double[]> u2;
    std::shared_ptr<double[]> Q;
    std::shared_ptr<std::shared_ptr<double[]>[]> ptr2_u0u1u2;

    std::shared_ptr<double[]> R_hat;

};


std::shared_ptr<potentialFunction> createPotentialFunction(const std::string& funcName, const std::string& coefsJsonStr)
{
    if (funcName == "V_BaTiO3")
    {
        return std::make_shared<V_BaTiO3>(coefsJsonStr);
    }

    else
    {
        throw std::invalid_argument("Unknown potential function type");
    }
}
