//
// Created by polya on 9/12/24.
//

#include <future>
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

        //read B100
        if (std::getline(iss, temp, ','))
        {
            this->B100_val = std::stod(temp);
        }

        //read B111
        if (std::getline(iss, temp, ','))
        {
            this->B111_val = std::stod(temp);
        }

        //read B412
        if (std::getline(iss, temp, ','))
        {
            this->B412_val = std::stod(temp);
        }
    } // end json2Coefs

    void init() override
    {
        this->json2Coefs(coefsInStr);

        this->gamma11_val = B11_val / 4.0;
        this->gamma12_val = B12_val / 8.0;
        this->gamma44_val = B44_val / 8.0;

        B211_val = B100_val;
        B322_val = B100_val;

        B122_val = B111_val;
        B200_val = B111_val;
        B222_val = B111_val;
        B300_val = B111_val;
        B311_val = B111_val;

        B421_val = B412_val;
        B502_val = B412_val;
        B520_val = B412_val;
        B601_val = B412_val;
        B610_val = B412_val;

        this->N_power_5=N*N*N*N*N;
        this->N_power_4=N*N*N*N;
        this->N_power_3=N*N*N;
        this->N_power_2=N*N;
        this->_3_power_2=3*3;

        u0 = std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());
        u1 = std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());
        u2 = std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());

        // ptr2_u0u1u2 = std::shared_ptr<std::shared_ptr<double[]>[]>(new std::shared_ptr<double[]>[3]{u0, u1, u2});
        ptr2_u0u1u2 = {u0, u1, u2};
        R_hat = std::shared_ptr<double[]>(new double[3], std::default_delete<double[]>());

        int QElemNum = static_cast<int>(std::pow(N, 6) * 9);
        Q = std::shared_ptr<double[]>(new double[QElemNum], std::default_delete<double[]>());

        this->coef_Ba=xi_Ba;
        this->coef_Ti=xi_Ti;
        this->coef_O_parallel=xi_O_parallel;
        this->coef_O_perpendicular=xi_O_perpendicular;

        this->fill_Q();

        std::cout << "kappa2=" << kappa2_val << ", alpha=" << alpha_val << ", gamma=" << gamma_val
            << ", j1=" << j1_val << ", j2=" << j2_val << ", j3=" << j3_val
            << ", j4=" << j4_val << ", j5=" << j5_val << ", j6=" << j6_val
            << ", j7=" << j7_val << ", B11=" << B11_val << ", B12=" << B12_val
            << ", B44=" << B44_val << ", B1xx=" << B1xx_val << ", B1yy=" << B1yy_val
            << ", B4yz=" << B4yz_val << ", ZStar=" << ZStar_val << ", epsilon_infty=" << epsilon_infty
            << ", xi_Ba=" << xi_Ba << ", xi_Ti=" << xi_Ti << ", xi_O_parallel=" << xi_O_parallel
            << ", xi_O_perpendicular=" << xi_O_perpendicular << ", N=" << N << ", lambda=" << lambda
            << ", gamma11_val=" << gamma11_val << ", gamma12_val=" << gamma12_val << ", gamma44_val=" << gamma44_val
            << ", B100_val=" << B100_val << ", B211_val=" << B211_val << ", B322_val" << B322_val
            << ", B111_val=" << B111_val << ", B122_val=" << B122_val << ", B200_val" << B200_val
            << ", B222_val=" << B222_val << ", B300_val=" << B300_val << ", B311_val" << B311_val
            << ", B412_val=" << B412_val << ", B421_val=" << B421_val << ", B502_val=" << B502_val
            << ", B520_val=" << B520_val << ", B601_val=" << B601_val << ", B610_val=" << B610_val
            << std::endl;
    }

    ///potential energy
    double operator()(const std::shared_ptr<double[]>& eta_H, const std::shared_ptr<double[]>& v0,
                      const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2) override
    {
        //v to u
        this->v2u(v0, u0);
        this->v2u(v1, u1);
        this->v2u(v2, u2);
        // std::cout<<"Q[1]="<<Q[1]<<std::endl;
        double energy_self = this->E_self(u0, u1, u2);
        // std::cout << "energy_self=" << energy_self << std::endl;

        double energy_dipole = E_dpl();
        // std::cout << "energy_dipole=" << energy_dipole << std::endl;

        double energy_short = E_short();
        // std::cout << "energy_short=" << energy_short << std::endl;

        double energy_elas = E_elas(eta_H, v0, v1, v2);
        // std::cout << "energy_elas=" << energy_elas << std::endl;

        double energy_elas_mode_int=E_elas_mode_int(eta_H,v0,v1,v2);
        // std::cout<<"energy_elas_mode_int="<<energy_elas_mode_int<<std::endl;

        double pot_energy=energy_self+energy_dipole+energy_short+energy_elas+energy_elas_mode_int;
        return pot_energy;
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

            // u[j] = xi_Ba * v_Ba + xi_Ti * v_Ti + xi_O_parallel * v_O1 + xi_O_perpendicular * (v_O2 + v_O3);
            u[j]=coef_Ba * v_Ba + coef_Ti * v_Ti + coef_O_parallel * v_O1 + coef_O_perpendicular * (v_O2 + v_O3);
        }
    }//end v2u

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
        int Q_elem_ind_part1;
        int Q_elem_ind_part2;
        int Q_elem_ind_part3;
        int Q_elem_ind_part4;
        int Q_elem_ind_part5;
        int Q_elem_ind_part6;
        int Q_elem_ind_part7;
        int Q_elem_ind;
        int u_left_ind;
        int u_right_ind;

        int u_left_ind_part1;
        int u_left_ind_part2;
        int u_right_ind_part1;
        int u_right_ind_part2;
        // std::shared_ptr<double[]> u_left_ptr;
        // std::shared_ptr<double[]> u_right_ptr;

        double Q_elem_val;
        double u_left_elem_val ;
        double u_right_elem_val;
        for (int alpha = 0; alpha < 3; alpha++)
    {
        Q_elem_ind_part7 = alpha * 3;
            u_left_ptr = ptr2_u0u1u2[alpha];

        for (int beta = 0; beta < 3; beta++)
        {u_right_ptr = ptr2_u0u1u2[beta];
            for (int i1 = 0; i1 < N; i1++)
            {
                Q_elem_ind_part1 = i1 * (N_power_5 * _3_power_2);
                u_left_ind_part1 = i1 * N_power_2;

                for (int j1 = 0; j1 < N; j1++)
                {
                    Q_elem_ind_part2 = j1 * (N_power_4 * _3_power_2);
                    u_left_ind_part2 = j1 * N;

                    for (int k1 = 0; k1 < N; k1++)
                    {
                        Q_elem_ind_part3 = k1 * (N_power_3 * _3_power_2);

                        for (int i2 = 0; i2 < N; i2++)
                        {
                            Q_elem_ind_part4 = i2 * (N_power_2 * _3_power_2);
                            u_right_ind_part1 = i2 * N_power_2;

                            for (int j2 = 0; j2 < N; j2++)
                            {
                                Q_elem_ind_part5 = j2 * (N * _3_power_2);
                                u_right_ind_part2 = j2 * N;

                                for (int k2 = 0; k2 < N; k2++)
                                {
                                    Q_elem_ind_part6 = k2 * (_3_power_2);

                                    // Final index computations
                                    Q_elem_ind = Q_elem_ind_part1 + Q_elem_ind_part2 + Q_elem_ind_part3 +
                                                 Q_elem_ind_part4 + Q_elem_ind_part5 + Q_elem_ind_part6 +
                                                 Q_elem_ind_part7 + beta;

                                    u_left_ind = u_left_ind_part1 + u_left_ind_part2 + k1;
                                    u_right_ind = u_right_ind_part1 + u_right_ind_part2 + k2;

                                    // Access values from Q and u_left, u_right
                                    Q_elem_val = Q[Q_elem_ind];

                                    // u_left_ptr = ptr2_u0u1u2[alpha];
                                    // u_right_ptr = ptr2_u0u1u2[beta];

                                    u_left_elem_val = u_left_ptr[u_left_ind];
                                    u_right_elem_val = u_right_ptr[u_right_ind];

                                    // Accumulate the result
                                    val += Q_elem_val * u_left_elem_val * u_right_elem_val;
                                } // end k2
                            } // end j2
                        } // end i2
                    } // end k1
                } // end j1
            } // end i1
        } // end beta
    } // end alpha

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
        double energy_short_1NN = E_short_1NN();
        // std::cout<<"energy_short_1NN="<<energy_short_1NN<<std::endl;

        double energy_short_2NN = E_short_2NN();
        // std::cout<<"energy_short_2NN="<<energy_short_2NN<<std::endl;

        double energy_short_3NN = E_short_3NN();
        // std::cout<<"energy_short_3NN="<<energy_short_3NN<<std::endl;

        return energy_short_1NN + energy_short_2NN + energy_short_3NN;
    }

    //short-range energy, 1NN term
    double E_short_1NN()
    {
        int ind_u_left ;
        int ind_u_right;

        int ind_u_left_part1;
        int ind_u_left_part2;

        int ind_u_right_part1;
        int ind_u_right_part2;

        // std::shared_ptr<double[]> u_left_ptr;
        // std::shared_ptr<double[]> u_right_ptr;
        double elem_left;
        double elem_right;

        double val1 = 0; //term 1NN1

        for (int alpha = 0; alpha < 3; alpha++) // Moving alpha loop to outermost
        {
            // Access the pointers for the current alpha
            u_left_ptr = ptr2_u0u1u2[alpha];
            u_right_ptr = ptr2_u0u1u2[alpha];

            for (int n0 = 0; n0 < N; n0++)
            {
                ind_u_left_part1 = n0 * N_power_2; // Calculated once for n0

                for (int n1 = 0; n1 < N; n1++)
                {
                    ind_u_left_part2 = n1 * N; // Calculated once for n1

                    for (int n2 = 0; n2 < N; n2++)
                    {
                        ind_u_left = ind_u_left_part1 + ind_u_left_part2 + n2;
                        ind_u_right = ind_u_left_part1 + ind_u_left_part2 + (n2 + 1) % N;

                        // Get the left and right elements
                        elem_left = u_left_ptr[ind_u_left];
                        elem_right = u_right_ptr[ind_u_right];

                        // Accumulate the value
                        val1 += j2_val * elem_left * elem_right;
                    } // end n2
                } // end n1
            } // end n0
        }//end alpha



        double val2 = 0; //term 1NN2
        for (int alpha = 0; alpha < 3; alpha++) // Move alpha loop to the outermost position
        {
            // Access the pointers for the current alpha
            u_left_ptr = ptr2_u0u1u2[alpha];
            u_right_ptr = ptr2_u0u1u2[alpha];

            for (int n0 = 0; n0 < N; n0++)
            {
                ind_u_left_part1 = n0 * N_power_2; // Calculated once per n0

                for (int n1 = 0; n1 < N; n1++)
                {
                    ind_u_left_part2 = n1 * N;               // Calculated once per n1
                    ind_u_right_part2 = ((n1 + 1) % N) * N;  // Wrapped index for u_right

                    for (int n2 = 0; n2 < N; n2++)
                    {
                        ind_u_left = ind_u_left_part1 + ind_u_left_part2 + n2;
                        ind_u_right = ind_u_left_part1 + ind_u_right_part2 + n2;

                        // Get the left and right elements for the current alpha
                        elem_left = u_left_ptr[ind_u_left];
                        elem_right = u_right_ptr[ind_u_right];

                        // Accumulate the value
                        val2 += j2_val * elem_left * elem_right;
                    } // end n2
                } // end n1
            } // end n0
        }//end alpha



        double val3 = 0; //term 1NN3
        for (int alpha = 0; alpha < 3; alpha++) // Move alpha loop to the outermost position
        {
            // Access the pointers for the current alpha
            u_left_ptr = ptr2_u0u1u2[alpha];
            u_right_ptr = ptr2_u0u1u2[alpha];

            for (int n0 = 0; n0 < N; n0++)
            {
                ind_u_left_part1 = n0 * N_power_2;           // Calculated once per n0
                ind_u_right_part1 = ((n0 + 1) % N) * N_power_2; // Wrap n0 for u_right

                for (int n1 = 0; n1 < N; n1++)
                {
                    ind_u_left_part2 = n1 * N; // Calculated once per n1

                    for (int n2 = 0; n2 < N; n2++)
                    {
                        ind_u_left = ind_u_left_part1 + ind_u_left_part2 + n2;
                        ind_u_right = ind_u_right_part1 + ind_u_left_part2 + n2;

                        // Get the left and right elements for the current alpha
                        elem_left = u_left_ptr[ind_u_left];
                        elem_right = u_right_ptr[ind_u_right];

                        // Accumulate the value
                        val3 += j2_val * elem_left * elem_right;
                    } // end n2
                } // end n1
            } // end n0
        }//end alpha

        return val1 + val2 + val3;
    } //end E_short_1NN

    double E_short_2NN()
    {
        double val1 = 0; //term 2NN1
        // std::shared_ptr<double[]> u_left_ptr;
        // std::shared_ptr<double[]> u_right_ptr;
        int m0 ;
        int left_ind;
        int right_ind;
        int left_ind_part1;
        int left_ind_part2;
        int right_ind_part1;
        int right_ind_part2;
        double elem_left ;
        double elem_right;
        double m1_double;
        double n1_double;
        double m2_double;
        double n2_double;
        double J;
        for (int alpha = 0; alpha < 3; alpha++)
        {
            u_left_ptr = ptr2_u0u1u2[alpha];
            for (int beta = 0; beta < 3; beta++)
            {
                u_right_ptr = ptr2_u0u1u2[beta];
                for (int n0 = 0; n0 < N; n0++)
                {
                    m0 = n0;
                    left_ind_part1=n0*N_power_2;
                    right_ind_part1=m0*N_power_2;
                    for (int n1 = 0; n1 < N; n1++)
                    {
                        left_ind_part2=n1*N;
                        for (int n2 = 0; n2 < N; n2++)
                        {
                            left_ind=left_ind_part1+left_ind_part2+n2;
                            for (int m1 : {python_mod(n1 - 1, N), python_mod(n1 + 1, N)})
                            {
                                right_ind_part2=m1*N;
                                for (int m2 : {python_mod(n2 - 1, N), python_mod(n2 + 1, N)})
                                {
                                     // u_left_ptr = ptr2_u0u1u2[alpha];
                                     // u_right_ptr = ptr2_u0u1u2[beta];


                                   // m0 = n0;
                                     // left_ind = n0 * N * N + n1 * N + n2;
                                    // left_ind=left_ind_part1+left_ind_part2+n2;
                                     // right_ind = m0 * N * N + m1 * N + m2;
                                    right_ind=right_ind_part1+right_ind_part2+m2;

                                     elem_left = u_left_ptr[left_ind];
                                     elem_right = u_right_ptr[right_ind];

                                    R_hat[0] = 0.0;

                                    m1_double = static_cast<double>(m1);
                                    n1_double = static_cast<double>(n1);

                                    m2_double = static_cast<double>(m2);
                                     n2_double = static_cast<double>(n2);

                                    R_hat[1] = (m1_double - n1_double) / std::sqrt(2.0);
                                    R_hat[2] = (m2_double - n2_double) / std::sqrt(2.0);

                                     J = (j4_val + std::sqrt(2.0) * (j3_val - j4_val) * std::abs(R_hat[alpha])) *
                                        delta(alpha, beta)
                                        + 2.0 * j5_val * R_hat[alpha] * R_hat[beta] * (1 - delta(alpha, beta));

                                    val1 += J * elem_left * elem_right;
                                } //end m2
                            } //end m1
                        } //end n2
                    } //end n1
                } //end n0
            } //end beta
        } //end alpha

        val1 *= 0.5;


        double val2 = 0; //term 2NN2
        int m1;
        double m0_double;
        double n0_double;

        for (int alpha = 0; alpha < 3; alpha++)
        {
            u_left_ptr = ptr2_u0u1u2[alpha];
            for (int beta = 0; beta < 3; beta++)
            {
                u_right_ptr = ptr2_u0u1u2[beta];
                for (int n0 = 0; n0 < N; n0++)
                {
                    left_ind_part1=n0*N_power_2;
                    for (int n1 = 0; n1 < N; n1++)
                    {
                        m1 = n1;
                        left_ind_part2=n1*N;
                        right_ind_part2=m1*N;
                        for (int n2 = 0; n2 < N; n2++)
                        {
                            left_ind = left_ind_part1+ left_ind_part2 + n2;
                            for (int m0 : {python_mod(n0 - 1, N), python_mod(n0 + 1, N)})
                            {
                                right_ind_part1=m0*N_power_2;
                                for (int m2 : {python_mod(n2 - 1, N), python_mod(n2 + 1, N)})
                                {
                                   // u_left_ptr = ptr2_u0u1u2[alpha];
                                     // u_right_ptr = ptr2_u0u1u2[beta];

                                    // m1 = n1;
                                    // left_ind = n0 * N * N + n1 * N + n2;
                                     // right_ind = m0 * N * N + m1 * N + m2;
                                    right_ind = right_ind_part1 + right_ind_part2 + m2;

                                    elem_left = u_left_ptr[left_ind];
                                     elem_right = u_right_ptr[right_ind];

                                     m0_double = static_cast<double>(m0);
                                     n0_double = static_cast<double>(n0);
                                    m2_double = static_cast<double>(m2);
                                    n2_double = static_cast<double>(n2);

                                    R_hat[0] = (m0_double - n0_double) / std::sqrt(2.0);
                                    R_hat[1] = 0.0;
                                    R_hat[2] = (m2_double - n2_double) / std::sqrt(2.0);

                                     J = (j4_val + std::sqrt(2.0) * (j3_val - j4_val) * std::abs(R_hat[alpha])) *
                                        delta(alpha, beta)
                                        + 2.0 * j5_val * R_hat[alpha] * R_hat[beta] * (1 - delta(alpha, beta));

                                    val2 += J * elem_left * elem_right;
                                } //end m2
                            } //end m0
                        } //end n2
                    } //end n1
                } //end n0
            } //end beta
        } //end alpha

        val2 *= 0.5;


        double val3 = 0; //term 2NN3
        int m2;

        for (int alpha = 0; alpha < 3; alpha++)
        {u_left_ptr = ptr2_u0u1u2[alpha];
            for (int beta = 0; beta < 3; beta++)
            {u_right_ptr = ptr2_u0u1u2[beta];
                for (int n0 = 0; n0 < N; n0++)
                {
                    left_ind_part1=n0*N_power_2;
                    for (int n1 = 0; n1 < N; n1++)
                    {
                        left_ind_part2=n1*N;
                        for (int n2 = 0; n2 < N; n2++)
                        {
                            left_ind =left_ind_part1 + left_ind_part2 + n2;
                            for (int m0 : {python_mod(n0 - 1, N), python_mod(n0 + 1, N)})
                            {
                                right_ind_part1=m0*N_power_2;
                                for (int m1 : {python_mod(n1 - 1, N), python_mod(n1 + 1, N)})
                                {
                                  // u_left_ptr = ptr2_u0u1u2[alpha];
                                     // u_right_ptr = ptr2_u0u1u2[beta];

                                    m2 = n2;
                                   // left_ind = n0 * N * N + n1 * N + n2;
                                   // right_ind = m0 * N * N + m1 * N + m2;
                                    right_ind = right_ind_part1 + m1 * N + m2;

                                     elem_left = u_left_ptr[left_ind];
                                     elem_right = u_right_ptr[right_ind];

                                     m0_double = static_cast<double>(m0);
                                     n0_double = static_cast<double>(n0);
                                     m1_double = static_cast<double>(m1);
                                    n1_double = static_cast<double>(n1);

                                    R_hat[0] = (m0_double - n0_double) / std::sqrt(2.0);
                                    R_hat[1] = (m1_double - n1_double) / std::sqrt(2.0);
                                    R_hat[2] = 0.0;

                                    double J = (j4_val + std::sqrt(2.0) * (j3_val - j4_val) * std::abs(R_hat[alpha])) *
                                        delta(alpha, beta)
                                        + 2.0 * j5_val * R_hat[alpha] * R_hat[beta] * (1 - delta(alpha, beta));


                                    val3 += J * elem_left * elem_right;
                                } //end m1
                            } //end m0
                        } //end n2
                    } //end n1
                } //end n0
            } //end beta
        } //end alpha

        val3 *= 0.5;

        return val1 + val2 + val3;
    } //end E_short_2NN


    double E_short_3NN()
    {
        double val = 0;
        // std::shared_ptr<double[]> u_left_ptr;
        // std::shared_ptr<double[]> u_right_ptr;
        int left_ind;
        int right_ind;
        double elem_left;
        double elem_right;
        double m0_double;
        double n0_double;
        double m1_double;
        double n1_double;
        double m2_double;
        double n2_double;
        double J ;
        int  left_ind_part1;
        int  left_ind_part2;
        int  right_ind_part1;
        int  right_ind_part2;


        for (int alpha = 0; alpha < 3; alpha++)
        {
            u_left_ptr = ptr2_u0u1u2[alpha];
            for (int beta = 0; beta < 3; beta++)
            {
                u_right_ptr = ptr2_u0u1u2[beta];
                for (int n0 = 0; n0 < N; n0++)
                {
                    left_ind_part1=n0*N_power_2;
                    for (int n1 = 0; n1 < N; n1++)
                    {
                        left_ind_part2=n1*N;
                        for (int n2 = 0; n2 < N; n2++)
                        {
                            left_ind = left_ind_part1 + left_ind_part2+ n2;
                            for (int m0 : {python_mod(n0 - 1, N), python_mod(n0 + 1, N)})
                            {
                                right_ind_part1=m0*N_power_2;
                                for (int m1 : {python_mod(n1 - 1, N), python_mod(n1 + 1, N)})
                                {
                                    right_ind_part2=m1*N;
                                    for (int m2 : {python_mod(n2 - 1, N), python_mod(n2 + 1, N)})
                                    {
                                         // u_left_ptr = ptr2_u0u1u2[alpha];
                                         // u_right_ptr = ptr2_u0u1u2[beta];
                                        // left_ind = n0 * N * N + n1 * N + n2;
                                         // right_ind = m0 * N * N + m1 * N + m2;
                                        right_ind = right_ind_part1 + right_ind_part2 + m2;

                                         elem_left = u_left_ptr[left_ind];
                                        elem_right = u_right_ptr[right_ind];

                                        m0_double = static_cast<double>(m0);
                                        n0_double = static_cast<double>(n0);
                                       m1_double = static_cast<double>(m1);
                                         n1_double = static_cast<double>(n1);
                                        m2_double = static_cast<double>(m2);
                                         n2_double = static_cast<double>(n2);

                                        R_hat[0] = (m0_double - n0_double) / std::sqrt(3.0);
                                        R_hat[1] = (m1_double - n1_double) / std::sqrt(3.0);
                                        R_hat[2] = (m2_double - n2_double) / std::sqrt(3.0);

                                         J = j6_val * delta(alpha, beta) + 3.0 * j7_val * R_hat[alpha] * R_hat[
                                            beta] * (1 - delta(alpha, beta));
                                        val += J * elem_left * elem_right;
                                    } //end m2
                                } //end m1
                            } //end m0
                        } //end n2
                    } //end n1
                } //end n0
            } //end beta
        } //end alpha
        val *= 0.5;

        return val;
    } //end E_short_3NN

    ///delta function
    double delta(const int& i, const int& j)
    {
        if (i == j)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }

    int python_mod(int a, int M)
    {
        return (a % M + M) % M;
    }


    //elastic energy

    double E_elas(const std::shared_ptr<double[]>& eta_H, const std::shared_ptr<double[]>& v0,
                  const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2)
    {
        //homogeneous part

        double eta_H1 = eta_H[0];
        double eta_H2 = eta_H[1];
        double eta_H3 = eta_H[2];

        double eta_H4 = eta_H[3];
        double eta_H5 = eta_H[4];
        double eta_H6 = eta_H[5];

        double N_double = static_cast<double>(N);

        double homogeneous_part = 0.5 * std::pow(N_double, 3) * B11_val * (std::pow(eta_H1, 2) + std::pow(eta_H2, 2) +
                std::pow(eta_H3, 2))
            + std::pow(N_double, 3) * B12_val * (eta_H1 * eta_H2 + eta_H2 * eta_H3 + eta_H3 * eta_H1)
            + 0.5 * std::pow(N_double, 3) * B44_val * (std::pow(eta_H4, 2) + std::pow(eta_H5, 2) + std::pow(eta_H6, 2));


        //index: Ba=0, Ti=1, O1=2, O2=3,O3=4

        //inhomogeneous part
        double E_elas_I1 = 0;
        int i_plus_1;
        int i_minus_1;
        int j_plus_1;
        int j_minus_1;

        for (int i = 0; i < N; i++)
        {
            i_plus_1 = python_mod(i + 1, N);
            i_minus_1 = python_mod(i - 1, N);
            for (int j = 0; j < N; j++)
            {
                j_plus_1 = python_mod(j + 1, N);
                j_minus_1 = python_mod(j - 1, N);
                for (int k = 0; k < N; k++)
                {
                     // i_plus_1 = python_mod(i + 1, N);
                     // i_minus_1 = python_mod(i - 1, N);
                    // j_plus_1 = python_mod(j + 1, N);
                    // j_minus_1 = python_mod(j - 1, N);

                    int ind_Ba = 0;

                    int ijk_Ba = flattened_ind_for_E_elas(i, j, k, ind_Ba);

                    int iPlus1jk_Ba = flattened_ind_for_E_elas(i_plus_1, j, k, ind_Ba);
                    int iMinus1jk_Ba = flattened_ind_for_E_elas(i_minus_1, j, k, ind_Ba);

                    int ijPlus1k_Ba = flattened_ind_for_E_elas(i, j_plus_1, k, ind_Ba);
                    int ijMinus1k_Ba = flattened_ind_for_E_elas(i, j_minus_1, k, ind_Ba);


                    //row 1
                    E_elas_I1 += gamma11_val * std::pow(v0[ijk_Ba] - v0[iPlus1jk_Ba], 2);

                    //row2
                    E_elas_I1 += gamma11_val * std::pow(v0[ijk_Ba] - v0[iMinus1jk_Ba], 2);

                    //row 3
                    E_elas_I1 += gamma12_val * (v0[ijk_Ba] - v0[iPlus1jk_Ba]) * (v1[ijk_Ba] - v1[ijPlus1k_Ba]);

                    //row 4
                    E_elas_I1 += gamma12_val * (v0[ijk_Ba] - v0[iMinus1jk_Ba]) * (v1[ijk_Ba] - v1[ijMinus1k_Ba]);


                    //row 5
                    E_elas_I1 += gamma44_val * std::pow(v0[ijk_Ba] - v0[ijPlus1k_Ba] + v1[ijk_Ba] - v1[iPlus1jk_Ba], 2);

                    //row 6
                    E_elas_I1 += gamma44_val * std::pow(v0[ijk_Ba] - v0[ijMinus1k_Ba] + v1[ijk_Ba] - v1[iMinus1jk_Ba],
                                                        2);
                } //end k
            } //end j
        } //end i

        double E_elas_I2 = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    int i_plus_1 = python_mod(i + 1, N);
                    int i_minus_1 = python_mod(i - 1, N);

                    int k_plus_1 = python_mod(k + 1, N);
                    int k_minus_1 = python_mod(k - 1, N);

                    int ind_Ba = 0;

                    int ijk_Ba = flattened_ind_for_E_elas(i, j, k, ind_Ba);

                    int ij_kPlus1 = flattened_ind_for_E_elas(i, j, k_plus_1, ind_Ba);
                    int ij_kMinus1 = flattened_ind_for_E_elas(i, j, k_minus_1, ind_Ba);

                    int iPlus1_jk = flattened_ind_for_E_elas(i_plus_1, j, k, ind_Ba);
                    int iMinus1_jk = flattened_ind_for_E_elas(i_minus_1, j, k, ind_Ba);

                    //row 1
                    E_elas_I2 += gamma11_val * std::pow(v2[ijk_Ba] - v2[ij_kPlus1], 2);

                    //row 2
                    E_elas_I2 += gamma11_val * std::pow(v2[ijk_Ba] - v2[ij_kMinus1], 2);

                    //row 3
                    E_elas_I2 += gamma12_val * (v2[ijk_Ba] - v2[ij_kPlus1]) * (v0[ijk_Ba] - v0[iPlus1_jk]);

                    //row 4
                    E_elas_I2 += gamma12_val * (v2[ijk_Ba] - v2[ij_kMinus1]) * (v0[ijk_Ba] - v0[iMinus1_jk]);

                    //row 5
                    E_elas_I2 += gamma44_val * std::pow(v2[ijk_Ba] - v2[iPlus1_jk] + v0[ijk_Ba] - v0[ij_kPlus1], 2);

                    //row 6
                    E_elas_I2 += gamma44_val * std::pow(v2[ijk_Ba] - v2[iMinus1_jk] + v0[ijk_Ba] - v0[ij_kMinus1], 2);
                } //end k
            } //end j
        } //end i
        double E_elas_I3 = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    int j_plus_1 = python_mod(j + 1, N);
                    int j_minus_1 = python_mod(j - 1, N);

                    int k_plus_1 = python_mod(k + 1, N);
                    int k_minus_1 = python_mod(k - 1, N);

                    int ind_Ba = 0;

                    int ijk_Ba = flattened_ind_for_E_elas(i, j, k, ind_Ba);

                    int i_jPlus1_k = flattened_ind_for_E_elas(i, j_plus_1, k, ind_Ba);
                    int i_jMinus1_k = flattened_ind_for_E_elas(i, j_minus_1, k, ind_Ba);

                    int ij_kPlus1 = flattened_ind_for_E_elas(i, j, k_plus_1, ind_Ba);
                    int ij_kMinus1 = flattened_ind_for_E_elas(i, j, k_minus_1, ind_Ba);

                    //row 1
                    E_elas_I3 += gamma11_val * std::pow(v1[ijk_Ba] - v1[i_jPlus1_k], 2);

                    //row 2
                    E_elas_I3 += gamma11_val * std::pow(v1[ijk_Ba] - v1[i_jMinus1_k], 2);

                    //row 3
                    E_elas_I3 += gamma12_val * (v1[ijk_Ba] - v1[i_jPlus1_k]) * (v2[ijk_Ba] - v2[ij_kPlus1]);

                    //row 4
                    E_elas_I3 += gamma12_val * (v1[ijk_Ba] - v1[i_jMinus1_k]) * (v2[ijk_Ba] - v2[ij_kMinus1]);

                    // row 5
                    E_elas_I3 += gamma44_val * std::pow(v1[ijk_Ba] - v1[ij_kPlus1] + v2[ijk_Ba] - v2[i_jPlus1_k], 2);

                    //row 6
                    E_elas_I3 += gamma44_val * std::pow(v1[ijk_Ba] - v1[ij_kMinus1] + v2[ijk_Ba] - v2[i_jMinus1_k], 2);
                } //end k
            } //end j
        } //end i

        return homogeneous_part + E_elas_I1 + E_elas_I2 + E_elas_I3;
    } //end E_elas

    int flattened_ind_for_E_elas(const int& i, const int& j, const int& k, const int& q)
    {
        //i,j,k take values from 0 to N-1, q takes values from 0 to 4
        int ind = q + 5 * k + 5 * j * N + 5 * i * N * N;

        return ind;
    }


    double E_elas_mode_int(const std::shared_ptr<double[]>& eta_H, const std::shared_ptr<double[]>& v0,
                           const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2)
    {
        double val = 0;
        //homogeneous part

        double eta_H1 = eta_H[0];
        double eta_H2 = eta_H[1];
        double eta_H3 = eta_H[2];

        double eta_H4 = eta_H[3];
        double eta_H5 = eta_H[4];
        double eta_H6 = eta_H[5];

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    int i_minus_1 = python_mod(i - 1, N);
                    int j_minus_1 = python_mod(j - 1, N);
                    int k_minus_1 = python_mod(k - 1, N);
                    int ind_Ba = 0;

                    //i,j,k
                    int ijk = flattened_ind_for_E_elas(i, j, k, ind_Ba);

                    //i-1,j,k
                    int iMinus1_jk = flattened_ind_for_E_elas_mode_int(i_minus_1, j, k, ind_Ba);


                    //i-1,j-1,k
                    int iMinus1_jMinus1_k = flattened_ind_for_E_elas(i_minus_1, j_minus_1, k, ind_Ba);

                    //i,j-1,k
                    int i_jMinus1_k = flattened_ind_for_E_elas(i, j_minus_1, k, ind_Ba);

                    //i-1,j,k-1
                    int iMinus1_j_kMinus1 = flattened_ind_for_E_elas(i_minus_1, j, k_minus_1, ind_Ba);

                    //i,j,k-1
                    int ij_kMinus1 = flattened_ind_for_E_elas(i, j, k_minus_1, ind_Ba);

                    //i-1,j-1,k-1
                    int iMinus1_jMinus1_kMinus1 = flattened_ind_for_E_elas(i_minus_1, j_minus_1, k_minus_1, ind_Ba);

                    //i,j-1,k-1
                    int i_jMinus1_kMinus1 = flattened_ind_for_E_elas(i, j_minus_1, k_minus_1, ind_Ba);

                    double Delta_vxx_ijk = v0[iMinus1_jk] - v0[ijk]
                        + v0[iMinus1_jMinus1_k] - v0[i_jMinus1_k]
                        + v0[iMinus1_j_kMinus1] - v0[ij_kMinus1]
                        + v0[iMinus1_jMinus1_kMinus1] - v0[i_jMinus1_kMinus1];

                    double etaI1ijk = Delta_vxx_ijk / 4.0;

                    double Delta_vyy_ijk = v1[i_jMinus1_k] - v1[ijk]
                        + v1[i_jMinus1_kMinus1] - v1[ij_kMinus1]
                        + v1[iMinus1_jMinus1_k] - v1[iMinus1_jk]
                        + v1[iMinus1_jMinus1_kMinus1] - v1[iMinus1_j_kMinus1];
                    double etaI2ijk = Delta_vyy_ijk / 4.0;


                    double Delta_vzz_ijk = v2[ij_kMinus1] - v2[ijk]
                        + v2[iMinus1_j_kMinus1] - v2[iMinus1_jk]
                        + v2[i_jMinus1_kMinus1] - v2[i_jMinus1_k]
                        + v2[iMinus1_jMinus1_kMinus1] - v2[iMinus1_jMinus1_k];

                    double etaI3ijk = Delta_vzz_ijk / 4.0;

                    double Delta_vxy_ijk = v1[iMinus1_jk] - v1[ijk]
                        + v1[iMinus1_jMinus1_k] - v1[i_jMinus1_k]
                        + v1[iMinus1_j_kMinus1] - v1[ij_kMinus1]
                        + v1[iMinus1_jMinus1_kMinus1] - v1[iMinus1_jMinus1_k];


                    double Delta_vyx_ijk = v0[i_jMinus1_k] - v0[ijk]
                        + v0[iMinus1_jMinus1_k] - v0[iMinus1_jk]
                        + v0[i_jMinus1_kMinus1] - v0[ij_kMinus1]
                        + v0[iMinus1_jMinus1_kMinus1] - v0[iMinus1_jMinus1_k];

                    double etaI6ijk = (Delta_vxy_ijk + Delta_vyx_ijk) / 4.0;


                    double Delta_vzx_ijk = v0[ij_kMinus1] - v0[ijk]
                        + v0[iMinus1_j_kMinus1] - v0[iMinus1_jk]
                        + v0[i_jMinus1_kMinus1] - v0[i_jMinus1_k]
                        + v0[iMinus1_jMinus1_kMinus1] - v0[iMinus1_j_kMinus1];


                    double Delta_vxz_ijk = v2[iMinus1_jk] - v2[ijk]
                        + v2[iMinus1_j_kMinus1] - v2[ij_kMinus1]
                        + v2[iMinus1_jMinus1_k] - v2[i_jMinus1_k]
                        + v2[iMinus1_jMinus1_kMinus1] - v2[iMinus1_j_kMinus1];

                    double etaI5ijk = (Delta_vzx_ijk + Delta_vxz_ijk) / 4.0;


                    double Delta_vyz_ijk = v2[i_jMinus1_k] - v2[ijk]
                        + v2[i_jMinus1_kMinus1] - v2[ij_kMinus1]
                        + v2[iMinus1_jMinus1_k] - v2[iMinus1_jk]
                        + v2[iMinus1_jMinus1_kMinus1] - v2[i_jMinus1_kMinus1];


                    double Delta_vzy_ijk = v1[ij_kMinus1] - v1[ijk]
                        + v1[i_jMinus1_kMinus1] - v1[i_jMinus1_k]
                        + v1[iMinus1_j_kMinus1] - v1[iMinus1_jk]
                        + v1[iMinus1_jMinus1_kMinus1] - v1[i_jMinus1_kMinus1];

                    double etaI4ijk = (Delta_vyz_ijk + Delta_vzy_ijk) / 4.0;


                    int u_ind_ijk = i * N * N + j * N + k;

                    double E_int_ijk = 0;
                    double eta1ijk=eta_H1+etaI1ijk;
                    double eta2ijk=eta_H2+etaI2ijk;
                    double eta3ijk=eta_H3+etaI3ijk;
                    double eta4ijk=eta_H4+etaI4ijk;
                    double eta5ijk=eta_H5+etaI5ijk;
                    double eta6ijk=eta_H6+etaI6ijk;

                    E_int_ijk += B100_val * eta1ijk * u0[u_ind_ijk] * u0[u_ind_ijk];

                    E_int_ijk += B111_val * eta1ijk * u1[u_ind_ijk] * u1[u_ind_ijk];

                    E_int_ijk += B122_val * eta1ijk * u2[u_ind_ijk] * u2[u_ind_ijk];

                    E_int_ijk += B200_val * eta2ijk * u0[u_ind_ijk] * u0[u_ind_ijk];

                    E_int_ijk += B211_val * eta2ijk * u1[u_ind_ijk] * u1[u_ind_ijk];

                    E_int_ijk += B222_val * eta2ijk * u2[u_ind_ijk] * u2[u_ind_ijk];

                    E_int_ijk += B300_val * eta3ijk * u0[u_ind_ijk] * u0[u_ind_ijk];

                    E_int_ijk += B311_val * eta3ijk * u1[u_ind_ijk] * u1[u_ind_ijk];

                    E_int_ijk += B322_val * eta3ijk * u2[u_ind_ijk] * u2[u_ind_ijk];

                    E_int_ijk += B412_val * eta4ijk * u1[u_ind_ijk] * u2[u_ind_ijk];

                    E_int_ijk += B421_val * eta4ijk * u2[u_ind_ijk] * u1[u_ind_ijk];

                    E_int_ijk += B502_val * eta5ijk * u0[u_ind_ijk] * u2[u_ind_ijk];

                    E_int_ijk += B520_val * eta5ijk * u2[u_ind_ijk] * u0[u_ind_ijk];

                    E_int_ijk += B601_val * etaI6ijk * u0[u_ind_ijk] * u1[u_ind_ijk];

                    E_int_ijk += B610_val * etaI6ijk * u1[u_ind_ijk] * u0[u_ind_ijk];


                    val += 0.5 * E_int_ijk;
                } //end k
            } //end j
        } //end i

        return val;
    } //end E_elas_mode_int

    int flattened_ind_for_E_elas_mode_int(const int& i, const int& j, const int& k, const int& q)
    {
        //i,j,k take values from 0 to N1, q takes values from 0 to 4
        int ind = q + 5 * k + 5 * j * N + 5 * i * N * N;

        return ind;
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

    double coef_Ba;
    double coef_Ti;
    double coef_O_parallel;
    double coef_O_perpendicular;

    std::shared_ptr<double[]> u0;
    std::shared_ptr<double[]> u1;
    std::shared_ptr<double[]> u2;
    std::shared_ptr<double[]> Q;
    std::vector<std::shared_ptr<double[]>> ptr2_u0u1u2;

    std::shared_ptr<double[]> R_hat;
    std::shared_ptr<double[]> u_left_ptr;
    std::shared_ptr<double[]> u_right_ptr;
};



std::shared_ptr<potentialFunction> createPotentialFunction(const std::string& funcName, const std::string& coefsJsonStr)
{
    if (funcName == "V_BaTiO3")
    {
        return std::make_shared<V_BaTiO3>(coefsJsonStr);
    }
    if(funcName=="V_BaTiO3_parallel")
    {
        return std::make_shared<V_BaTiO3_parallel>(coefsJsonStr);
    }

    else
    {
        throw std::invalid_argument("Unknown potential function type");
    }
}






/////////////////////parallel code
void V_BaTiO3_parallel::v2u(const std::shared_ptr<double[]>& v, std::shared_ptr<double[]>& u)
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







 void V_BaTiO3_parallel::init()
    {
        this->json2Coefs(coefsInStr);

        this->gamma11_val = B11_val / 4.0;
        this->gamma12_val = B12_val / 8.0;
        this->gamma44_val = B44_val / 8.0;

        B211_val = B100_val;
        B322_val = B100_val;

        B122_val = B111_val;
        B200_val = B111_val;
        B222_val = B111_val;
        B300_val = B111_val;
        B311_val = B111_val;

        B421_val = B412_val;
        B502_val = B412_val;
        B520_val = B412_val;
        B601_val = B412_val;
        B610_val = B412_val;
        this->N_power_6=static_cast<uint64_t>(N)*N*N*N*N*N;
        this->N_power_5=static_cast<uint64_t>(N)*N*N*N*N;
        this->N_power_4=static_cast<uint64_t>(N)*N*N*N;
        this->N_power_3=static_cast<uint64_t>(N)*N*N;
        this->N_power_2=static_cast<uint64_t>(N)*N;
        this->_3_power_2=3*3;

        u0 = std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());
        u1 = std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());
        u2 = std::shared_ptr<double[]>(new double[elemNumTot_u], std::default_delete<double[]>());

        // ptr2_u0u1u2 = std::shared_ptr<std::shared_ptr<double[]>[]>(new std::shared_ptr<double[]>[3]{u0, u1, u2});
        ptr2_u0u1u2 = {u0, u1, u2};
        R_hat = std::shared_ptr<double[]>(new double[3], std::default_delete<double[]>());

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
            << ", gamma11_val=" << gamma11_val << ", gamma12_val=" << gamma12_val << ", gamma44_val=" << gamma44_val
            << ", B100_val=" << B100_val << ", B211_val=" << B211_val << ", B322_val" << B322_val
            << ", B111_val=" << B111_val << ", B122_val=" << B122_val << ", B200_val" << B200_val
            << ", B222_val=" << B222_val << ", B300_val=" << B300_val << ", B311_val" << B311_val
            << ", B412_val=" << B412_val << ", B421_val=" << B421_val << ", B502_val=" << B502_val
            << ", B520_val=" << B520_val << ", B601_val=" << B601_val << ", B610_val=" << B610_val
            << std::endl;
    }



void V_BaTiO3_parallel::json2Coefs(const std::string& coefsStr)
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

        //read B100
        if (std::getline(iss, temp, ','))
        {
            this->B100_val = std::stod(temp);
        }

        //read B111
        if (std::getline(iss, temp, ','))
        {
            this->B111_val = std::stod(temp);
        }

        //read B412
        if (std::getline(iss, temp, ','))
        {
            this->B412_val = std::stod(temp);
        }
    } // end json2Coefs



 void V_BaTiO3_parallel::fill_Q()
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



double V_BaTiO3_parallel::operator()(const std::shared_ptr<double[]>& eta_H,const std::shared_ptr<double[]>& v0,const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2)
{
    this->v2u(v0, u0);
    this->v2u(v1, u1);
    this->v2u(v2, u2);
    // std::cout<<"Q[1]="<<Q[1]<<std::endl;
    double energy_self = this->E_self(u0, u1, u2);
    // std::cout << "energy_self=" << energy_self << std::endl;

    double energy_dipole = E_dpl();
    // std::cout << "energy_dipole=" << energy_dipole << std::endl;

    // double energy_short = E_short();
    // std::cout << "energy_short=" << energy_short << std::endl;

    // double energy_elas = E_elas(eta_H, v0, v1, v2);
    // std::cout << "energy_elas=" << energy_elas << std::endl;

    // double energy_elas_mode_int=E_elas_mode_int(eta_H,v0,v1,v2);
    // std::cout<<"energy_elas_mode_int="<<energy_elas_mode_int<<std::endl;

    double pot_energy=energy_self+energy_dipole;//+energy_short+energy_elas+energy_elas_mode_int;
    return pot_energy;



}






double V_BaTiO3_parallel::E_self(const std::shared_ptr<double[]>& u0, const std::shared_ptr<double[]>& u1,
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


}//end E_self


// Function to compute indices from a flattened index
void V_BaTiO3_parallel::compute_indices(uint64_t idx, int N, int& i1, int& j1, int& k1, int& i2, int& j2, int& k2) {
    k2 = idx % N;
    idx /= N;
    j2 = idx % N;
    idx /= N;
    i2 = idx % N;
    idx /= N;
    k1 = idx % N;
    idx /= N;
    j1 = idx % N;
    idx /= N;
    i1 = idx % N;
}



// Dipole energy computation using smart pointers and C++11 threads
/// Function to calculate a partial sum for a range of indices
void V_BaTiO3_parallel::compute_partial_sum(uint64_t start_idx, uint64_t end_idx, int N, double& partial_val, int alpha, int beta) {
    partial_val = 0.0;

    int Q_elem_ind_part1, Q_elem_ind_part2, Q_elem_ind_part3, Q_elem_ind_part4;
    int Q_elem_ind_part5, Q_elem_ind_part6, Q_elem_ind_part7, Q_elem_ind;
    int u_left_ind_part1, u_left_ind_part2, u_right_ind_part1, u_right_ind_part2;
    int u_left_ind, u_right_ind;

    double Q_elem_val, u_left_elem_val, u_right_elem_val;

    Q_elem_ind_part7 = alpha * 3;
    u_left_ptr = ptr2_u0u1u2[alpha];
    u_right_ptr = ptr2_u0u1u2[beta];

    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {
        int i1, j1, k1, i2, j2, k2;
        compute_indices(idx, N, i1, j1, k1, i2, j2, k2);

        // Compute parts of the index for Q and u_left, u_right
        Q_elem_ind_part1 = i1 * (N_power_5 * _3_power_2);
        Q_elem_ind_part2 = j1 * (N_power_4 * _3_power_2);
        Q_elem_ind_part3 = k1 * (N_power_3 * _3_power_2);
        Q_elem_ind_part4 = i2 * (N_power_2 * _3_power_2);
        Q_elem_ind_part5 = j2 * (N * _3_power_2);
        Q_elem_ind_part6 = k2 * (_3_power_2);
        Q_elem_ind = Q_elem_ind_part1 + Q_elem_ind_part2 + Q_elem_ind_part3 +
                     Q_elem_ind_part4 + Q_elem_ind_part5 + Q_elem_ind_part6 +
                     Q_elem_ind_part7 + beta;

        u_left_ind_part1 = i1 * N_power_2;
        u_left_ind_part2 = j1 * N;
        u_left_ind = u_left_ind_part1 + u_left_ind_part2 + k1;

        u_right_ind_part1 = i2 * N_power_2;
        u_right_ind_part2 = j2 * N;
        u_right_ind = u_right_ind_part1 + u_right_ind_part2 + k2;

        // Access values from Q and u_left, u_right
        Q_elem_val = Q[Q_elem_ind];
        u_left_elem_val = u_left_ptr[u_left_ind];
        u_right_elem_val = u_right_ptr[u_right_ind];

        // Accumulate the result
        partial_val += Q_elem_val * u_left_elem_val * u_right_elem_val;
    }
}


/// Main dipole energy function using std::thread
double V_BaTiO3_parallel::E_dpl()
{
    double val = 0.0;
    // uint64_t N_power_6 = N * N * N * N * N * N;

    // Number of threads to use
    const unsigned num_threads = std::thread::hardware_concurrency();
    std::vector<double> partial_vals(num_threads, 0.0);

    // Outer loops over alpha and beta
    for (int alpha = 0; alpha < 3; alpha++)
    {
        for (int beta = 0; beta < 3; beta++)
        {
            // Clear threads before usage
            threads.clear();

            // Divide the total range (0 to N^6) among threads
            uint64_t chunk_size = N_power_6 / num_threads;

            for (unsigned t = 0; t < num_threads; ++t)
            {
                uint64_t start_idx = t * chunk_size;
                uint64_t end_idx = (t == num_threads - 1) ? N_power_6 : (t + 1) * chunk_size;

                // Create threads and assign tasks
                threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum, this, start_idx, end_idx, N, std::ref(partial_vals[t]), alpha, beta));
            }

            // Join threads and collect results
            for (auto& thread : threads)
            {
                thread.join();
            }

            // Accumulate the results from each thread
            for (const auto& partial_val : partial_vals)
            {
                val += partial_val;
            }

            // Clear threads after usage
            threads.clear();
        }
    }

    return val;
}

// E_short
double V_BaTiO3_parallel::E_short()
{
    double energy_short_1NN = E_short_1NN();
    // std::cout<<"energy_short_1NN="<<energy_short_1NN<<std::endl;

    double energy_short_2NN = E_short_2NN();
    // std::cout<<"energy_short_2NN="<<energy_short_2NN<<std::endl;

    double energy_short_3NN = E_short_3NN();
    // std::cout<<"energy_short_3NN="<<energy_short_3NN<<std::endl;

    return energy_short_1NN + energy_short_2NN + energy_short_3NN;

}

/// Function to compute indices from a flattened index
void V_BaTiO3_parallel::compute_indices_1NN(uint64_t idx, int& n0, int& n1, int& n2) {
    n2 = idx % N;
    idx /= N;
    n1 = idx % N;
    idx /= N;
    n0 = idx % N;
}

/// Function to calculate a partial sum for the 1NN1 term (shift in n2)
void V_BaTiO3_parallel::compute_partial_sum_1NN1(uint64_t start_idx, uint64_t end_idx, double& partial_val, int alpha) {
    partial_val = 0.0;

    int ind_u_left, ind_u_right;
    int ind_u_left_part1, ind_u_left_part2;
    double elem_left, elem_right;

    u_left_ptr = ptr2_u0u1u2[alpha];
    u_right_ptr = ptr2_u0u1u2[alpha];

    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {
        int n0, n1, n2;
        compute_indices_1NN(idx, n0, n1, n2);

        ind_u_left_part1 = n0 * N_power_2;
        ind_u_left_part2 = n1 * N;

        ind_u_left = ind_u_left_part1 + ind_u_left_part2 + n2;
        ind_u_right = ind_u_left_part1 + ind_u_left_part2 + (n2 + 1) % N;

        // Get the left and right elements
        elem_left = u_left_ptr[ind_u_left];
        elem_right = u_right_ptr[ind_u_right];

        // Accumulate the value
        partial_val += j2_val * elem_left * elem_right;
    }
}

/// Function to calculate a partial sum for the 1NN2 term (shift in n1)
void V_BaTiO3_parallel::compute_partial_sum_1NN2(uint64_t start_idx, uint64_t end_idx, double& partial_val, int alpha) {
    partial_val = 0.0;

    int ind_u_left, ind_u_right;
    int ind_u_left_part1, ind_u_left_part2, ind_u_right_part2;
    double elem_left, elem_right;

    u_left_ptr = ptr2_u0u1u2[alpha];
    u_right_ptr = ptr2_u0u1u2[alpha];

    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {
        int n0, n1, n2;
        compute_indices_1NN(idx, n0, n1, n2);

        ind_u_left_part1 = n0 * N_power_2;
        ind_u_left_part2 = n1 * N;
        ind_u_right_part2 = ((n1 + 1) % N) * N;

        ind_u_left = ind_u_left_part1 + ind_u_left_part2 + n2;
        ind_u_right = ind_u_left_part1 + ind_u_right_part2 + n2;

        // Get the left and right elements
        elem_left = u_left_ptr[ind_u_left];
        elem_right = u_right_ptr[ind_u_right];

        // Accumulate the value
        partial_val += j2_val * elem_left * elem_right;
    }
}

/// Function to calculate a partial sum for the 1NN3 term (shift in n0)
void V_BaTiO3_parallel::compute_partial_sum_1NN3(uint64_t start_idx, uint64_t end_idx, double& partial_val, int alpha) {
    partial_val = 0.0;

    int ind_u_left, ind_u_right;
    int ind_u_left_part1, ind_u_left_part2, ind_u_right_part1;
    double elem_left, elem_right;

    u_left_ptr = ptr2_u0u1u2[alpha];
    u_right_ptr = ptr2_u0u1u2[alpha];

    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {
        int n0, n1, n2;
        compute_indices_1NN(idx, n0, n1, n2);

        ind_u_left_part1 = n0 * N_power_2;
        ind_u_right_part1 = ((n0 + 1) % N) * N_power_2;
        ind_u_left_part2 = n1 * N;

        ind_u_left = ind_u_left_part1 + ind_u_left_part2 + n2;
        ind_u_right = ind_u_right_part1 + ind_u_left_part2 + n2;

        // Get the left and right elements
        elem_left = u_left_ptr[ind_u_left];
        elem_right = u_right_ptr[ind_u_right];

        // Accumulate the value
        partial_val += j2_val * elem_left * elem_right;
    }
}


/// Main short-range energy 1NN function using std::thread
double V_BaTiO3_parallel::E_short_1NN()
{
    double val1 = 0.0, val2 = 0.0, val3 = 0.0;

    // Number of threads to use
    const unsigned num_threads = std::thread::hardware_concurrency();
    std::vector<double> partial_vals(num_threads, 0.0);

    // Loop over alpha
    for (int alpha = 0; alpha < 3; alpha++)
    {
        // Clear threads before usage
        threads.clear();

        // Divide the total range (0 to N^3) among threads
        uint64_t chunk_size = N_power_3 / num_threads;

        // First nearest neighbor term (1NN1)
        for (unsigned t = 0; t < num_threads; ++t)
        {
            uint64_t start_idx = t * chunk_size;
            uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

            threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_1NN1, this, start_idx, end_idx, std::ref(partial_vals[t]), alpha));
        }
        for (auto& thread : threads) thread.join();
        for (const auto& partial_val : partial_vals) val1 += partial_val;

        // Clear threads before next term
        threads.clear();

        // Second nearest neighbor term (1NN2)
        for (unsigned t = 0; t < num_threads; ++t)
        {
            uint64_t start_idx = t * chunk_size;
            uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

            threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_1NN2, this, start_idx, end_idx, std::ref(partial_vals[t]), alpha));
        }
        for (auto& thread : threads) thread.join();
        for (const auto& partial_val : partial_vals) val2 += partial_val;

        // Clear threads before next term
        threads.clear();

        // Third nearest neighbor term (1NN3)
        for (unsigned t = 0; t < num_threads; ++t)
        {
            uint64_t start_idx = t * chunk_size;
            uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

            threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_1NN3, this, start_idx, end_idx, std::ref(partial_vals[t]), alpha));
        }
        for (auto& thread : threads) thread.join();
        for (const auto& partial_val : partial_vals) val3 += partial_val;

        // Clear threads after usage
        threads.clear();
    }

    return val1 + val2 + val3;
}

/// Function to compute indices from a flattened index for 2NN term
void V_BaTiO3_parallel::compute_indices_2NN(uint64_t idx, int& n0, int& n1, int& n2) {
    n2 = idx % N;
    idx /= N;
    n1 = idx % N;
    idx /= N;
    n0 = idx % N;
}


/// Function to calculate a partial sum for 2NN1 term
void V_BaTiO3_parallel::compute_partial_sum_2NN1(uint64_t start_idx, uint64_t end_idx, double& partial_val, int alpha, int beta) {
    partial_val = 0.0;

    int m0, left_ind, right_ind;
    int left_ind_part1, left_ind_part2, right_ind_part1, right_ind_part2;
    double elem_left, elem_right;
    double m1_double, n1_double, m2_double, n2_double, J;

    u_left_ptr = ptr2_u0u1u2[alpha];
    u_right_ptr = ptr2_u0u1u2[beta];

    for (uint64_t idx = start_idx; idx < end_idx; ++idx) {
        int n0, n1, n2;
        compute_indices_2NN(idx, n0, n1, n2);

        m0 = n0;
        left_ind_part1 = n0 * N_power_2;
        right_ind_part1 = m0 * N_power_2;
        left_ind_part2 = n1 * N;

        left_ind = left_ind_part1 + left_ind_part2 + n2;

        for (int m1 : {python_mod(n1 - 1, N), python_mod(n1 + 1, N)}) {
            right_ind_part2 = m1 * N;
            for (int m2 : {python_mod(n2 - 1, N), python_mod(n2 + 1, N)}) {
                right_ind = right_ind_part1 + right_ind_part2 + m2;

                elem_left = u_left_ptr[left_ind];
                elem_right = u_right_ptr[right_ind];

                R_hat[0] = 0.0;
                m1_double = static_cast<double>(m1);
                n1_double = static_cast<double>(n1);
                m2_double = static_cast<double>(m2);
                n2_double = static_cast<double>(n2);

                R_hat[1] = (m1_double - n1_double) / std::sqrt(2.0);
                R_hat[2] = (m2_double - n2_double) / std::sqrt(2.0);

                J = (j4_val + std::sqrt(2.0) * (j3_val - j4_val) * std::abs(R_hat[alpha])) *
                    delta(alpha, beta) +
                    2.0 * j5_val * R_hat[alpha] * R_hat[beta] * (1 - delta(alpha, beta));

                partial_val += J * elem_left * elem_right;
            }
        }
    }
}


void V_BaTiO3_parallel::compute_partial_sum_2NN2(uint64_t start_idx, uint64_t end_idx, double& partial_val, int alpha, int beta)
{
    partial_val = 0.0;

    int  m1, left_ind, right_ind;
    int left_ind_part1, left_ind_part2, right_ind_part1, right_ind_part2;
    double elem_left, elem_right;
    double m0_double, n0_double, m2_double, n2_double, J;

    u_left_ptr = ptr2_u0u1u2[alpha];
    u_right_ptr = ptr2_u0u1u2[beta];

    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {
        int n0, n1, n2;
        compute_indices_2NN(idx, n0, n1, n2);
        left_ind_part1=n0*N_power_2;
        m1 = n1;
        left_ind_part2=n1*N;
        right_ind_part2=m1*N;
        left_ind = left_ind_part1+ left_ind_part2 + n2;
        for (int m0 : {python_mod(n0 - 1, N), python_mod(n0 + 1, N)})
        {
            right_ind_part1=m0*N_power_2;
            for (int m2 : {python_mod(n2 - 1, N), python_mod(n2 + 1, N)})
            {
                // u_left_ptr = ptr2_u0u1u2[alpha];
                // u_right_ptr = ptr2_u0u1u2[beta];

                // m1 = n1;
                // left_ind = n0 * N * N + n1 * N + n2;
                // right_ind = m0 * N * N + m1 * N + m2;
                right_ind = right_ind_part1 + right_ind_part2 + m2;

                elem_left = u_left_ptr[left_ind];
                elem_right = u_right_ptr[right_ind];

                m0_double = static_cast<double>(m0);
                n0_double = static_cast<double>(n0);
                m2_double = static_cast<double>(m2);
                n2_double = static_cast<double>(n2);

                R_hat[0] = (m0_double - n0_double) / std::sqrt(2.0);
                R_hat[1] = 0.0;
                R_hat[2] = (m2_double - n2_double) / std::sqrt(2.0);

                J = (j4_val + std::sqrt(2.0) * (j3_val - j4_val) * std::abs(R_hat[alpha])) *
                   delta(alpha, beta)
                   + 2.0 * j5_val * R_hat[alpha] * R_hat[beta] * (1 - delta(alpha, beta));

                partial_val += J * elem_left * elem_right;
            } //end m2
        } //end m0

    }//end idx

}


void V_BaTiO3_parallel::compute_partial_sum_2NN3(uint64_t start_idx, uint64_t end_idx, double& partial_val, int alpha, int beta)
{

    partial_val = 0.0;

    int m2, left_ind, right_ind;
    int left_ind_part1, left_ind_part2, right_ind_part1;
    double elem_left, elem_right;
    double m0_double, n0_double, m1_double, n1_double, J;

    u_left_ptr = ptr2_u0u1u2[alpha];
    u_right_ptr = ptr2_u0u1u2[beta];
    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {
        int n0, n1, n2;
        compute_indices_2NN(idx, n0, n1, n2);
        left_ind_part1=n0*N_power_2;
        left_ind_part2=n1*N;
        left_ind =left_ind_part1 + left_ind_part2 + n2;
        for (int m0 : {python_mod(n0 - 1, N), python_mod(n0 + 1, N)})
        {
            right_ind_part1=m0*N_power_2;
            for (int m1 : {python_mod(n1 - 1, N), python_mod(n1 + 1, N)})
            {
                // u_left_ptr = ptr2_u0u1u2[alpha];
                // u_right_ptr = ptr2_u0u1u2[beta];

                m2 = n2;
                // left_ind = n0 * N * N + n1 * N + n2;
                // right_ind = m0 * N * N + m1 * N + m2;
                right_ind = right_ind_part1 + m1 * N + m2;

                elem_left = u_left_ptr[left_ind];
                elem_right = u_right_ptr[right_ind];

                m0_double = static_cast<double>(m0);
                n0_double = static_cast<double>(n0);
                m1_double = static_cast<double>(m1);
                n1_double = static_cast<double>(n1);

                R_hat[0] = (m0_double - n0_double) / std::sqrt(2.0);
                R_hat[1] = (m1_double - n1_double) / std::sqrt(2.0);
                R_hat[2] = 0.0;

                 J = (j4_val + std::sqrt(2.0) * (j3_val - j4_val) * std::abs(R_hat[alpha])) *
                    delta(alpha, beta)
                    + 2.0 * j5_val * R_hat[alpha] * R_hat[beta] * (1 - delta(alpha, beta));


                partial_val += J * elem_left * elem_right;
            } //end m1
        } //end m0

    }//end idx

}


/// Main short-range energy 2NN function using std::thread
double V_BaTiO3_parallel::E_short_2NN() {
    double val1 = 0.0, val2 = 0.0, val3 = 0.0;

    // Number of threads to use
    const unsigned num_threads = std::thread::hardware_concurrency();
    std::vector<double> partial_vals(num_threads, 0.0);

    // Loop over alpha and beta
    for (int alpha = 0; alpha < 3; alpha++) {
        for (int beta = 0; beta < 3; beta++) {

            // Clear threads before usage
            threads.clear();

            // Divide the total range (0 to N^3) among threads
            uint64_t chunk_size = N_power_3 / num_threads;

            // First nearest neighbor term (2NN1)
            for (unsigned t = 0; t < num_threads; ++t) {
                uint64_t start_idx = t * chunk_size;
                uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

                threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_2NN1, this, start_idx, end_idx, std::ref(partial_vals[t]), alpha, beta));
            }
            for (auto& thread : threads) thread.join();
            for (const auto& partial_val : partial_vals) val1 += partial_val;

            // Clear threads before next term
            threads.clear();

            // Second nearest neighbor term (2NN2)
            for (unsigned t = 0; t < num_threads; ++t) {
                uint64_t start_idx = t * chunk_size;
                uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

                threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_2NN2, this, start_idx, end_idx, std::ref(partial_vals[t]), alpha, beta));
            }
            for (auto& thread : threads) thread.join();
            for (const auto& partial_val : partial_vals) val2 += partial_val;

            // Clear threads before next term
            threads.clear();

            // Third nearest neighbor term (2NN3)
            for (unsigned t = 0; t < num_threads; ++t) {
                uint64_t start_idx = t * chunk_size;
                uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

                threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_2NN3, this, start_idx, end_idx, std::ref(partial_vals[t]), alpha, beta));
            }
            for (auto& thread : threads) thread.join();
            for (const auto& partial_val : partial_vals) val3 += partial_val;

            // Clear threads after usage
            threads.clear();
        }
    }

    // Combine values
    val1 *= 0.5;
    val2 *= 0.5;
    val3 *= 0.5;

    return val1 + val2 + val3;
}

int V_BaTiO3_parallel::python_mod(int a, int M)
{
    return (a % M + M) % M;
}

///delta function
double V_BaTiO3_parallel::delta(const int& i, const int& j)
{
    if (i == j)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}



void V_BaTiO3_parallel::compute_indices_3NN(uint64_t idx, int& n0, int& n1, int& n2) {
    n2 = idx % N;
    idx /= N;
    n1 = idx % N;
    idx /= N;
    n0 = idx % N;
}

void V_BaTiO3_parallel::compute_partial_sum_3NN(uint64_t start_idx, uint64_t end_idx, double& partial_val, int alpha, int beta)
{
    partial_val = 0.0;

    int left_ind, right_ind;
    int left_ind_part1, left_ind_part2, right_ind_part1, right_ind_part2;
    double elem_left, elem_right;
    double m0_double, n0_double, m1_double, n1_double, m2_double, n2_double, J;

    u_left_ptr = ptr2_u0u1u2[alpha];
    u_right_ptr = ptr2_u0u1u2[beta];
    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {
        int n0, n1, n2;
        compute_indices_3NN(idx, n0, n1, n2);
        left_ind_part1=n0*N_power_2;
        left_ind_part2=n1*N;
        left_ind = left_ind_part1 + left_ind_part2+ n2;
        for (int m0 : {python_mod(n0 - 1, N), python_mod(n0 + 1, N)})
        {
            right_ind_part1=m0*N_power_2;
            for (int m1 : {python_mod(n1 - 1, N), python_mod(n1 + 1, N)})
            {
                right_ind_part2=m1*N;
                for (int m2 : {python_mod(n2 - 1, N), python_mod(n2 + 1, N)})
                {
                    // u_left_ptr = ptr2_u0u1u2[alpha];
                    // u_right_ptr = ptr2_u0u1u2[beta];
                    // left_ind = n0 * N * N + n1 * N + n2;
                    // right_ind = m0 * N * N + m1 * N + m2;
                    right_ind = right_ind_part1 + right_ind_part2 + m2;

                    elem_left = u_left_ptr[left_ind];
                    elem_right = u_right_ptr[right_ind];

                    m0_double = static_cast<double>(m0);
                    n0_double = static_cast<double>(n0);
                    m1_double = static_cast<double>(m1);
                    n1_double = static_cast<double>(n1);
                    m2_double = static_cast<double>(m2);
                    n2_double = static_cast<double>(n2);

                    R_hat[0] = (m0_double - n0_double) / std::sqrt(3.0);
                    R_hat[1] = (m1_double - n1_double) / std::sqrt(3.0);
                    R_hat[2] = (m2_double - n2_double) / std::sqrt(3.0);

                    J = j6_val * delta(alpha, beta) + 3.0 * j7_val * R_hat[alpha] * R_hat[
                       beta] * (1 - delta(alpha, beta));
                    partial_val += J * elem_left * elem_right;
                } //end m2
            } //end m1
        } //end m0
    }//end idx
}


/// Main short-range energy 3NN function using std::thread
double V_BaTiO3_parallel::E_short_3NN() {
    double val = 0.0;

    // Number of threads to use
    const unsigned num_threads = std::thread::hardware_concurrency();
    std::vector<double> partial_vals(num_threads, 0.0);

    // Loop over alpha and beta
    for (int alpha = 0; alpha < 3; alpha++) {
        for (int beta = 0; beta < 3; beta++) {

            // Clear threads before usage
            threads.clear();

            // Divide the total range (0 to N^3) among threads
            uint64_t chunk_size = N_power_3 / num_threads;

            // Third nearest neighbor term (3NN)
            for (unsigned t = 0; t < num_threads; ++t) {
                uint64_t start_idx = t * chunk_size;
                uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

                threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_3NN, this, start_idx, end_idx, std::ref(partial_vals[t]), alpha, beta));
            }
            for (auto& thread : threads) thread.join();
            for (const auto& partial_val : partial_vals) val += partial_val;

            // Clear threads after usage
            threads.clear();
        }
    }

    val *= 0.5;
    return val;
}


/// Function to compute indices from a flattened index for E_elas
void V_BaTiO3_parallel::compute_indices_elas(uint64_t idx, int& i, int& j, int& k) {
    k = idx % N;
    idx /= N;
    j = idx % N;
    idx /= N;
    i = idx % N;
}



void V_BaTiO3_parallel::compute_partial_sum_elas_I1(uint64_t start_idx, uint64_t end_idx, double& partial_val,
                                                    const std::shared_ptr<double[]>& v0,
                                                    const std::shared_ptr<double[]>& v1)
{
    partial_val = 0.0;
    int i_plus_1, i_minus_1, j_plus_1, j_minus_1;
    int ind_Ba;
    int ijk_Ba;
    int iPlus1jk_Ba;
    int iMinus1jk_Ba;
    int ijPlus1k_Ba;
    int ijMinus1k_Ba;
    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {
        int i, j, k;
        compute_indices_elas(idx, i, j, k);
        i_plus_1 = python_mod(i + 1, N);
        i_minus_1 = python_mod(i - 1, N);
        j_plus_1 = python_mod(j + 1, N);
        j_minus_1 = python_mod(j - 1, N);
         ind_Ba = 0;

         ijk_Ba = flattened_ind_for_E_elas(i, j, k, ind_Ba);

         iPlus1jk_Ba = flattened_ind_for_E_elas(i_plus_1, j, k, ind_Ba);
        iMinus1jk_Ba = flattened_ind_for_E_elas(i_minus_1, j, k, ind_Ba);

         ijPlus1k_Ba = flattened_ind_for_E_elas(i, j_plus_1, k, ind_Ba);
         ijMinus1k_Ba = flattened_ind_for_E_elas(i, j_minus_1, k, ind_Ba);
        // row 1 and row 2
        partial_val += gamma11_val * (std::pow(v0[ijk_Ba] - v0[iPlus1jk_Ba], 2) +
                                      std::pow(v0[ijk_Ba] - v0[iMinus1jk_Ba], 2));

        // row 3 and row 4
        partial_val += gamma12_val * ((v0[ijk_Ba] - v0[iPlus1jk_Ba]) * (v1[ijk_Ba] - v1[ijPlus1k_Ba]) +
                                      (v0[ijk_Ba] - v0[iMinus1jk_Ba]) * (v1[ijk_Ba] - v1[ijMinus1k_Ba]));

        // row 5 and row 6
        partial_val += gamma44_val * (std::pow(v0[ijk_Ba] - v0[ijPlus1k_Ba] + v1[ijk_Ba] - v1[iPlus1jk_Ba], 2) +
                                      std::pow(v0[ijk_Ba] - v0[ijMinus1k_Ba] + v1[ijk_Ba] - v1[iMinus1jk_Ba], 2));
    }//end idx

}



void V_BaTiO3_parallel::compute_partial_sum_elas_I2(uint64_t start_idx, uint64_t end_idx, double& partial_val,
                                                    const std::shared_ptr<double[]>& v0,
                                                    const std::shared_ptr<double[]>& v2)
{
    partial_val = 0.0;
    int i, j, k;
    int i_plus_1;
    int i_minus_1 ;
    int k_plus_1;
    int k_minus_1;
    int ind_Ba;
    int ijk_Ba;
    int ij_kPlus1;
    int ij_kMinus1;
    int iPlus1_jk ;
    int iMinus1_jk;
    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {

        compute_indices_elas(idx, i, j, k);
        i_plus_1 = python_mod(i + 1, N);
        i_minus_1 = python_mod(i - 1, N);

       k_plus_1 = python_mod(k + 1, N);
        k_minus_1 = python_mod(k - 1, N);
         ind_Ba = 0;

        ijk_Ba = flattened_ind_for_E_elas(i, j, k, ind_Ba);

         ij_kPlus1 = flattened_ind_for_E_elas(i, j, k_plus_1, ind_Ba);
         ij_kMinus1 = flattened_ind_for_E_elas(i, j, k_minus_1, ind_Ba);

        iPlus1_jk = flattened_ind_for_E_elas(i_plus_1, j, k, ind_Ba);
         iMinus1_jk = flattened_ind_for_E_elas(i_minus_1, j, k, ind_Ba);
        // row 1 and row 2
        partial_val += gamma11_val * (std::pow(v2[ijk_Ba] - v2[ij_kPlus1], 2) +
                                      std::pow(v2[ijk_Ba] - v2[ij_kMinus1], 2));

        // row 3 and row 4
        partial_val += gamma12_val * ((v2[ijk_Ba] - v2[ij_kPlus1]) * (v0[ijk_Ba] - v0[iPlus1_jk]) +
                                      (v2[ijk_Ba] - v2[ij_kMinus1]) * (v0[ijk_Ba] - v0[iMinus1_jk]));

        // row 5 and row 6
        partial_val += gamma44_val * (std::pow(v2[ijk_Ba] - v2[iPlus1_jk] + v0[ijk_Ba] - v0[ij_kPlus1], 2) +
                                      std::pow(v2[ijk_Ba] - v2[iMinus1_jk] + v0[ijk_Ba] - v0[ij_kMinus1], 2));
    }//end idx

}



/// Function to calculate the inhomogeneous part E_elas_I3
void V_BaTiO3_parallel::compute_partial_sum_elas_I3(uint64_t start_idx, uint64_t end_idx, double& partial_val,
                                                    const std::shared_ptr<double[]>& v1,
                                                    const std::shared_ptr<double[]>& v2)
{
    partial_val = 0.0;
    int j_plus_1, j_minus_1, k_plus_1, k_minus_1;
    int ind_Ba ;
    int ijk_Ba;
    int i_jPlus1_k ;
    int i_jMinus1_k;
    int ij_kPlus1;
    int ij_kMinus1;
    for (uint64_t idx = start_idx; idx < end_idx; ++idx)
    {
        int i, j, k;
        compute_indices_elas(idx, i, j, k);
        j_plus_1 = python_mod(j + 1, N);
        j_minus_1 = python_mod(j - 1, N);

        k_plus_1 = python_mod(k + 1, N);
        k_minus_1 = python_mod(k - 1, N);
        ind_Ba = 0;

        ijk_Ba = flattened_ind_for_E_elas(i, j, k, ind_Ba);

        i_jPlus1_k = flattened_ind_for_E_elas(i, j_plus_1, k, ind_Ba);
        i_jMinus1_k = flattened_ind_for_E_elas(i, j_minus_1, k, ind_Ba);

        ij_kPlus1 = flattened_ind_for_E_elas(i, j, k_plus_1, ind_Ba);
        ij_kMinus1 = flattened_ind_for_E_elas(i, j, k_minus_1, ind_Ba);
        // row 1 and row 2
        partial_val += gamma11_val * (std::pow(v1[ijk_Ba] - v1[i_jPlus1_k], 2) +
                                      std::pow(v1[ijk_Ba] - v1[i_jMinus1_k], 2));

        // row 3 and row 4
        partial_val += gamma12_val * ((v1[ijk_Ba] - v1[i_jPlus1_k]) * (v2[ijk_Ba] - v2[ij_kPlus1]) +
                                      (v1[ijk_Ba] - v1[i_jMinus1_k]) * (v2[ijk_Ba] - v2[ij_kMinus1]));

        // row 5 and row 6
        partial_val += gamma44_val * (std::pow(v1[ijk_Ba] - v1[ij_kPlus1] + v2[ijk_Ba] - v2[i_jPlus1_k], 2) +
                                      std::pow(v1[ijk_Ba] - v1[ij_kMinus1] + v2[ijk_Ba] - v2[i_jMinus1_k], 2));
    }//end idx
}




int V_BaTiO3_parallel::flattened_ind_for_E_elas(const int& i, const int& j, const int& k, const int& q)
{
    //i,j,k take values from 0 to N-1, q takes values from 0 to 4
    int ind = q + 5 * k + 5 * j * N + 5 * i * N * N;

    return ind;
}

/// Main elastic energy function using std::thread
double V_BaTiO3_parallel::E_elas(const std::shared_ptr<double[]>& eta_H, const std::shared_ptr<double[]>& v0,
                                 const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2) {
    // Homogeneous part
    double eta_H1 = eta_H[0], eta_H2 = eta_H[1], eta_H3 = eta_H[2];
    double eta_H4 = eta_H[3], eta_H5 = eta_H[4], eta_H6 = eta_H[5];
    double N_double = static_cast<double>(N);

    double homogeneous_part = 0.5 * std::pow(N_double, 3) * B11_val * (std::pow(eta_H1, 2) + std::pow(eta_H2, 2) + std::pow(eta_H3, 2)) +
                              std::pow(N_double, 3) * B12_val * (eta_H1 * eta_H2 + eta_H2 * eta_H3 + eta_H3 * eta_H1) +
                              0.5 * std::pow(N_double, 3) * B44_val * (std::pow(eta_H4, 2) + std::pow(eta_H5, 2) + std::pow(eta_H6, 2));

    // Inhomogeneous parts
    double E_elas_I1 = 0.0, E_elas_I2 = 0.0, E_elas_I3 = 0.0;
    uint64_t N_power_3 = N * N * N;

    // Number of threads to use
    const unsigned num_threads = std::thread::hardware_concurrency();
    std::vector<double> partial_vals(num_threads, 0.0);

    // Parallelize for E_elas_I1
    threads.clear();
    uint64_t chunk_size = N_power_3 / num_threads;

    for (unsigned t = 0; t < num_threads; ++t) {
        uint64_t start_idx = t * chunk_size;
        uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

        threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_elas_I1, this, start_idx, end_idx,
                                      std::ref(partial_vals[t]), v0, v1));
    }
    for (auto& thread : threads) thread.join();
    for (const auto& partial_val : partial_vals) E_elas_I1 += partial_val;

    // Parallelize for E_elas_I2
    threads.clear();
    std::fill(partial_vals.begin(), partial_vals.end(), 0.0);

    for (unsigned t = 0; t < num_threads; ++t) {
        uint64_t start_idx = t * chunk_size;
        uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

        threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_elas_I2, this, start_idx, end_idx,
                                      std::ref(partial_vals[t]), v0, v2));
    }
    for (auto& thread : threads) thread.join();
    for (const auto& partial_val : partial_vals) E_elas_I2 += partial_val;

    // Parallelize for E_elas_I3
    threads.clear();
    std::fill(partial_vals.begin(), partial_vals.end(), 0.0);

    for (unsigned t = 0; t < num_threads; ++t) {
        uint64_t start_idx = t * chunk_size;
        uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

        threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_elas_I3, this, start_idx, end_idx,
                                      std::ref(partial_vals[t]), v1, v2));
    }
    for (auto& thread : threads) thread.join();
    for (const auto& partial_val : partial_vals) E_elas_I3 += partial_val;

    return homogeneous_part + E_elas_I1 + E_elas_I2 + E_elas_I3;
}



void V_BaTiO3_parallel::compute_indices_mode_int(uint64_t idx, int& i, int& j, int& k) {
    k = idx % N;
    idx /= N;
    j = idx % N;
    idx /= N;
    i = idx % N;
}


/// Function to calculate the inhomogeneous part of E_elas_mode_int
void V_BaTiO3_parallel::compute_partial_sum_elas_mode_int(uint64_t start_idx, uint64_t end_idx, double& partial_val,
                                                          const std::shared_ptr<double[]>& eta_H,
                                                          const std::shared_ptr<double[]>& v0,
                                                          const std::shared_ptr<double[]>& v1,
                                                          const std::shared_ptr<double[]>& v2) {
    partial_val = 0.0;
    int i, j, k;
    int i_minus_1;
    int j_minus_1;
    int k_minus_1;
    int ind_Ba;
    int ijk;
    int iMinus1_jk;
    int iMinus1_jMinus1_k ;
    int i_jMinus1_k;
    int iMinus1_j_kMinus1;
    int ij_kMinus1;
    int iMinus1_jMinus1_kMinus1;
    int i_jMinus1_kMinus1;
    double Delta_vxx_ijk;
    double etaI1ijk;
    double Delta_vyy_ijk;
    double etaI2ijk;
    double Delta_vzz_ijk;
    double etaI3ijk;
    double Delta_vxy_ijk;
    double Delta_vyx_ijk;
    double etaI6ijk;
    double Delta_vzx_ijk;
    double Delta_vxz_ijk;
    double etaI5ijk;
    double Delta_vyz_ijk;
    double Delta_vzy_ijk;
    double eta1ijk;
    double eta2ijk ;
    double eta3ijk ;
    double eta4ijk;
    double eta5ijk;
    double eta6ijk;
    double eta_H1 = eta_H[0];
    double eta_H2 = eta_H[1];
    double eta_H3 = eta_H[2];

    double eta_H4 = eta_H[3];
    double eta_H5 = eta_H[4];
    double eta_H6 = eta_H[5];
    int u_ind_ijk;
    double E_int_ijk;
    for (uint64_t idx = start_idx; idx < end_idx; ++idx) {

        compute_indices_mode_int(idx, i, j, k);

        i_minus_1 = python_mod(i - 1, N);
         j_minus_1 = python_mod(j - 1, N);
         k_minus_1 = python_mod(k - 1, N);
         ind_Ba = 0;

        // i,j,k
         ijk = flattened_ind_for_E_elas(i, j, k, ind_Ba);
        iMinus1_jk = flattened_ind_for_E_elas_mode_int(i_minus_1, j, k, ind_Ba);
         iMinus1_jMinus1_k = flattened_ind_for_E_elas(i_minus_1, j_minus_1, k, ind_Ba);
         i_jMinus1_k = flattened_ind_for_E_elas(i, j_minus_1, k, ind_Ba);
         iMinus1_j_kMinus1 = flattened_ind_for_E_elas(i_minus_1, j, k_minus_1, ind_Ba);
         ij_kMinus1 = flattened_ind_for_E_elas(i, j, k_minus_1, ind_Ba);
         iMinus1_jMinus1_kMinus1 = flattened_ind_for_E_elas(i_minus_1, j_minus_1, k_minus_1, ind_Ba);
         i_jMinus1_kMinus1 = flattened_ind_for_E_elas(i, j_minus_1, k_minus_1, ind_Ba);

        // Calculate the Delta_vxx, Delta_vyy, Delta_vzz, etc.
         Delta_vxx_ijk = v0[iMinus1_jk] - v0[ijk] + v0[iMinus1_jMinus1_k] - v0[i_jMinus1_k] +
                               v0[iMinus1_j_kMinus1] - v0[ij_kMinus1] + v0[iMinus1_jMinus1_kMinus1] - v0[i_jMinus1_kMinus1];
         etaI1ijk = Delta_vxx_ijk / 4.0;

         Delta_vyy_ijk = v1[i_jMinus1_k] - v1[ijk] + v1[i_jMinus1_kMinus1] - v1[ij_kMinus1] +
                               v1[iMinus1_jMinus1_k] - v1[iMinus1_jk] + v1[iMinus1_jMinus1_kMinus1] - v1[iMinus1_j_kMinus1];
         etaI2ijk = Delta_vyy_ijk / 4.0;

         Delta_vzz_ijk = v2[ij_kMinus1] - v2[ijk] + v2[iMinus1_j_kMinus1] - v2[iMinus1_jk] +
                               v2[i_jMinus1_kMinus1] - v2[i_jMinus1_k] + v2[iMinus1_jMinus1_kMinus1] - v2[iMinus1_jMinus1_k];
        etaI3ijk = Delta_vzz_ijk / 4.0;

       Delta_vxy_ijk = v1[iMinus1_jk] - v1[ijk] + v1[iMinus1_jMinus1_k] - v1[i_jMinus1_k] +
                               v1[iMinus1_j_kMinus1] - v1[ij_kMinus1] + v1[iMinus1_jMinus1_kMinus1] - v1[iMinus1_jMinus1_k];
         Delta_vyx_ijk = v0[i_jMinus1_k] - v0[ijk] + v0[iMinus1_jMinus1_k] - v0[iMinus1_jk] +
                               v0[i_jMinus1_kMinus1] - v0[ij_kMinus1] + v0[iMinus1_jMinus1_kMinus1] - v0[iMinus1_jMinus1_k];
       etaI6ijk = (Delta_vxy_ijk + Delta_vyx_ijk) / 4.0;

        Delta_vzx_ijk = v0[ij_kMinus1] - v0[ijk] + v0[iMinus1_j_kMinus1] - v0[iMinus1_jk] +
                               v0[i_jMinus1_kMinus1] - v0[i_jMinus1_k] + v0[iMinus1_jMinus1_kMinus1] - v0[iMinus1_j_kMinus1];
         Delta_vxz_ijk = v2[iMinus1_jk] - v2[ijk] + v2[iMinus1_j_kMinus1] - v2[ij_kMinus1] +
                               v2[iMinus1_jMinus1_k] - v2[i_jMinus1_k] + v2[iMinus1_jMinus1_kMinus1] - v2[iMinus1_j_kMinus1];
         etaI5ijk = (Delta_vzx_ijk + Delta_vxz_ijk) / 4.0;

        Delta_vyz_ijk = v2[i_jMinus1_k] - v2[ijk] + v2[i_jMinus1_kMinus1] - v2[ij_kMinus1] +
                               v2[iMinus1_jMinus1_k] - v2[iMinus1_jk] + v2[iMinus1_jMinus1_kMinus1] - v2[i_jMinus1_kMinus1];
         Delta_vzy_ijk = v1[ij_kMinus1] - v1[ijk] + v1[i_jMinus1_kMinus1] - v1[i_jMinus1_k] +
                               v1[iMinus1_j_kMinus1] - v1[iMinus1_jk] + v1[iMinus1_jMinus1_kMinus1] - v1[i_jMinus1_kMinus1];
        double etaI4ijk = (Delta_vyz_ijk + Delta_vzy_ijk) / 4.0;

        // Calculate eta values
         eta1ijk = eta_H1 + etaI1ijk;
         eta2ijk = eta_H2 + etaI2ijk;
        eta3ijk = eta_H3+ etaI3ijk;
         eta4ijk = eta_H4 + etaI4ijk;
        eta5ijk = eta_H5 + etaI5ijk;
         eta6ijk = eta_H6 + etaI6ijk;

        // Calculate energy contribution at this point
         u_ind_ijk = i * N_power_2 + j * N + k;
       E_int_ijk = 0.0;

        E_int_ijk += B100_val * eta1ijk * u0[u_ind_ijk] * u0[u_ind_ijk];
        E_int_ijk += B111_val * eta1ijk * u1[u_ind_ijk] * u1[u_ind_ijk];
        E_int_ijk += B122_val * eta1ijk * u2[u_ind_ijk] * u2[u_ind_ijk];
        E_int_ijk += B200_val * eta2ijk * u0[u_ind_ijk] * u0[u_ind_ijk];
        E_int_ijk += B211_val * eta2ijk * u1[u_ind_ijk] * u1[u_ind_ijk];
        E_int_ijk += B222_val * eta2ijk * u2[u_ind_ijk] * u2[u_ind_ijk];
        E_int_ijk += B300_val * eta3ijk * u0[u_ind_ijk] * u0[u_ind_ijk];
        E_int_ijk += B311_val * eta3ijk * u1[u_ind_ijk] * u1[u_ind_ijk];
        E_int_ijk += B322_val * eta3ijk * u2[u_ind_ijk] * u2[u_ind_ijk];
        E_int_ijk += B412_val * eta4ijk * u1[u_ind_ijk] * u2[u_ind_ijk];
        E_int_ijk += B421_val * eta4ijk * u2[u_ind_ijk] * u1[u_ind_ijk];
        E_int_ijk += B502_val * eta5ijk * u0[u_ind_ijk] * u2[u_ind_ijk];
        E_int_ijk += B520_val * eta5ijk * u2[u_ind_ijk] * u0[u_ind_ijk];
        E_int_ijk += B601_val * eta6ijk * u0[u_ind_ijk] * u1[u_ind_ijk];
        E_int_ijk += B610_val * eta6ijk * u1[u_ind_ijk] * u0[u_ind_ijk];

        partial_val += 0.5 * E_int_ijk;
    }
}



/// Main elastic energy with mode interaction function using std::thread
double V_BaTiO3_parallel::E_elas_mode_int(const std::shared_ptr<double[]>& eta_H, const std::shared_ptr<double[]>& v0,
                                         const std::shared_ptr<double[]>& v1, const std::shared_ptr<double[]>& v2) {
    double val = 0.0;

    // Number of threads to use
    const unsigned num_threads = std::thread::hardware_concurrency();
    std::vector<double> partial_vals(num_threads, 0.0);

    // Parallelize the computation for E_elas_mode_int
    threads.clear();
    uint64_t chunk_size = N_power_3 / num_threads;

    for (unsigned t = 0; t < num_threads; ++t) {
        uint64_t start_idx = t * chunk_size;
        uint64_t end_idx = (t == num_threads - 1) ? N_power_3 : (t + 1) * chunk_size;

        threads.push_back(std::thread(&V_BaTiO3_parallel::compute_partial_sum_elas_mode_int, this, start_idx, end_idx,
                                      std::ref(partial_vals[t]), eta_H, v0, v1, v2));
    }

    for (auto& thread : threads) thread.join();
    for (const auto& partial_val : partial_vals) val += partial_val;

    return val;
}

int V_BaTiO3_parallel::flattened_ind_for_E_elas_mode_int(const int& i, const int& j, const int& k, const int& q)
{
    //i,j,k take values from 0 to N1, q takes values from 0 to 4
    int ind = q + 5 * k + 5 * j * N + 5 * i * N * N;

    return ind;
}