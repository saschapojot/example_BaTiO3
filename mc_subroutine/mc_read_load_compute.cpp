//
// Created by polya on 9/12/24.
//

#include "mc_read_load_compute.hpp"

void mc_computation::load_pickle_data(const std::string& filename, std::shared_ptr<double[]> data_ptr, std::size_t size)
{

    // Initialize Python and NumPy
    Py_Initialize();
    np::initialize();

    try {
        // Use Python's 'io' module to open the file directly in binary mode
        py::object io_module = py::import("io");
        py::object file = io_module.attr("open")(filename, "rb");  // Open file in binary mode

        // Import the 'pickle' module
        py::object pickle_module = py::import("pickle");

        // Use pickle.load to deserialize from the Python file object
        py::object loaded_data = pickle_module.attr("load")(file);

        // Close the file
        file.attr("close")();

        // Check if the loaded object is a NumPy array
        if (py::extract<np::ndarray>(loaded_data).check()) {
            np::ndarray np_array = py::extract<np::ndarray>(loaded_data);

            // Convert the NumPy array to a Python list using tolist()
            py::object py_list = np_array.attr("tolist")();

            // Ensure the list size matches the expected size
            ssize_t list_size = py::len(py_list);
            if (static_cast<std::size_t>(list_size) > size) {
                throw std::runtime_error("The provided shared_ptr array size is smaller than the list size.");
            }

            // Copy the data from the Python list to the shared_ptr array
            for (ssize_t i = 0; i < list_size; ++i) {
                data_ptr[i] = py::extract<double>(py_list[i]);
            }
        } else {
            throw std::runtime_error("Loaded data is not a NumPy array.");
        }
    }
    catch (py::error_already_set&) {
        PyErr_Print();
        throw std::runtime_error("Python error occurred.");
    }

}


/// load data by flushNum
/// @param flushNum
    void mc_computation::load_data(const int& flushNum)
    {
        std::string name;
    if (flushNum==-1)
    {
        name="init.pkl";
    }else
    {
        name="flushEnd"+std::to_string(flushNum)+".pkl";
    }

    //load v0

    std::string v0FileName=U_dist_dataDir+"/v0/v0_"+name;

    load_pickle_data(v0FileName,v0_init,elemNumTot_v);

    // load v1
    std::string v1FileName=U_dist_dataDir+"/v1/v1_"+name;
    load_pickle_data(v1FileName,v1_init,elemNumTot_v);

    //load v2
    std::string v2FileName=U_dist_dataDir+"/v2/v2_"+name;
    load_pickle_data(v2FileName,v2_init,elemNumTot_v);

    //load eta_H
    std::string eta_HFileName=U_dist_dataDir+"/eta_H/eta_H_"+name;
    load_pickle_data(eta_HFileName,eta_H_init,elemNumTot_v);



    }



double mc_computation::generate_nearby_normal(const double & x,const double &sigma)
{

    std::normal_distribution<> dist_normal(x,sigma);

    double xNext=dist_normal(e2);

    return xNext;

}

void mc_computation::proposal(const std::shared_ptr<double[]> & vecCurr,std::shared_ptr<double[]>&vecNext,const int&pos, const int &vecLength)
{
    double elem_next=generate_nearby_normal(vecCurr[pos],h);

    std::memcpy(vecNext.get(),vecCurr.get(),vecLength*sizeof(double));
    vecNext[pos]=elem_next;
}



double mc_computation::acceptanceRatio(const double &xCurr, const double&UCurr, const double& xNext, const double& UNext)
{
    double numerator = -this->beta*UNext;

    double denominator=-this->beta*UCurr;

    double ratio = std::exp(numerator - denominator);

    return std::min(1.0, ratio);

}


void mc_computation::execute_mc_one_sweep(std::shared_ptr<double[]>& v0_Curr,std::shared_ptr<double[]>&v1_Curr,std::shared_ptr<double[]>&v2_Curr,
        std::shared_ptr<double[]>& eta_H_Curr,double &UCurr,
        std::shared_ptr<double[]>&v0_Next,std::shared_ptr<double[]>& v1_Next, std::shared_ptr<double[]>&v2_Next,std::shared_ptr<double[]>&eta_H_Next
        ,const int &fls, const int& swp)
{

    //next U
    double UNext;

    //first update eta_H

    for(int j=0;j<6;j++)
    {

        int pos_eta_H=randint_0_5(e2);
        double eta_H_elem_Curr=eta_H_Curr[pos_eta_H];
        double eta_H_elem_Next;

    }




}
