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