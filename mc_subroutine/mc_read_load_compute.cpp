//
// Created by polya on 9/12/24.
//

#include "mc_read_load_compute.hpp"

void mc_computation::load_pickle_data(const std::string& filename, std::shared_ptr<double[]>& data_ptr, std::size_t size)
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




void mc_computation::initialize_v0_v1_v2_eta_H()
{
    std::string name;
    std::string eta_H_inFileName,v0_inFileName, v1_inFileName, v2_inFileName;
    if(this->flushLastFile==-1)
    {
       name= "init";
        eta_H_inFileName=out_eta_H_path+"/eta_H_"+name+".pkl";
        v0_inFileName=out_v0_path+"/v0_"+name+".pkl";
        v1_inFileName=out_v1_path+"/v1_"+name+".pkl";
        v2_inFileName=out_v2_path+"/v2_"+name+".pkl";
    }
    else
    {
        name="flushEnd"+std::to_string(this->flushLastFile);
        eta_H_inFileName=out_eta_H_path+"/"+name+".eta_H.pkl";
        // std::cout<<"eta_H_inFileName="<<eta_H_inFileName<<std::endl;
        v0_inFileName=out_v0_path+"/"+name+".v0.pkl";
        v1_inFileName=out_v1_path+"/"+name+".v1.pkl";
        v2_inFileName=out_v2_path+"/"+name+".v2.pkl";
    }

    //load eta_H
    this->load_pickle_data(eta_H_inFileName,eta_H_data_ptr,6*sweepToWrite);
    std::memcpy(eta_H_init.get(),eta_H_data_ptr.get()+6*(sweepToWrite-1),6*sizeof(double));

    //load v0
    this->load_pickle_data(v0_inFileName,v0_data_ptr,elemNumTot_v*sweepToWrite);
    std::memcpy(v0_init.get(),v0_data_ptr.get()+elemNumTot_v*(sweepToWrite-1),elemNumTot_v*sizeof(double));

    //load v1
    this->load_pickle_data(v1_inFileName,v1_data_ptr,elemNumTot_v*sweepToWrite);
    std::memcpy(v1_init.get(),v1_data_ptr.get()+elemNumTot_v*(sweepToWrite-1),elemNumTot_v*sizeof(double));

    //load v2
    this->load_pickle_data(v2_inFileName,v2_data_ptr,elemNumTot_v*sweepToWrite);
    std::memcpy(v2_init.get(),v2_data_ptr.get()+elemNumTot_v*(sweepToWrite-1),elemNumTot_v*sizeof(double));





}

double mc_computation::generate_nearby_normal(const double & x,const double &sigma)
{

    std::normal_distribution<> dist_normal(x,sigma);

    double xNext=dist_normal(e2);

    return xNext;

}

void mc_computation::proposal(const std::shared_ptr<double[]> & vecCurr,std::shared_ptr<double[]>&vecNext,const int&pos, const int &vecLength,const double& sigma)
{
    double elem_next=generate_nearby_normal(vecCurr[pos],sigma);

    std::memcpy(vecNext.get(),vecCurr.get(),vecLength*sizeof(double));
    vecNext[pos]=elem_next;
}



double mc_computation::acceptanceRatio( const double&UCurr, const double& UNext)
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

    for(int j=0;j<elemNumTot_eta_H;j++)
    {

        int pos_eta_H=randint_0_5(e2);
        this->proposal(eta_H_Curr,eta_H_Next,pos_eta_H,elemNumTot_eta_H,h_eta_H);

        UCurr=(*potFuncPtr)(eta_H_Curr,v0_Curr,v1_Curr,v2_Curr);
        UNext=(*potFuncPtr)(eta_H_Next,v0_Curr,v1_Curr,v2_Curr);

        double r=this->acceptanceRatio(UCurr,UNext);
        double u = distUnif01(e2);
        if(u<=r){
            UCurr=UNext;
            std::memcpy(eta_H_Curr.get(),eta_H_Next.get(),elemNumTot_eta_H*sizeof(double));
        }




    }//end updating eta_H

    //update v0
    for(int l=0;l<elemNumTot_v;l++){
        int i=randint_0_N_minus1(e2);
        int j=randint_0_N_minus1(e2);
        int k=randint_0_N_minus1(e2);
        int q=randint_0_5(e2);
        int ind= flattened_ind_for_v(i,j,k,q);

        this->proposal(v0_Curr,v0_Next,ind,elemNumTot_v,h_v);
        UCurr=(*potFuncPtr)(eta_H_Curr,v0_Curr,v1_Curr,v2_Curr);
        UNext=(*potFuncPtr)(eta_H_Curr,v0_Next,v1_Curr,v2_Curr);
        double r=this->acceptanceRatio(UCurr,UNext);
        double u = distUnif01(e2);
        if(u<=r){
            UCurr=UNext;
            std::memcpy(v0_Curr.get(),v0_Next.get(),elemNumTot_v*sizeof(double ));
        }



    }//end updating v0

    //update v1
    for(int l=0;l<elemNumTot_v;l++){
        int i=randint_0_N_minus1(e2);
        int j=randint_0_N_minus1(e2);
        int k=randint_0_N_minus1(e2);
        int q=randint_0_5(e2);
        int ind= flattened_ind_for_v(i,j,k,q);
        this->proposal(v1_Curr,v1_Next,ind,elemNumTot_v,h_v);
        UCurr=(*potFuncPtr)(eta_H_Curr,v0_Curr,v1_Curr,v2_Curr);
        UNext=(*potFuncPtr)(eta_H_Curr,v0_Curr,v1_Next,v2_Curr);
        double r=this->acceptanceRatio(UCurr,UNext);
        double u = distUnif01(e2);
        if(u<=r) {
            UCurr = UNext;
            std::memcpy(v1_Curr.get(),v1_Next.get(),elemNumTot_v*sizeof(double ));
        }

    }//end updating v1


//update v2
for(int l=0;l<elemNumTot_v;l++){
    int i=randint_0_N_minus1(e2);
    int j=randint_0_N_minus1(e2);
    int k=randint_0_N_minus1(e2);
    int q=randint_0_5(e2);
    int ind= flattened_ind_for_v(i,j,k,q);
    this->proposal(v2_Curr,v2_Next,ind,elemNumTot_v,h_v);
    UCurr=(*potFuncPtr)(eta_H_Curr,v0_Curr,v1_Curr,v2_Curr);
    UNext=(*potFuncPtr)(eta_H_Curr,v0_Curr,v1_Curr,v2_Next);
    double r=this->acceptanceRatio(UCurr,UNext);
    double u = distUnif01(e2);
    if(u<=r) {
        UCurr = UNext;
        std::memcpy(v2_Curr.get(),v2_Next.get(),elemNumTot_v*sizeof(double ));
    }

}//end updating v1

}


int mc_computation::flattened_ind_for_v(const int& i, const int& j, const int& k, const int& q){
    //i,j,k take values from 0 to N-1, q takes values from 0 to 4
    int ind = q + 5 * k + 5 * j * N + 5 * i * N * N;

    return ind;

}


void mc_computation::execute_mc(const std::shared_ptr<double[]>& v0Vec,const std::shared_ptr<double[]>& v1Vec, const std::shared_ptr<double[]>& v2Vec,const std::shared_ptr<double[]>& eta_HVec, const int & flushNum){

    std::shared_ptr<double[]> v0_Curr=std::shared_ptr<double[]>(new double[elemNumTot_v], std::default_delete<double[]>());
    std::shared_ptr<double[]> v0_Next=std::shared_ptr<double[]>(new double[elemNumTot_v], std::default_delete<double[]>());
    std::memcpy(v0_Curr.get(),v0Vec.get(),elemNumTot_v*sizeof(double ));

    std::shared_ptr<double[]> v1_Curr=std::shared_ptr<double[]>(new double[elemNumTot_v], std::default_delete<double[]>());
    std::shared_ptr<double[]> v1_Next=std::shared_ptr<double[]>(new double[elemNumTot_v], std::default_delete<double[]>());
    std::memcpy(v1_Curr.get(),v1Vec.get(),elemNumTot_v*sizeof(double ));

    std::shared_ptr<double[]> v2_Curr=std::shared_ptr<double[]>(new double[elemNumTot_v], std::default_delete<double[]>());
    std::shared_ptr<double[]> v2_Next=std::shared_ptr<double[]>(new double[elemNumTot_v], std::default_delete<double[]>());
    std::memcpy(v2_Curr.get(),v2Vec.get(),elemNumTot_v*sizeof(double ));

    std::shared_ptr<double[]> eta_H_Curr=std::shared_ptr<double[]>(new double[elemNumTot_eta_H], std::default_delete<double[]>());
    std::shared_ptr<double[]> eta_H_Next=std::shared_ptr<double[]>(new double[elemNumTot_eta_H], std::default_delete<double[]>());
    std::memcpy(eta_H_Curr.get(),eta_HVec.get(),elemNumTot_eta_H*sizeof(double ));

    // int sweepStart = sweepInit;
    double UCurr=0;
    int flushThisFileStart=this->flushLastFile+1;
    int sweepStart =flushThisFileStart*sweepToWrite*sweep_multiple;

    for (int fls = 0; fls < flushNum; fls++) {
        const auto tMCStart{std::chrono::steady_clock::now()};

        for (int swp = 0; swp < sweepToWrite*sweep_multiple; swp++) {
            // std::cout<<"swp="<<swp<<std::endl;
            execute_mc_one_sweep(v0_Curr,v1_Curr,v2_Curr,eta_H_Curr,UCurr,v0_Next,v1_Next,v2_Next,eta_H_Next,fls,swp);
            if(swp%sweep_multiple==0)
            {
                int swp_out=swp/sweep_multiple;
                U_data_ptr[swp_out]=UCurr;
                std::memcpy(v0_data_ptr.get()+swp_out*elemNumTot_v,v0_Curr.get(),elemNumTot_v*sizeof(double));
                std::memcpy(v1_data_ptr.get()+swp_out*elemNumTot_v,v1_Curr.get(),elemNumTot_v*sizeof(double));
                std::memcpy(v2_data_ptr.get()+swp_out*elemNumTot_v,v2_Curr.get(),elemNumTot_v*sizeof(double));
                std::memcpy(eta_H_data_ptr.get()+swp_out*elemNumTot_eta_H,eta_H_Curr.get(),elemNumTot_eta_H*sizeof(double ));
            }//end swp mod

        }//end sweep for
        int sweepEnd = sweepStart + sweepToWrite*sweep_multiple - 1;
        int flushEnd=flushThisFileStart+fls;
        std::string fileNameMiddle =  "flushEnd" + std::to_string(flushEnd);

        std::string out_U_PickleFileName = out_U_path+"/" + fileNameMiddle + ".U.pkl";

        std::string out_v0_PickleFileName=out_v0_path+"/"+fileNameMiddle+".v0.pkl";
        std::string out_v1_PickleFileName=out_v1_path+"/"+fileNameMiddle+".v1.pkl";
        std::string out_v2_PickleFileName=out_v2_path+"/"+fileNameMiddle+".v2.pkl";

        std::string out_eta_H_PickleFileName=out_eta_H_path+"/"+fileNameMiddle+".eta_H.pkl";

        save_array_to_pickle(U_data_ptr.get(),sweepToWrite,out_U_PickleFileName);

        save_array_to_pickle(v0_data_ptr.get(),sweepToWrite * elemNumTot_v,out_v0_PickleFileName);
        save_array_to_pickle(v1_data_ptr.get(),sweepToWrite * elemNumTot_v,out_v1_PickleFileName);
        save_array_to_pickle(v2_data_ptr.get(),sweepToWrite * elemNumTot_v,out_v2_PickleFileName);

        save_array_to_pickle(eta_H_data_ptr.get(),sweepToWrite * elemNumTot_eta_H,out_eta_H_PickleFileName);
        const auto tMCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
        std::cout << "flush " + std::to_string(flushEnd)  + ": "
                  << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;

    }//end flush for loop
    std::cout << "mc executed for " << flushNum << " flushes." << std::endl;

}
// void mc_computation::save_array_to_pickle(double *ptr,const int& size,const std::string& filename){
//     using namespace boost::python;
//     try {
//         Py_Initialize();  // Initialize the Python interpreter
//         if (!Py_IsInitialized()) {
//             throw std::runtime_error("Failed to initialize Python interpreter");
//         }
//
//         // Debug output
//         //        std::cout << "Python interpreter initialized successfully." << std::endl;
//
//         // Import the pickle module
//         object pickle = import("pickle");
//         object pickle_dumps = pickle.attr("dumps");
//
//         // Create a Python list from the C++ array
//         list py_list;
//         for (std::size_t i = 0; i < size; i++) {
//             py_list.append(ptr[i]);
//         }
//
//         // Serialize the list using pickle.dumps
//         object serialized_array = pickle_dumps(py_list);
//
//         // Extract the serialized data as a string
//         std::string serialized_str = extract<std::string>(serialized_array);
//
//         // Write the serialized data to a file
//         std::ofstream file(filename, std::ios::binary);
//         if (!file) {
//             throw std::runtime_error("Failed to open file for writing");
//         }
//         file.write(serialized_str.data(), serialized_str.size());
//         file.close();
//
//         // Debug output
//         //        std::cout << "Array serialized and written to file successfully." << std::endl;
//     } catch (const error_already_set&) {
//         PyErr_Print();
//         std::cerr << "Boost.Python error occurred." << std::endl;
//     } catch (const std::exception& e) {
//         std::cerr << "Exception: " << e.what() << std::endl;
//     }
//
//     if (Py_IsInitialized()) {
//         Py_Finalize();  // Finalize the Python interpreter
//     }
//
//
//
// }

void mc_computation::save_array_to_pickle(double *ptr, const int& size, const std::string& filename) {
    using namespace boost::python;
    namespace np = boost::python::numpy;

    // Initialize Python interpreter if not already initialized
    if (!Py_IsInitialized()) {
        Py_Initialize();
        if (!Py_IsInitialized()) {
            throw std::runtime_error("Failed to initialize Python interpreter");
        }
        np::initialize(); // Initialize NumPy
    }

    try {
        // Import the pickle module
        object pickle = import("pickle");
        object pickle_dumps = pickle.attr("dumps");

        // Convert C++ array to NumPy array
        np::ndarray numpy_array = np::from_data(
            ptr,                                  // Pointer to the C++ array data
            np::dtype::get_builtin<double>(),      // NumPy data type (double)
            boost::python::make_tuple(size),       // Shape of the array (1D array)
            boost::python::make_tuple(sizeof(double)),  // Strides
            object()                               // Optional base object
        );

        // Serialize the NumPy array using pickle.dumps
        object serialized_array = pickle_dumps(numpy_array);

        // Extract the serialized data as a string
        std::string serialized_str = extract<std::string>(serialized_array);

        // Write the serialized data to a file
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to open file for writing");
        }
        file.write(serialized_str.data(), serialized_str.size());
        file.close();

        // Debug output (optional)
        // std::cout << "Array serialized and written to file successfully." << std::endl;

    } catch (const error_already_set&) {
        PyErr_Print();
        std::cerr << "Boost.Python error occurred." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    // Finalize the Python interpreter at the end of the program, or manage separately
    // Py_Finalize();
}
void mc_computation::init_and_run()
{
    this->initialize_v0_v1_v2_eta_H();

    this->execute_mc(v0_init,v1_init,v2_init,eta_H_init,newFlushNum);


}
