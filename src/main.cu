#include <iostream>
#include <cstddef> //size_t
#include <algorithm> //swap
#include <map>
#include <fstream>
#include "kernels.cuh"
#include "functions.cuh"
#include "forces.cuh"
#include "DEFINITIONS.cuh" //BLOCKS, THREADS, REAL

int main() {

    std::ifstream inFile;

    std::map<std::string, REAL> constants;
    std::string key;
    REAL value;

    /*********** R E A D  I N S H A P E *****************/
    std::string line;
    size_t N = 0;

    inFile.open("../shape_cpp.txt");
    while (std::getline(inFile, line) && line!="")
        ++N;
    std::cout<<"Number of nodes (N) = "<<N<<std::endl;
    std::cout<<"|----------"<<std::endl;
    inFile.close();


    auto *var_cpu = new REAL[4*N];
    for (int i = 0; i<4*N; ++i)
        var_cpu[i] = 0.0;

    REAL x,y;
    //std::ifstream inFile; // number of filled lines must be N
    inFile.open("../shape_cpp.txt");
    for (int i=0; i<N; ++i) {
        inFile >> x >> y,
                get_x(i, var_cpu, N) = x;
        get_y(i, var_cpu, N) = y;
        //std::cout<<x<<" "<<y<<std::endl;
    }

    inFile.close();
    /****************************************************/

    /************** R E A D  C O N S T A N T S *******************/
    inFile.open("../constants.txt");

    while (not inFile.eof()) {
        inFile >> key >> value;
        constants.emplace(key, value);
        //std::cout<<"key = "<<key<<"; value = "<<value<<std::endl;
    }

    inFile.close();

    std::cout<<"CONSTANTS :"<<std::endl;

    REAL mass = constants.at("mass"); std::cout << "|  mass = " << mass << std::endl;
    REAL betta = constants.at("betta"); std::cout << "|  betta = " << betta << std::endl;
    REAL l0 = constants.at("l0"); std::cout << "|  l0 = " << l0 << std::endl;
    REAL kl = constants.at("kl"); std::cout << "|  kl = " << kl << std::endl;
    REAL s0 = constants.at("s0"); std::cout << "|  s0 = " << s0 << std::endl;
    REAL ks = constants.at("ks"); std::cout << "|  ks = " << ks << std::endl;
    REAL kb = constants.at("kb"); std::cout << "|  kb = " << kb << std::endl;
    REAL lr = constants.at("lr"); std::cout << "|  lr = " << lr << std::endl;
    REAL kr = constants.at("kr"); std::cout << "|  kr = " << kr << std::endl;
    REAL lp = constants.at("lp"); std::cout << "|  lp = " << lp << std::endl;
    REAL kp = constants.at("kp"); std::cout << "|  kp = " << kp << std::endl;

    std::cout<<"|----------"<<std::endl;
    /**************************************************************/



    /*
    get_x(0, var_cpu, N) = 0.0;
    get_y(0, var_cpu, N) = 0.0;
    get_x(1, var_cpu, N) = 1.0;
    get_y(1, var_cpu, N) = 0.0;
    get_x(2, var_cpu, N) = 1.1;
    get_y(2, var_cpu, N) = 1.1;
    get_x(3, var_cpu, N) = 0.0;
    get_y(3, var_cpu, N) = 1.0;
    */


    /********** M E M O R Y  A L L O C A T I O N ********/
    REAL *var = 0, *var_temp = 0;

    cudaMalloc( (void**)&var, 4*N*sizeof(REAL) );
    cudaMemcpy( var, var_cpu, 4*N*sizeof(REAL),
                cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&var_temp, 4*N*sizeof(REAL) );

    //REAL *temp;
    //cudaMalloc( (void**)&temp, sizeof(REAL) );


    REAL *area = 0;

    cudaMalloc( (void**)&area, sizeof(REAL) );


    Lock lock; // MUTEX for reduce kernels


    REAL *fl_x = 0, *fl_y = 0;
    REAL *fs_x = 0, *fs_y = 0;
    REAL *fb_x = 0, *fb_y = 0;

    cudaMalloc( (void**)&fl_x, N*sizeof(REAL) );
    cudaMalloc( (void**)&fl_y, N*sizeof(REAL) );
    cudaMalloc( (void**)&fs_x, N*sizeof(REAL) );
    cudaMalloc( (void**)&fs_y, N*sizeof(REAL) );
    cudaMalloc( (void**)&fb_x, N*sizeof(REAL) );
    cudaMalloc( (void**)&fb_y, N*sizeof(REAL) );

    REAL *f_x = 0, *f_y = 0;

    cudaMalloc( (void**)&f_x, N*sizeof(REAL) );
    cudaMalloc( (void**)&f_y, N*sizeof(REAL) );
    /****************************************************/


    /****************** R E A D  T I M E S **************/
    inFile.open("../times.txt");

    while (not inFile.eof()) {
        inFile >> key >> value;
        constants.emplace(key, value);
        //std::cout<<"key = "<<key<<"; value = "<<value<<std::endl;
    }

    inFile.close();

    std::cout<<"TIMES :"<<std::endl;

    REAL t_step = constants.at("t_step"); std::cout << "|  t_step = " << t_step << std::endl;
    REAL t_end = constants.at("t_end"); std::cout << "|  t_end = " << t_end << std::endl;
    REAL t_print = constants.at("t_print"); std::cout << "|  t_print = " << t_print << std::endl;

    std::cout<<"|----------"<<std::endl;
    /****************************************************/


    /// First Frame
    std::ofstream outFile;
    outFile.open("../out.txt");

    for (int i=0; i<2*N; ++i)
        outFile << var_cpu[i] << " ";
    outFile << "\n";

    REAL t = 0.0, t_p = 0.0;
    while (t <= t_end) {

        std::cout<<std::endl<<"step: "<<t<<" / "<<t_end<<std::endl;

        /****** F O R C E  C A L C U L A T I O N ********************************/
        calculate_fl<<<BLOCKS,THREADS>>>(var, N, fl_x, fl_y, l0, kl);

        fill<<<BLOCKS,THREADS>>>(area, 0.0, 1);
        calculate_area<<<BLOCKS,THREADS>>>(lock, area, var, N);
        calculate_fs<<<BLOCKS,THREADS>>>(var, area, N, fs_x, fs_y, s0, ks);

        calculate_fb<<<BLOCKS,THREADS>>>(var, N, fb_x, fb_y, kb);

        fill<<<BLOCKS,THREADS>>>(f_x, 0.0, N);
        fill<<<BLOCKS,THREADS>>>(f_y, 0.0, N);

        add<<<BLOCKS,THREADS>>>(f_x, f_x, fl_x, N); //f_x = f_x + fl_x
        add<<<BLOCKS,THREADS>>>(f_x, f_x, fs_x, N);
        add<<<BLOCKS,THREADS>>>(f_x, f_x, fb_x, N);

        add<<<BLOCKS,THREADS>>>(f_y, f_y, fl_y, N);
        add<<<BLOCKS,THREADS>>>(f_y, f_y, fs_y, N);
        add<<<BLOCKS,THREADS>>>(f_y, f_y, fb_y, N);
        /************************************************************************/


        /**************** S T E P P E R *****************************************/
        rhs<<<BLOCKS,THREADS>>>(var_temp, var, N, f_x, f_y, mass, betta);

        // WHOLE STEP - E U L E R !!!
        //TODO : change to do_step(...)
        mul<<<BLOCKS,THREADS>>>(var_temp, var_temp, t_step, 4*N);
        add<<<BLOCKS,THREADS>>>(var_temp, var_temp, var, 4*N);

        std::swap(var, var_temp);
        /************************************************************************/

        t += t_step;
        t_p += t_step;

        /*************** M A K E  S N A P S H O T *******************************/
        if (t_p >= t_print) {

            cudaMemcpy(var_cpu, var, 4 * N * sizeof(REAL),
                       cudaMemcpyDeviceToHost);

            for (int i=0; i<2*N; ++i)
                outFile << var_cpu[i] << " ";
            outFile << "\n";

            t_p = 0.0;
        }
        /************************************************************************/

    }

    outFile.close();


    /************************* M E M O R Y  F R E E I N G ***********************/
    cudaFree(var);
    cudaFree(var_temp);

    cudaFree(area);
    //cudaFree(temp);

    cudaFree(fl_x);
    cudaFree(fl_y);

    cudaFree(fs_x);
    cudaFree(fs_y);

    cudaFree(fb_x);
    cudaFree(fb_y);

    cudaFree(f_x);
    cudaFree(f_y);

    delete [] var_cpu;
    /*****************************************************************************/

    //std::cout<<std::endl;
    return 0;

}