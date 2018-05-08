#include <iostream>
#include <cstddef> //size_t
#include <algorithm> //swap
#include <fstream>
#include "kernels.cuh"
#include "functions.cuh"
#include "forces.cuh"
#include "DEFINITIONS.cuh" //BLOCKS, THREADS, REAL

int main() {

    size_t N = 113;
    auto *var_cpu = new REAL[4*N];
    for (int i = 0; i<4*N; ++i)
        var_cpu[i] = 0.0;


    /*********** R E A D  I N S H A P E *****************/
    std::ifstream inFile; // number of lines must be N
    inFile.open("shape_cpp.txt");
    REAL x,y;

    for (int i=0; i<N; ++i) {
        inFile >> x >> y,
        get_x(i, var_cpu, N) = x;
        get_y(i, var_cpu, N) = y;
        //std::cout<<x<<" "<<y<<std::endl;
    }

    inFile.close();
    /****************************************************/


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


    REAL t = 0.0, t_step = 0.01, t_end = 100.0;
    REAL t_print = t_end / 100.0, t_p = 0.0;


    /// First Frame
    std::ofstream outFile;
    outFile.open("out.txt");

    for (int i=0; i<2*N; ++i)
        outFile << var_cpu[i] << " ";
    outFile << "\n";


    while (t <= t_end) {

        //std::cout<<std::endl<<"step: "<<t<<" / "<<t_end<<std::endl;

        /****** F O R C E  C A L C U L A T I O N ********************************/
        calculate_fl<<<BLOCKS,THREADS>>>(var, N, fl_x, fl_y, 0.66, 0.3);

        fill<<<BLOCKS,THREADS>>>(area, 0.0, 1);
        calculate_area<<<BLOCKS,THREADS>>>(lock, area, var, N);
        calculate_fs<<<BLOCKS,THREADS>>>(var, area, N, fs_x, fs_y, 314.0, 2000.0);

        calculate_fb<<<BLOCKS,THREADS>>>(var, N, fb_x, fb_y, 0.05);

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
        rhs<<<BLOCKS,THREADS>>>(var_temp, var, N, f_x, f_y, 1.0, 0.5);

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