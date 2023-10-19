#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>
#include <chrono>
#include <fstream> //to export vectors

#include "strassen_mm.hpp"

// to compile use g++ -std=c++14 -o xstrass strassen_mm.cpp strass_test.cpp

int trials = 3; //establish the number of trials
int matrix_size = 128;  //final size of matrix
int ii = 0; //index for flops vector
int step = 2;  //step size for NN


// Set up a random number generator.  Normal distribution between 1 and 0
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dist(0.0, 1.0); // range [0,1]

//code to perform the FLOP calculation of daxpy
int main() {


       // initialize the vector for flops 
       std::vector<double> vec_flops(matrix_size/step -1);  // has one less element that size of matrix
       
        for (int NN=2; NN <= matrix_size; NN += step) {

            int m = NN;  //number of rows
            int p = NN;  //number of columns
            int n = NN;

            std::cout << "NN " << NN << std::endl;
  

            // initialize the matrices
            std::vector<std::vector<double>> A(NN, std::vector<double> (NN));
            std::vector<std::vector<double>> B(NN, std::vector<double> (NN));
            std::vector<std::vector<double>> C(NN, std::vector<double> (NN));

            // place random values into the matrices
            for (int i = 0; i < NN; i++) {
                for (int j = 0; j < NN; j++) {
                    A[i][j] = dist(gen);
                    B[i][j] = dist(gen);
                }
            }

            double a = dist(gen); //generate scalar 'a'
            double b = dist(gen); //generate scalar 'a'


            // create FLOP cout, (7/4)m^3+(9/2)m^2
            long double flop_count = static_cast<long double>((7/4)*m*m*m+(9/2)*m*m);
            long double nano_to_secondL = 1.e-9L;

            // create parameter sum duration in for loop trials
            long double time_trials = 0;
              

              //loop to perform multiple trials of the function
              for (int ii=0; ii < trials; ii++) {


                 

                     // start time
                     auto start = std::chrono::high_resolution_clock::now();
                     
                     // perform algebraic operation by calling function
                     C = strassen_mm (A, B);

                     
                     // stop time and compute duration
                     auto stop = std::chrono::high_resolution_clock::now();
                     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds> (stop-start);
                     long double elapsed_time = duration.count()*1e-9; //convert to seconds
                     time_trials += elapsed_time; //capture the sum of time for all the trials

                    // std::cout << "duration in loop = " << elapsed_time << std::endl;
                    // std::cout << "total flops " << flop_count << std::endl;

              }

              // calculate Mega flops per second
              long double flops = flop_count*trials;
              long double mflops = flops/time_trials/1e6;
            
              vec_flops[ii] = mflops; // place value into a vector
              ii += 1; //update index for flops vector


        // check to see if matrix C is being generated
/*
        std::cout << "C matrix is " << std::endl;
        for (int i = 0; i < NN; i++) {
            for (int j = 0; j < NN; j++) {
            // print the value at (i, j) position
                std::cout << C[i][j] << " ";
            }
        std::cout << std::endl;     
        }
*/
       //       std::cout <<"final y =  " << final_y << std::endl;
       //       std::cout << "Mflop count = " << flops/1e6 << std::endl;
       //       std::cout << "ii = " << ii << std::endl;
       //       std::cout << "vec_flops[ii] " << vec_flops[ii] << std::endl;
       //       std::cout << "MFLOPS = " << mflops << std::endl;
       //       std::cout << "size of vector x " << x.size() << std::endl;
        }

       // std::cout << "size of flops vector " << vec_flops.size() << std::endl;
       
       //code to export vec_flops
    std::ofstream outfile("strassen_flops.csv");
    if (outfile.is_open()) {
        for (int i = 0; i < vec_flops.size(); i++) {
            outfile << vec_flops[i] << std::endl;
        }
        outfile.close();
        std::cout << "Vector has been exported to 'strassen_flops.csv'" << std::endl;
    } else {
        std::cout << "Unable to open file" << std::endl;
    }

       // print out the flops vector
       std::cout << "flops vector " << std::endl;
       for (int i =0; i< ii ; i++ ) {
            std::cout << vec_flops[i] << std::endl;

       }


        return 0;
}