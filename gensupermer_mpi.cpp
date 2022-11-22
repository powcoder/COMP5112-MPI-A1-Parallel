https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
// *********************************************************************
//     [NAME]:    [STUDENT ID]: 
//     [EMAIL]: 
//     NOTICE: Write your code only in the specified section.
// *********************************************************************
// 7 MAR (update2.1), 28 FEB (update1): UPDATES IN read2supermers(...)
#define _in_
#define _out_
#define _MPI_TEST_
// #define DEBUG
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include "utilities.hpp"
#ifdef _MPI_TEST_
#include "mpi.h"
#endif
using namespace std;

// void read2supermers(const char* _read, int read_len, int k, int p, _out_ char* &supermers, _out_ /**/int* &supermer_offs, int &n_supermers);
void read2supermers(const char* _read, int read_len, int k, int p, _out_ vector<string> &supermers);

const int MAX_PROCESS = 64;
int K, P;

int main(int argc, char **argvs) {
    #ifdef _MPI_TEST_
    MPI_Init(&argc, &argvs);
    MPI_Comm comm;
    int num_process; // number of processors
    int my_rank;     // my global rank
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &num_process);
    MPI_Comm_rank(comm, &my_rank);
    #endif

    int correctness_checking = 0;
    string output_path;
    string read_file_path;
    ArgParser(argc, argvs, K, P, read_file_path, correctness_checking, output_path);
    vector<string> reads;
    
    // Input data (the reads in CSR format)
    int num_of_reads = 0;
    char* reads_CSR;
    /**/int* reads_CSR_offs;

    // Output data, each supermers should be a string in the vector
    // you need to save all the supermers to the vector below in the root(0) process
    vector<string> all_supermers;

    #ifdef _MPI_TEST_
    if (my_rank == 0) {
        cout<<"MPI running with "<<num_process<<" threads."<<endl<<endl;
    #endif
        LoadReadsFromFile(read_file_path.c_str(), reads);
        Vector2CSR(reads, num_of_reads, reads_CSR, reads_CSR_offs);
        cout << reads.size() << " reads loaded from "<< read_file_path << endl << endl;
    #ifdef _MPI_TEST_
    }
    #endif

    // time measurement starts
    auto start_time = chrono::high_resolution_clock::now();
    
    #ifdef _MPI_TEST_
    // hint: Input: "num_of_reads", "reads_CSR", "reads_CSR_offs", "K", "P"
    //       You need to save all the generated supermers in the vector "all_supermers" in Process 0.
    // you need to do:
    //       1. Scatter the read data to each MPI processes.
    //       2. Perform the super-mer generation in each process. 
    //          (you can refer to the sequential version to know the usage of the function read2supermers(...))
    //       3. Gather all the super-mers to the root process and store in the vector "all_supermers". (The order in the vector doesn't matter.)
    
    // ==============================================================
    // ==============================================================
    // ====    Write your implementation only below this line    ====
    // ==============================================================
    
    
    // ==============================================================
    // ====    Write your implementation only above this line    ====
    // ==============================================================
    // ==============================================================
    #endif
    
    #ifdef _MPI_TEST_
    if (my_rank == 0) {
    #endif
        // time measurement ends
        auto end_time = chrono::high_resolution_clock::now();
        auto duration_sec = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count()/1000.0;
        cout << "Your algorithm finished in " << duration_sec << " sec." << endl << endl;
        
        // output to text file and correctness checking
        delete reads_CSR;
        delete reads_CSR_offs;
        if (correctness_checking) CorrectnessChecking(reads, K, P, all_supermers);
        if (!output_path.empty()) {
            if (!correctness_checking) sort(all_supermers.begin(), all_supermers.end());
            SaveSupermers(output_path, all_supermers);
        }
    #ifdef _MPI_TEST_
    }
    MPI_Barrier(comm);
    // cout<<"Thread "<<my_rank<<" ends."<<endl;
    
    MPI_Finalize();
    #endif
    
    return 0;
}

/*
This function receives a C-style read string, the length of the read,
k (length of k-mer), p (length of miniizer), and output the supermers 
which can be generated from this read to a vector<string>.
*/
void read2supermers(const char* _read, int read_len, int k, int p, _out_ vector<string> &supermers) {
    string prev_minimizer, minimizer, new_minimizer;
    string read(_read, read_len); // from-buffer init
    int i, j;
    char base;
    int skm_begin_pos, skm_end_pos, mm_begin_pos;
    
    // Generate the first k-mer's minimizer:
    skm_begin_pos = 0;
    skm_end_pos = k;
    mm_begin_pos = 0;
    minimizer = new_minimizer = read.substr(0, p);
    for (i=p; i<k; i++) {
        new_minimizer = new_minimizer.substr(1, p-1) + read[i]; // UPDATE1
        if (new_minimizer <= minimizer) minimizer = new_minimizer, mm_begin_pos = i-p+1;
    }

    // Continue generating minimizers:
    for (i=1; i<read_len-k+1; i++) { // i: the beginning position of the current k-mer
        if (i > mm_begin_pos) {
            // new minimizer required
            prev_minimizer = minimizer;
            minimizer = new_minimizer = read.substr(i, p);
            for (j=i+p; j<i+k; j++) {
                new_minimizer = new_minimizer.substr(1, p-1) + read[j]; // UPDATE1
                if (new_minimizer <= minimizer) minimizer = new_minimizer, mm_begin_pos = j-p+1;
            }
            // if the new minimizer equals to the previous one, we can continue
            if (minimizer != prev_minimizer) {
                skm_end_pos = i-1+k;
                supermers.push_back(read.substr(skm_begin_pos, skm_end_pos-skm_begin_pos)); // save the supermer
                skm_begin_pos = i;
            }
        }
        else {
            new_minimizer = read.substr(i+k-p, p); // UPDATE1
            if (new_minimizer < minimizer) { // save the supermer
                skm_end_pos = i-1+k;
                supermers.push_back(read.substr(skm_begin_pos, skm_end_pos-skm_begin_pos));
                skm_begin_pos = i;
                minimizer = new_minimizer, mm_begin_pos = i+k-1-p+1;
            }
            if (new_minimizer == minimizer) mm_begin_pos = i+k-1-p+1; // UPDATE1
        }
    } // UPDATE 2.1
    skm_end_pos = read_len;
    supermers.push_back(read.substr(skm_begin_pos, skm_end_pos-skm_begin_pos));
}
