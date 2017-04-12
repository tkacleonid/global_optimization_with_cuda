#ifndef INTERVAL_H
#define INTERVAL_H

#include <iostream>
#include <stdio.h>
#include <exception>
#include <string>
#include <vector>
#include <algorithm>

#define DEVICE 0
#define MAX_NUM_RUNS (100)

struct Interval {
    double bottom;
    double up;   
};

struct Box {
    Interval *i;
    Interval funBounds;
    double funValue;
    
    Box(int inRank) {
        i = new Interval[inRank];
    }
    ~Box(){
        delete [] i;
    }
};


int const BLOCK_SIZE = 128;

#define CHECKED_CALL(func)                                     \
    do {                                                       \
        cudaError_t err = (func);                              \
        if (err != cudaSuccess) {                              \
            printf("%s(%d): ERROR: %s returned %s (err#%d)\n", \
                   __FILE__, __LINE__, #func,                  \
                   cudaGetErrorString(err), err);              \
            exit(EXIT_FAILURE);                                \
        }                                                      \
    } while (0)

#endif
