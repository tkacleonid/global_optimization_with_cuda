
#include "CPUGlobalOptimization.h"
#include "interval.h"
#include <time.h>
#include <stdio.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <exception>

int main() {
    int inRank = 6;

    double *inBox = new double[inRank * 2];
    double *outBox = new double[inRank * 2];
    double outMin = 0.0;
    double inEps = 0.001;
    double outEps = 0.0;
    double inMaxIter = 200000000;
    int inNumBoxesSplitCoeff = 2;
    int status = -1;


    for (int i = 0; i < inRank; i++) {
        inBox[i * 2] = -2000.0;
        inBox[i * 2 + 1] = 2000.0;
    }

    /*
            t_start =clock();
            //fnGetOptValueOnCPU(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsRozenbroke, outBox,&outMin, &outEps, &status);
            fnGetOptValueOnCPU(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsAluffiPentini2, outBox,&outMin, &outEps, &status);
            t_end =clock();


            std::cout << "Result: ";
            for(int i = 0; i < inRank; i++)
            {
                    std::cout << "[" << outBox[i*2] << "; " << outBox[i*2 + 1]  << "]\t";
            }

            std::cout << "\n";
            std::cout << "min = " << outMin << "\t";
            std::cout << "eps = " << outEps;
            std::cout << "\n";
            std::cout << "time = " << (t_end-t_start) << "\t";
            std::cout << "\n";
     */
    auto start = std::chrono::high_resolution_clock::now();


    fnGetOptValueOnCPUSortDFS(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsRozenbroke, outBox,&outMin, &outEps, &status);
    //fnGetOptValueOnCPUSortDFS(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsAluffiPentini2, outBox, &outMin, &outEps, &status);
    
    //fnGetOptValueOnCPUSort(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsRozenbroke, outBox,&outMin, &outEps, &status);
    //fnGetOptValueOnCPUSort(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsAluffiPentini2, outBox, &outMin, &outEps, &status);
    
    
    
    auto end = std::chrono::high_resolution_clock::now();


    std::cout << "Result: ";
    for (int i = 0; i < inRank; i++) {
        std::cout << "[" << outBox[i * 2] << "; " << outBox[i * 2 + 1] << "]\t";
    }

    std::cout << "\n";
    std::cout << "min = " << outMin << "\t";
    std::cout << "eps = " << outEps;
    std::cout << "\n";
    std::cout << "time in millisecs: " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - start)).count() << "\t";
    std::cout << "\n";

    delete [] inBox;
    delete [] outBox;

    return 0;
}




