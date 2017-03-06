
#include "CUDAGlobalOptimization.h"
#include "CPUGlobalOptimization.h"
#include "interval.h"
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <exception>


void testFunc(double *inBox, int inRank, int inNumBoxesSplitCoeffCPU, int inNumBoxesSplitCoeffGPU, double inEps, int inMaxIter, int inFun, double * outBox, double * outMin, double *outEps, int *status, std::ofstream &fout){

	time_t t_startCPU, t_endCPU,t_startGPU, t_endGPU;
	try{
		switch(inFun){
		case 1:
			fout << "Произведение параметров (" << inRank << " переменных)" << "\t";
			time(&t_startCPU);
			fnGetOptValueOnCPU(inBox, inRank, inNumBoxesSplitCoeffCPU, inEps, inMaxIter, fnCalcFunLimitsMultiple2, outBox,outMin, outEps, status);
			time(&t_endCPU);
			
			time(&t_startGPU);
			fnGetOptValueWithCUDA(inBox, inRank, inNumBoxesSplitCoeffGPU, inEps, inMaxIter, 1, outBox,outMin, outEps, status);
			time(&t_endGPU);
			break;
		case 2:
			fout << "Гиперболическая функция (" << inRank << " переменных)" << "\t";
			time(&t_startCPU);
			fnGetOptValueOnCPU(inBox, inRank, inNumBoxesSplitCoeffCPU, inEps, inMaxIter, fnCalcFunLimitsHypebolic2, outBox,outMin, outEps, status);
			time(&t_endCPU);
			
			time(&t_startGPU);
			fnGetOptValueWithCUDA(inBox, inRank, inNumBoxesSplitCoeffGPU, inEps, inMaxIter, 2, outBox,outMin, outEps, status);
			time(&t_endGPU);
			break;
		case 3:
			fout << "Функция Алуффи-Пентини (" << inRank << " переменных)" << "\t";
			time(&t_startCPU);
			fnGetOptValueOnCPU(inBox, inRank, inNumBoxesSplitCoeffCPU, inEps, inMaxIter, fnCalcFunLimitsAluffiPentini2, outBox,outMin, outEps, status);
			time(&t_endCPU);
			
			time(&t_startGPU);
			fnGetOptValueWithCUDA(inBox, inRank, inNumBoxesSplitCoeffGPU, inEps, inMaxIter, 3, outBox,outMin, outEps, status);
			time(&t_endGPU);
			break;
		case 4:
			fout << "Функция розенброка (" << inRank << " переменных)" << "\t";
			time(&t_startCPU);
			fnGetOptValueOnCPU(inBox, inRank, inNumBoxesSplitCoeffCPU, inEps, inMaxIter, fnCalcFunLimitsRozenbroke, outBox,outMin, outEps, status);
			time(&t_endCPU);
			
			time(&t_startGPU);
			fnGetOptValueWithCUDA(inBox, inRank, inNumBoxesSplitCoeffGPU, inEps, inMaxIter, 4, outBox,outMin, outEps, status);
			time(&t_endGPU);
			break;
		}
	}
	catch(std::exception ex){
		std::cout << "Произошла ошибка при решении задачи оптимизации: " << ex.what() << "\n";
		exit(EXIT_FAILURE); 
	}
	fout << inEps << "\t";
	fout << inNumBoxesSplitCoeffCPU << "\t";
	fout << inNumBoxesSplitCoeffGPU << "\t";
	fout << (t_endCPU-t_startCPU) << "\t";
	fout << (t_endGPU-t_startGPU) << "\n";
}

int main()
{
    int inRank = 10;
	time_t t_start, t_end;

	std::ofstream fout;
	fout.open("test_function.txt",std::ios::app);

	

	double *inBox = new double[inRank*2];
	double *outBox = new double[inRank*2];
	double outMin = 0.0;
	double inEps = 0.0001;
	double outEps = 0.0;
	double inMaxIter = 100;
	int inNumBoxesSplitCoeff = 2;
	int status = -1;
	int numTestCycles = 10;


	inRank = 3;
	for(int i = 0; i < inRank; i++)
	{
		inBox[i*2] = -20.0;
		inBox[i*2+1] = 20.0;
	}

	time(&t_start);
	//fnGetOptValueOnCPU(inBox, inRank, 2, inEps, inMaxIter, fnCalcFunLimitsMultiple2, outBox,&outMin, &outEps, &status);
	fnGetOptValueWithCUDA(inBox, inRank, 8, inEps, inMaxIter, 4, outBox,&outMin, &outEps, &status);
	time(&t_end);
			
	//time(&t_start);
	//fnGetOptValueWithCUDA(inBox, inRank, inNumBoxesSplitCoeffGPU, inEps, inMaxIter, 1, outBox,outMin, outEps, status);
	//time(&t_end);

	//testFunc(inBox,2,2,16,inEps,inMaxIter,1,outBox,&outMin,&outEps,&status,fout);
	//std::cout << "Test #1 passed!\n";

	//testFunc(inBox,2,2,16,inEps,inMaxIter,2,outBox,&outMin,&outEps,&status,fout);
	//std::cout << "Test #2 passed!\n";

	//testFunc(inBox,2,2,16,inEps,inMaxIter,3,outBox,&outMin,&outEps,&status,fout);
	//std::cout << "Test #3 passed!\n";

	//testFunc(inBox,2,2,16,inEps,inMaxIter,4,outBox,&outMin,&outEps,&status,fout);
	//std::cout << "Test #4 passed!\n";

	//testFunc(inBox,3,2,4,inEps,inMaxIter,4,outBox,&outMin,&outEps,&status,fout);
	//std::cout << "Test #5 passed!\n";

	//testFunc(inBox,4,2,4,inEps,inMaxIter,4,outBox,&outMin,&outEps,&status,fout);
	//std::cout << "Test #6 passed!\n";

	//testFunc(inBox,5,2,2,inEps,inMaxIter,4,outBox,&outMin,&outEps,&status,fout);
	//std::cout << "Test #7 passed!\n";

	//testFunc(inBox,6,2,2,inEps,inMaxIter,4,outBox,&outMin,&outEps,&status,fout);
	//std::cout << "Test #8 passed!\n";

	//testFunc(inBox,7,2,2,inEps,inMaxIter,4,outBox,&outMin,&outEps,&status,fout);
	//std::cout << "Test #9 passed!\n";



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

	std::cin.get();

	fout.close();

	delete [] inBox;
	delete [] outBox;

	return 0;
}




