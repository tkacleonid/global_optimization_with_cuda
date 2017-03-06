#ifndef __CUDAGLOBALOPTIMIZATION_H__
#define __CUDAGLOBALOPTIMIZATION_H__

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cuda_runtime_api.h>
#include <cuda_fp16.h>
#include <cuda_device_runtime_api.h>
#include <math_functions_dbl_ptx3.h>
#include "interval.h"
#include <stdio.h>
#include <stdlib.h>

// Send Data to GPU to calculate limits
void sendDataToCuda(double *outLimits, const double *inBox, int inRank, int inFunc, int numBoxes);

/**
*	Calculus Interval for Multiple function on GPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
__device__ void fnCalcFunLimitsMultiple2_CUDA(const double *inBox, int inRank, double *outLimits);
/**
*	Calculus Interval for Hyperbolic function on GPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
__device__ void fnCalcFunLimitsHypebolic2_CUDA(const double *inBox, int inRank, double *outLimits);
/**
*	Calculus Interval for AluffiPentini function
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
__device__ void fnCalcFunLimitsAluffiPentini2_CUDA(const double *inBox, int inRank, double *outLimits);
/**
*	Calculus Valuefor Multiple function on GPU
*	@param inbox pointer to point
*	@param inRank number of variables
*	@return value of function
*/
__device__ double fnCalcFunRozenbroke_CUDA(double *inBox, int inRank);
/**
*	Calculus Interval for Rozenbroke function
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
__device__ void fnCalcFunLimitsRozenbroke_CUDA(const double *inBox, int inRank, double *outLimits);
/**
*	Calculus minimum value for function on GPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param inNumBoxesSplitCoeff number of parts for each dimension
*	@param inEps required accuracy
*	@param inMaxIter maximum count of iterations
*	@param inFun number of optimazing function
*	@param outBox pointer to optimal box
*	@param outMin pointer to optimal value
*	@param outEps pointer to reached accuracy
*	@param outEps pointer to status of solving optimization problem
*/
void fnGetOptValueWithCUDA(double *inBox, int inRank, int inNumBoxesSplitCoeff, double inEps, int inMaxIter, int inFunc, double *outBox, double*outMin, double *outEps,int *status);
/**
*	Send data into GPU
*	@param outLimits pointer to calculated limits
*	@param inBox pointer to boxes
*	@param inRank demnsion
*	@param inFunc number of optimazion function
*	@param inNumBoxes number of boxes
*/
void sendDataToCuda(double *outLimits, const double *inBox, int inRank, int inFunc, int inNumBoxes);

__global__ void calculateLimitsOnCUDA(double *outLimits, const double *inBox,int inRank,int inFunc);

__device__ double fnCalcFunMultiple2_CUDA(const double *inBox, int inRank)
{
	return inBox[0]*inBox[1];
}

/**
*	Calculus Interval for Multiple function on GPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
__device__ void fnCalcFunLimitsMultiple2_CUDA(const double *inBox, int inRank, double *outLimits)
{
	
	double x1 = (inBox[0]+inBox[1])/2;
	double x2 = (inBox[2]+inBox[3])/2;

	double var1 = inBox[0]*inBox[2];
	double var2 = inBox[0]*inBox[3];
	double var3 = inBox[1]*inBox[2];
	double var4 = inBox[1]*inBox[3];

	outLimits[0] = fmin(fmin(var1,var2),fmin(var3,var4));
	outLimits[1] = fmax(fmax(var1,var2),fmax(var3,var4));
	outLimits[2] = x1*x2;
}


/**
*	Calculus Interval for Hyperbolic function on GPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
__device__ void fnCalcFunLimitsHypebolic2_CUDA(const double *inBox, int inRank, double *outLimits)
{
	double limits[2];
	double limits2[2];

	double x1 = (inBox[0]+inBox[1])/2;
	double x2 = (inBox[2]+inBox[3])/2;

	double var1 = inBox[0]*inBox[0];
	double var2 = inBox[0]*inBox[1];
	double var3 = inBox[1]*inBox[1];

	limits[0] = var2 < 0 ? 0 : fmin(fmin(var1,var2),var3);
	limits[1] = fmax(fmax(var1,var2),var3);

	var1 = inBox[2]*inBox[2];
	var2 = inBox[2]*inBox[3];
	var3 = inBox[3]*inBox[3];

	limits2[0] = var2 < 0 ? 0 : fmin(fmin(var1,var2),var3);
	limits2[1] = fmax(fmax(var1,var2),var3);

	outLimits[0] = limits[0] - limits2[1];
	outLimits[1] = limits[1] - limits2[0];
	outLimits[2] = x1*x1-x2*x2;
}

/**
*	Calculus Interval for AluffiPentini function
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/

__device__ void fnCalcFunLimitsAluffiPentini2_CUDA(const double *inBox, int inRank, double *outLimits)
{

	double limits[2];
	double limits2[2];

	double x1 = (inBox[0]+inBox[1])/2;
	double x2 = (inBox[2]+inBox[3])/2;

	double var1 = inBox[0]*inBox[0];
	double var2 = inBox[0]*inBox[1];
	double var3 = inBox[1]*inBox[1];

	limits[0] = var2 < 0 ? 0 : fmin(fmin(var1,var2),var3);
	limits[1] = fmax(fmax(var1,var2),var3);

	var1 = inBox[2]*inBox[2];
	var2 = inBox[2]*inBox[3];
	var3 = inBox[3]*inBox[3];

	limits2[0] = limits[0]*limits[0];
	limits2[1] = limits[1]*limits[1];

	outLimits[0] = 0.25*limits2[0] - 0.5*limits[1] + 0.1*inBox[0] + 0.5*(var2 < 0 ? 0 : fmin(fmin(var1,var2),var3));
	outLimits[1] = 0.25*limits2[1] - 0.5*limits[0] + 0.1*inBox[1] + 0.5*(fmax(fmax(var1,var2),var3));
	outLimits[2] = 0.25*pow(x1,4.0)-0.5*pow(x1,2.0) + 0.1*x1 + 0.5*pow(x2,2.0);
}


__device__ double fnCalcFunRozenbroke_CUDA(double *inBox, int inRank)
{
	int i;
	double val = 0;;
	for(i = 0; i < inRank - 1; i++)
	{
		val += ((1 - inBox[i])*(1 - inBox[i]) + 100.0*(inBox[i+1]-inBox[i]*inBox[i])*(inBox[i+1]-inBox[i]*inBox[i]));
	}
	return val;
}



/**
*	Calculus Interval for Rozenbroke function
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
__device__ void fnCalcFunLimitsRozenbroke_CUDA(const double *inBox, int inRank, double *outLimits)
{
	double sup = 0;
	double sub = 0;
	double sup1,sub1,sup2,sub2,a,b,val = 0,var1,var2,var3,x1,x2;
	int i,j;

	for(i = 0; i < inRank - 1; i++)
	{
		sub1 = 1 - inBox[i*2 + 1];
		sup1 = 1 - inBox[i*2];

		var1 = sup1*sup1;
		var2 = sup1*sub1;
		var3 = sub1*sub1;
		
		sub1 = (sub1*sup1 < 0) ? 0 : fmin(fmin(var1,var2),var3);
		sup1 = fmax(fmax(var1,var2),var3);

		var1 = inBox[i*2 + 1]*inBox[i*2 + 1];
		var2 = inBox[i*2 + 1]*inBox[i*2];
		var3 = inBox[i*2]*inBox[i*2];

		a = (inBox[i*2 + 1]*inBox[i*2] < 0) ? 0 : fmin(fmin(var1,var2),var3);
		b = fmax(fmax(var1,var2),var3);

		sub2 = inBox[(i+1)*2] - b;
		sup2 = inBox[(i+1)*2 + 1] - a;

		var1 = sup2*sup2;
		var2 = sup2*sub2;
		var3 = sub2*sub2;

		sub2 = (sub2*sup2 < 0) ? 0 : 100*fmin(fmin(var1,var2),var3);
		sup2 = 100*fmax(fmax(var1,var2),var3);

		sub += sub1 + sub2;
		sup += sup1 + sup2;

		x1 = (inBox[i*2 + 1] + inBox[i*2])/2;
		x2 = (inBox[(i+1)*2 + 1] + inBox[(i+1)*2])/2;
		val += ((1 - x1)*(1 - x1) + 100*(x2-x1*x1)*(x2-x1*x1));
	}

	
	//double *x;
	//x = (double *) malloc(inRank*sizeof(double));
	//double minFun;
	//int index;
	//double md;
	//double pw;
	//double a1;
	//double b1;
	//int numFunValues = __double2int_rd(pow(2.0,__int2double_rn(inRank)));

	//for(j = 0; j < inRank; j++){
	//		x[j] = (inBox[j*2]+inBox[j*2+1])/2.0;
	//}
	//minFun = fnCalcFunRozenbroke_CUDA(x, inRank);
	//for(i = 0; i < numFunValues; i++){
	//	for(j = 0; j < inRank; j++){
	//		a1 = __int2double_rn(i);
	//		b1 = __int2double_rn(j+1);
	//		md = fmod(a1,pow(2.0,b1));
	//		pw = pow(2.0,__int2double_rn(b1-1.0));
	//		index = __double2int_rd(md / pw);
	//		x[j] = inBox[j*2+index];
	//	}
	//	val = fnCalcFunRozenbroke_CUDA(x, inRank);
	//	if(minFun > val) minFun = val;
	//}

	/*free(x);
*/
	outLimits[0] = inBox[0];
	outLimits[1] = sup;
	outLimits[2] = inBox[1];
}



void fnGetOptValueWithCUDA(double *inBox, int inRank, int inNumBoxesSplitCoeff, double inEps, int inMaxIter, int inFunc, double *outBox, double*outMin, double *outEps,int *status)
{
	int curNumBoxesSplitCoeff = inNumBoxesSplitCoeff;
	int numBoxes = pow((double) curNumBoxesSplitCoeff,inRank);
	double *boxes =  new double[numBoxes*inRank*2];
	double *boxesResult = new double[numBoxes*3];
	double *restBoxes = new double[inRank*2];
	double *tempRestBoxes = NULL;
	double *h = new double[inRank];
	int numNewBoxes = 0;

	memcpy(restBoxes,inBox,inRank*2*sizeof(double));

	int numRestBoxes = 1;
	int index = 0;
	int i,j,k,n;
	double temp;

	int countIter = 0;
	double curEps = inEps*10;
	
	double funcMin = 0;
	int boxMinIndex = 0;


	*status = 1;
	while((countIter < inMaxIter) && (curEps >= inEps))
	{
		/*curNumBoxesSplitCoeff = (int) (inNumBoxesSplitCoeff/pow(numRestBoxes,1.0/inRank)) + 2;
		numBoxes = pow((double) curNumBoxesSplitCoeff,inRank);
		boxes = new double[numRestBoxes*numBoxes*inRank*2];
		boxesResult = new double[numRestBoxes*numBoxes*3];
		for(k = 0; k < numRestBoxes; k++)
		{
			for(i = 0; i < inRank; i++)
			{
				h[i] = (restBoxes[(k*inRank+i)*2 + 1] - restBoxes[(k*inRank+i)*2])/curNumBoxesSplitCoeff;
			}

			for(n = 0; n < numBoxes; n++)
			{
				for(i = 0; i < inRank; i++)
				{
					index = ((n % numBoxes) % (long) pow((double)curNumBoxesSplitCoeff,i+1))/((long)pow((double)curNumBoxesSplitCoeff,i));
					boxes[((k*numBoxes + n)*inRank+i)*2] = restBoxes[(k*inRank+i)*2] + h[i]*index;
					boxes[((k*numBoxes + n)*inRank+i)*2 + 1] = restBoxes[(k*inRank+i)*2] + h[i]*(index+1);
				}
			}
		}*/

				curNumBoxesSplitCoeff = 8;//(int) (inNumBoxesSplitCoeff/pow(numRestBoxes,1.0/inRank)) + 2;
		numBoxes = pow((double) curNumBoxesSplitCoeff,inRank);
		boxes = new double[numRestBoxes*numBoxes*inRank*2];
		boxesResult = new double[numRestBoxes*numBoxes*3];
		for(k = 0; k < numRestBoxes; k++)
		{
			for(i = 0; i < inRank; i++)
			{
				h[i] = (restBoxes[(k*inRank+i)*2 + 1] - restBoxes[(k*inRank+i)*2])/curNumBoxesSplitCoeff;
			}
			#pragma omp parallel for num_threads(5) private(index,i)
			for(n = 0; n < numBoxes; n++)
			{
				for(i = 0; i < inRank; i++)
				{
					index = ((n % numBoxes) % (long) pow((double)curNumBoxesSplitCoeff,i+1))/((long)pow((double)curNumBoxesSplitCoeff,i));
					boxes[((k*numBoxes + n)*inRank+i)*2] = restBoxes[(k*inRank+i)*2] + h[i]*index;
					boxes[((k*numBoxes + n)*inRank+i)*2 + 1] = restBoxes[(k*inRank+i)*2] + h[i]*(index+1);
					//std::cout << "[" << boxes[((k*numBoxes + n)*inRank+i)*2] << ";" << boxes[((k*numBoxes + n)*inRank+i)*2 + 1] << "]\t";
				}
				//inFun(&boxes[((k*numBoxes + n)*inRank)*2],inRank,&boxesResult[(k*numBoxes + n)*3]);
				//std::cout <<  "\n";
			}
		}




		sendDataToCuda(boxesResult, boxes, inRank, inFunc, numRestBoxes*numBoxes);

		for(n = 0; n < numRestBoxes*numBoxes; n++)
		{
			for(j=0; j < inRank; j++)
			{
				std::cout << "[" << boxes[(n*inRank+j)*2] << ";" << boxes[(n*inRank+j)*2 + 1] << "] ";
			}
			
			std::cout << boxesResult[n*3] << " " << boxesResult[n*3+2] << " " << "\n";
		}
		funcMin = boxesResult[2];
		boxMinIndex = 0;
		for(n = 0; n < numRestBoxes*numBoxes; n++)
		{
			for(i = n + 1; i < numRestBoxes*numBoxes; i++)
			{
				if(boxesResult[n*3] > boxesResult[i*3])
				{
					temp = boxesResult[n*3];
					boxesResult[n*3] = boxesResult[i*3];
					boxesResult[i*3] = temp;

					temp = boxesResult[n*3+1];
					boxesResult[n*3+1] = boxesResult[i*3+1];
					boxesResult[i*3+1] = temp;

					temp = boxesResult[n*3+2];
					boxesResult[n*3+2] = boxesResult[i*3+2];
					boxesResult[i*3+2] = temp;

					for(j=0; j < inRank; j++)
					{
						temp = boxes[(n*inRank+j)*2];
						boxes[(n*inRank+j)*2] = boxes[(i*inRank+j)*2];
						boxes[(i*inRank+j)*2] = temp;

						temp = boxes[(n*inRank+j)*2+1];
						boxes[(n*inRank+j)*2+1] = boxes[(i*inRank+j)*2+1];
						boxes[(i*inRank+j)*2+1] = temp;
					}
				}
				if(funcMin > boxesResult[n*3 + 2] ) {funcMin = boxesResult[n*3+2];boxMinIndex = n;}
			}
			if(funcMin < boxesResult[n*3]) break;
		}

		std::cout << boxesResult[n*3] << "\t" << boxesResult[n*3+2] << "\t" << funcMin << "\n\n";

		curEps = std::abs(boxesResult[0] - funcMin);
		*outEps = curEps;
		*outMin = funcMin;
		memcpy(outBox,boxes + n*inRank*2,inRank*2*sizeof(double));
		if(curEps < inEps)
		{
			*status = 0;
			return;
		}
		numNewBoxes = n;

		tempRestBoxes = new double[numNewBoxes*inRank*2];
		memcpy(tempRestBoxes,boxes,numNewBoxes*inRank*2*sizeof(double));
		if(countIter > 0) delete [] restBoxes;
		restBoxes = tempRestBoxes;

		delete [] boxes;
		delete [] boxesResult;

		numRestBoxes = numNewBoxes;
		countIter++;
	}
	delete [] h;
}


__global__ void calculateLimitsOnCUDA(double *outLimits, const double *inBox,int inRank,int inFunc)
{
    int thread_id = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
	switch(inFunc){
	case 1:
			fnCalcFunLimitsMultiple2_CUDA(&inBox[thread_id*inRank*2], inRank, &outLimits[thread_id*3]);
			break;
	case 2:
			fnCalcFunLimitsHypebolic2_CUDA(&inBox[thread_id*inRank*2], inRank, &outLimits[thread_id*3]);
			break;
	case 3:
			fnCalcFunLimitsAluffiPentini2_CUDA(&inBox[thread_id*inRank*2], inRank, &outLimits[thread_id*3]);
			break;
	case 4:
			fnCalcFunLimitsRozenbroke_CUDA(&inBox[thread_id*inRank*2], inRank, &outLimits[thread_id*3]);
			break;
	}
}


// Send Data to GPU to calculate limits
void sendDataToCuda(double *outLimits, const double *inBox, int inRank, int inFunc, int numBoxes)
{
    double *dev_outLimits = 0;
    double *dev_inBox = 0;
	int GridSize = ((numBoxes%BLOCK_SIZE == 0) ? numBoxes/BLOCK_SIZE : numBoxes/BLOCK_SIZE + 1);
	int numThreads = GridSize*BLOCK_SIZE;
	int sizeOutLimits = numThreads*3*sizeof(double);
	int sizeInBox = numThreads*inRank*2*sizeof(double);

	cudaEvent_t start, stop;

	CHECKED_CALL(cudaSetDevice(DEVICE));
    CHECKED_CALL(cudaMalloc((void **)&dev_outLimits, sizeOutLimits));
    CHECKED_CALL(cudaMalloc((void **)&dev_inBox, sizeInBox));
    //CHECKED_CALL(cudaEventCreate(&start));
   //CHECKED_CALL(cudaEventCreate(&stop));
	CHECKED_CALL(cudaMemcpy(dev_inBox, inBox, numBoxes*2*sizeof(double), cudaMemcpyHostToDevice));

	//CHECKED_CALL(cudaEventRecord(start, 0));
	calculateLimitsOnCUDA<<<GridSize, BLOCK_SIZE>>>(dev_outLimits, dev_inBox, inRank, inFunc);
    CHECKED_CALL(cudaGetLastError());

    //CHECKED_CALL(cudaEventRecord(stop, 0));
    CHECKED_CALL(cudaDeviceSynchronize());

    CHECKED_CALL(cudaMemcpy(outLimits, dev_outLimits, numBoxes*3*sizeof(double), cudaMemcpyDeviceToHost));

	const double *boxes = inBox;
	for(int n = 0; n < numBoxes; n++)
		{
			for(int j=0; j < inRank; j++)
			{
				std::cout << "[" << boxes[(n*inRank+j)*2] << ";" << boxes[(n*inRank+j)*2 + 1] << "] ";
			}
			
			std::cout<< "\n";
		}

	float time;
   // CHECKED_CALL(cudaEventElapsedTime(&time, start, stop));

   // CHECKED_CALL(cudaEventDestroy(start));
    //CHECKED_CALL(cudaEventDestroy(stop));
    CHECKED_CALL(cudaFree(dev_outLimits));
    CHECKED_CALL(cudaFree(dev_inBox));

}









#endif