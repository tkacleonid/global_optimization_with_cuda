#ifndef __CPUGlobalOptimization_H__
#define __CPUGlobalOptimization_H__

/**
*	Calculus Interval for Multiple function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
void fnCalcFunLimitsMultiple2(double *inBox, int inRank, double *outLimits);
/**
*	Calculus Interval for Hyperbolic function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
void fnCalcFunLimitsHypebolic2(double *inBox, int inRank, double *outLimits);
/**
*	Calculus Interval for AluffiPentini function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
void fnCalcFunLimitsAluffiPentini2(double *inBox, int inRank, double *outLimits);
/**
*	Calculus Interval for Rozenbroke function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
void fnCalcFunLimitsRozenbroke(double *inBox, int inRank, double *outLimits);
/**
*	Calculus minimum value for function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param inNumBoxesSplitCoeff number of parts for each dimension
*	@param inEps required accuracy
*	@param inMaxIter maximum count of iterations
*	@param inFun pointer to optimazing function
*	@param outBox pointer to optimal box
*	@param outMin pointer to optimal value
*	@param outEps pointer to reached accuracy
*	@param outEps pointer to status of solving optimization problem
*/
void fnGetOptValueOnCPU(double *inBox, int inRank, int inNumBoxesSplitCoeff, double inEps, double inMaxIter, void (*inFun)(double *,int,double *), double *outBox, double*outMin, double *outEps,int *outStatus);



#endif