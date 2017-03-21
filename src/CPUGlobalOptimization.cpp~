#include "CPUGlobalOptimization.h"
#include "interval.h"


/**
*	Calculus Interval for Multiple function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
void fnCalcFunLimitsMultiple2(double *inBox, int inRank, double *outLimits)
{
	double x1 = (inBox[0]+inBox[1])/2;
	double x2 = (inBox[2]+inBox[3])/2;

	double var1 = inBox[0]*inBox[2];
	double var2 = inBox[0]*inBox[3];
	double var3 = inBox[1]*inBox[2];
	double var4 = inBox[1]*inBox[3];

	outLimits[0] = std::min(std::min(var1,var2),std::min(var3,var4));
	outLimits[1] = std::max(std::max(var1,var2),std::max(var3,var4));
	outLimits[2] = x1*x2;
}

/**
*	Calculus Interval for Hyperbolic function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
void fnCalcFunLimitsHypebolic2(double *inBox, int inRank, double *outLimits)
{
	double limits[2];
	double limits2[2];

	double x1 = (inBox[0]+inBox[1])/2;
	double x2 = (inBox[2]+inBox[3])/2;

	double var1 = inBox[0]*inBox[0];
	double var2 = inBox[0]*inBox[1];
	double var3 = inBox[1]*inBox[0];
	double var4 = inBox[1]*inBox[1];

	limits[0] = var2 < 0 ? 0 : std::min(std::min(var1,var2),std::min(var3,var4));
	limits[1] = std::max(std::max(var1,var2),std::max(var3,var4));

	var1 = inBox[2]*inBox[2];
	var2 = inBox[2]*inBox[3];
	var3 = inBox[3]*inBox[2];
	var4 = inBox[3]*inBox[3];

	limits2[0] = var2 < 0 ? 0 : std::min(std::min(var1,var2),std::min(var3,var4));
	limits2[1] = std::max(std::max(var1,var2),std::max(var3,var4));

	outLimits[0] = limits[0] - limits2[1];
	outLimits[1] = limits[1] - limits2[0];

	outLimits[2] = x1*x1-x2*x2;
}

/**
*	Calculus Interval for AluffiPentini function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
void fnCalcFunLimitsAluffiPentini2(double *inBox, int inRank, double *outLimits)
{
		double limits[2];
	double limits2[2];

	double x1 = (inBox[0]+inBox[1])/2;
	double x2 = (inBox[2]+inBox[3])/2;

	double var1 = inBox[0]*inBox[0];
	double var2 = inBox[0]*inBox[1];
	double var3 = inBox[1]*inBox[0];
	double var4 = inBox[1]*inBox[1];

	limits[0] = var2 < 0 ? 0 : std::min(std::min(var1,var2),std::min(var3,var4));
	limits[1] = std::max(std::max(var1,var2),std::max(var3,var4));

	var1 = inBox[2]*inBox[2];
	var2 = inBox[2]*inBox[3];
	var3 = inBox[3]*inBox[2];
	var4 = inBox[3]*inBox[3];

	limits2[0] = limits[0]*limits[0];
	limits2[1] = limits[1]*limits[1];

	outLimits[0] = 0.25*limits2[0] - 0.5*limits[1] + 0.1*inBox[0] + 0.5*(var2 < 0 ? 0 : std::min(std::min(var1,var2),std::min(var3,var4)));
	outLimits[1] = 0.25*limits2[1] - 0.5*limits[0] + 0.1*inBox[1] + 0.5*(std::max(std::max(var1,var2),std::max(var3,var4)));
	outLimits[2] = 0.25*pow(x1,4)-0.5*pow(x1,2) + 0.1*x1 + 0.5*pow(x2,2);
}


 double fnCalcFunRozenbroke(double *inBox, int inRank)
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
*	Calculus Interval for Rozenbroke function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param outlimits pointer to estimated function limits
*/
void fnCalcFunLimitsRozenbroke(double *inBox, int inRank, double *outLimits)
{

	double sup = 0;
	double sub = 0;

	double sup1,sub1,sup2,sub2, a,b, val = 0,var1,var2,var3,x1,x2;


	for(int i = 0; i < inRank - 1; i++)
	{
		sub1 = 1 - inBox[i*2 + 1];
		sup1 = 1 - inBox[i*2];

		var1 = sup1*sup1;
		var2 = sup1*sub1;
		var3 = sub1*sub1;
		
		if(sub1*sup1 < 0) sub1 = 0;
		else sub1 = std::min(std::min(var1,var2),var3);
		sup1 = std::max(std::max(var1,var2),var3);


		var1 = inBox[i*2 + 1]*inBox[i*2 + 1];
		var2 = inBox[i*2 + 1]*inBox[i*2];
		var3 = inBox[i*2]*inBox[i*2];

		if(inBox[i*2 + 1]*inBox[i*2] < 0) a = 0;
		else a = std::min(std::min(var1,var2),var3);
		b = std::max(std::max(var1,var2),var3);

		sub2 = inBox[(i+1)*2] - b;
		sup2 = inBox[(i+1)*2 + 1] - a;

		var1 = sup2*sup2;
		var2 = sup2*sub2;
		var3 = sub2*sub2;

		if(sub2*sup2 < 0) sub2 = 0;
		else sub2 = 100*std::min(std::min(var1,var2),var3);
		sup2 = 100*std::max(std::max(var1,var2),var3);

		sub += sub1 + sub2;
		sup += sup1 + sup2;

		x1 = (inBox[i*2 + 1] + inBox[i*2])/2;
		x2 = (inBox[(i+1)*2 + 1] + inBox[(i+1)*2])/2;
		val += ((1 - x1)*(1 - x1) + 100*(x2-x1*x1)*(x2-x1*x1));

	}

	outLimits[0] = sub;
	outLimits[1] = sup;
	outLimits[2] = val;


	double *x;
	x = (double *) malloc(inRank*sizeof(double));
	double minFun;
	int index;
	double md;
	double pw;
	double a1;
	double b1;
	int numFunValues = int(pow(2.0,double(inRank)));

	
	int i,j;
	
	for(j = 0; j < inRank; j++){
			x[j] = (inBox[j*2]+inBox[j*2+1])/2.0;
	}
	minFun = fnCalcFunRozenbroke(x, inRank);
	for(i = 0; i < numFunValues; i++){
		for(j = 0; j < inRank; j++){
			a1 = double(i);
			b1 = double(j+1);
			md = fmod(a1,pow(2.0,b1));
			pw = pow(2.0,(b1-1.0));
			index = int(md / pw);
			x[j] = inBox[j*2+index];
		}
		val = fnCalcFunRozenbroke(x, inRank);
		if(minFun > val) minFun = val;
	}

	free(x);

	outLimits[0] = sub;
	outLimits[1] = sup;
	outLimits[2] = minFun;



}


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
void fnGetOptValueOnCPU(double *inBox, int inRank, int inNumBoxesSplitCoeff, double inEps, double inMaxIter, void (*inFun)(double *,int,double *), double *outBox, double*outMin, double *outEps,int *outStatus)
{
	int numBoxes = pow((double) inNumBoxesSplitCoeff,inRank);
	double *boxes =  new double[numBoxes*inRank*2];
	double *boxesResult = new double[numBoxes*3];
	double *restBoxes = new double[inRank*2];
	double *tempRestBoxes = NULL;
	double *h = new double[inRank];
	int numNewBoxes = 0;

	memcpy(restBoxes,inBox,inRank*2*sizeof(double));

	int numRestBoxes = 1;
	int index = 0;
	double temp;


	int countIter = 0;
	double curEps = inEps*10;
	
	inFun(inBox,inRank,boxesResult);
	double funcMin = boxesResult[2];
	int boxMinIndex;
	boxMinIndex = 0;

	int i,j,k,n;

	*outStatus = 1;

	

	while((countIter < inMaxIter) && (curEps >= inEps))
	{
		boxes = new double[numRestBoxes*numBoxes*inRank*2];
		boxesResult = new double[numRestBoxes*numBoxes*3];

		for(k = 0; k < numRestBoxes; k++)
		{
			for(i = 0; i < inRank; i++)
			{
				h[i] = (restBoxes[(k*inRank+i)*2 + 1] - restBoxes[(k*inRank+i)*2])/inNumBoxesSplitCoeff;
			}
			for(n = 0; n < numBoxes; n++)
			{
				for(i = 0; i < inRank; i++)
				{
					index = ((n % numBoxes) % (long) pow((double)inNumBoxesSplitCoeff,i+1))/((long)pow((double)inNumBoxesSplitCoeff,i));
					boxes[((k*numBoxes + n)*inRank+i)*2] = restBoxes[(k*inRank+i)*2] + h[i]*index;
					boxes[((k*numBoxes + n)*inRank+i)*2 + 1] = restBoxes[(k*inRank+i)*2] + h[i]*(index+1);
				}
				inFun(&boxes[((k*numBoxes + n)*inRank)*2],inRank,&boxesResult[(k*numBoxes + n)*3]);
			}
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
			if(funcMin < boxesResult[(n)*3] + inEps) break;	
		}
		curEps = abs(boxesResult[0] - funcMin);
		*outEps = curEps;
		*outMin = funcMin;
		memcpy(outBox,boxes + n*inRank*2,inRank*2*sizeof(double));
		std::cout << boxesResult[0] << "\t" << boxesResult[n*3] << "\t" << funcMin << "\t" << curEps << "\t" << n << "\n\n";
		if(curEps < inEps)
		{
			*outStatus = 0;
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


