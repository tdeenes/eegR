// This file is part of the program ept_TFCE.
// ept_TFCE is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ept_TFCE is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with ept_TFCE.  If not, see <http://www.gnu.org/licenses/>.

/*
 * Christian Gaser
 *
 * Pau Coma 09.12.2010
 *
 * -Independent grow_i, grow_j, grow_k variables to increase understandability
 *
 * Armand Mensen 20.12.2010
 *
 * - Added the input which indicates channels neighbours
 * - Searches the 2D surface coordinates for neighbours
 *
 */

#include <stdlib.h>
#include <R.h>

//#ifdef OPENMP
//#include "omp.h"
//#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

void tfce_thread(double *inData, double *outData, double *ChN, double thresh, double delta, const int *dims, const int *dims2, double *EH)
{
   double valToAdd;
   double E, H;
   int i, t, ti, tt, maxt, mint, temp, growingInd, growingCur, ChCurr, idx;
   int numVoxels = dims[0] * dims[1];
   char* flagUsed;
   //#-short* growing;
   short* grow_i;
   short* grow_t;
   //-#
   
   E = EH[0];
   H = EH[1];
   
   
   flagUsed = (char*)malloc(numVoxels*sizeof(char));
   //#-growing  = (short*)malloc(numVoxels*3*sizeof(short));
   grow_i  = (short*)malloc(numVoxels*sizeof(short));
   grow_t  = (short*)malloc(numVoxels*sizeof(short));
   //-#
   

    for (temp = 0; temp < numVoxels; ++temp) flagUsed[temp] = 0;
    
	for (t = 0; t < dims[1]; ++t)
	{
		for (i = 0; i < dims[0]; ++i)
		{
			// temp is the current point in the grid
			temp = (t*dims[0]) + i;
			// Check if this point has been seen (flagUsed) and whether its over the threshold
			if (!flagUsed[temp] && inData[temp] >= thresh)
			{
				// make the flagUsed so that algorithm doesn't visit this point again
				flagUsed[temp] = 1;
				growingInd = 1;
				growingCur = 0;
				
				   // mexPrintf("\nCurrent temp value is %f \n", temp);
				   // mexEvalString("drawnow;"); // to dump string.
				
				// Define the current coordinates to continue with
				grow_i[0]=i;
				grow_t[0]=t;

				while (growingCur < growingInd)
				{
				   //This just limits not to overrun borders <-- in our case the zero padding
				   //And creates a 3x3 windows of scanning for valid points (one point around the current point)
				   //E.g. In 2 dimensions... if the current point is 2,j then maxi = Min of dimension size or 2+2
				   //thus maxi = 4... mini = max of 0 or 2-1... so 1...
				   //so in the next set of loops ti will look at i=1,i=2,and i=3 (not 4 because its set to < maxi
				   maxt = MIN(dims[1], grow_t[growingCur] + 2);
				   mint = MAX(0, grow_t[growingCur] - 1);
					
					// start of the smaller scanning window
					
					for (tt = mint; tt < maxt; ++tt)
					{
						for (ti = 0; ti < dims2[1]; ++ti)
						{
							idx = (ti*dims2[0]) + grow_i[growingCur];
							ChCurr = ChN[idx];
							
							// mexPrintf("\nCurrent Point is (%d,%d) \n", ChCurr, tt+1);
							// mexEvalString("drawnow;"); // to dump string.
							
							if (ChCurr == 0)
							{
								break;
							}
							
							ChCurr = ChCurr - 1;
							
							temp = (tt*dims[0]) + ChCurr;
							
							// mexPrintf("\nLooking at value %f \n", inData[temp]);
							// mexEvalString("drawnow;"); // to dump string.
							
							if (!flagUsed[temp] && inData[temp] >= thresh)
							{
								flagUsed[temp] = 1;
								grow_i[growingInd] = ChCurr;
								grow_t[growingInd] = tt;
								growingInd++;
								//Here the growing index increases everytime we find a value above the threshhold
								//This grows depending on the origin of our supporting weight which is the outer loop
								//-#
							}
						}	
					}
				   // GrowingCur increases one and thus looks at the next point that was found in the previous loop
				   growingCur++;

				}
				// Reset back to 0 so that when the next "while loop" runs it adds the value to each point found
				growingCur = 0;
				
							// mexPrintf("\nCluster size is %d \n", growingInd);
							// mexEvalString("drawnow;"); // to dump string.
							
				valToAdd = pow(growingInd, E) * pow(thresh, H) * delta;
				
				// Adds the valToAdd to the points that it found in that cluster
				while (growingCur < growingInd)
				{
				   outData[(grow_t[growingCur]*dims[0]) + grow_i[growingCur]] += valToAdd;
				   // GrowingCur adds one so that next coordinate is added until growingCurr is less than the number of points found
				   growingCur++;
				}
			}
		}
	}
   
   free(flagUsed);
   
   free(grow_i);
   free(grow_t);
}

void tfce(double *inData, double *outData, double *ChN, const int *dims, const int *dims2, double *EH)
{
   double fmax = 0.0, thresh0, thresh, delta;
   int i, numSteps = 50;
   int numVoxels = dims[0] * dims[1];
   
   for (i = 0; i < numVoxels; ++i)
   {
      if (inData[i] > fmax) fmax = inData[i];
      outData[i] = 0.0;
   }

   delta = fmax/numSteps;
   thresh0 = delta/2.0;
   
  {
	//#ifdef OPENMP
	//	/* this is faster than the default setting to the # of processors */
	//	omp_set_num_threads(numSteps);
	//	# pragma omp parallel for private(thresh) shared(outData) 
	//#endif
    for (i = 0; i < numSteps; i++)
    {
       thresh = thresh0 + (double)i*delta;
	   
	    // mexPrintf("\n Current Threshold is %f... \n", thresh);
		// mexEvalString("drawnow;"); // to dump string.
		
		// mexPrintf("\n Current Delta is %f... \n", delta);
		// mexEvalString("drawnow;"); // to dump string.
			   
       tfce_thread(inData, outData, ChN, thresh, delta, dims, dims2, EH);
    }
  }   
}

