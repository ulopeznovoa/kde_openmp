/**************

computePDF.c

Functions used to compute the PDFs for several dimensionalities.


Copyright (c) 2014, Unai Lopez-Novoa, Jon Saenz, Alexander Mendiburu 
and Jose Miguel-Alonso  (from Universidad del Pais Vasco/Euskal 
		    Herriko Unibertsitatea)

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Universidad del Pais Vasco/Euskal 
      Herriko Unibertsitatea  nor the names of its contributors may be 
      used to endorse or promote products derived from this software 
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

***************/

#include "computePDF.h"
#include "linalg.h"

double volumeConstant(int dim)
{
	if(dim == 1)
		return 2.;
	else if(dim == 2)
		return acos(-1.);
	else if (dim == 3)
		return acos(-1.)*4./3.;	
	else
		return unit_sphere_volume(dim);
}
 
/**** Functions to calculate the PDF of a defined 2D space (box) for a given sample. ****/ 
 
//Compute the density in the bounding box of a sample - Function for 2D spaces
void compute2DBox_2D(PDFPtr pdf, double * PC, double * lower, int * tot_ev_per_dim, double * gridpoint, size_t * dif_pos, 
	double * x0, double * dx, double h2, double cd, MAT * eigenvectors, double * restrict densValues, int * restrict densPosition)
{
	int u,v,l; //Loop variables
	double temp; //Will contain the absolute distance value from gridpoint to sample.
	double PCdot[2] __attribute__((aligned(64)));
	
	for(gridpoint[0] = lower[0], u = 0; u < tot_ev_per_dim[0]; gridpoint[0] += dx[0], u++)
	{  
		int HalfPosition = (((gridpoint[0] - x0[0])/ dx[0]) * pdf->pdfcumsize[0]);
				
		//Compiler flag to inform about structure alignment		
		__assume_aligned(densValues,64);
		__assume_aligned(densPosition,64);				
				
		#ifdef _OPENMP
		#pragma simd private(PCdot,temp) assert
		#endif
		for(v = 0; v < tot_ev_per_dim[1]; v++)
		{  	
			//Conversion to PC space			
			PCdot[0] = (eigenvectors->me[0][0] * gridpoint[0]) + (eigenvectors->me[0][1] * (lower[1] + (dx[1] * v)));
		    	PCdot[1] = (eigenvectors->me[1][0] * gridpoint[0]) + (eigenvectors->me[1][1] * (lower[1] + (dx[1] * v)));		    
		
			//Absolute distance calculation
			temp = (((PC[0] - PCdot[0]) * (PC[0] - PCdot[0])) + ((PC[1] - PCdot[1]) * (PC[1] - PCdot[1])) ) / h2;
					
			//If OpenMP version, store the density value in an auxiliar vector densValues, previous to storing in the final PDF structure
			//Vector densPosition will contain the position of the gridpoint in the final PDF structure
			#ifdef _OPENMP

			//PDFposition			
			densPosition[v] = HalfPosition + ((((lower[1] + (dx[1] * v)) - x0[1])/ dx[1]) * pdf->pdfcumsize[1]);		  
			densValues[v] = (0.5/cd*(2+2.)*(1.-temp)) * (fabs(temp)<1.);

			//If serial version, store the density value of the sample over the gridpoint in the PDF structure
			#else

			gridpoint[1] = (lower[1] + (dx[1] * v));

	        	dif_pos[0] = (gridpoint[0] - x0[0])/ dx[0];
			dif_pos[1] = (gridpoint[1] - x0[1])/ dx[1];

			*PDFitem(pdf ,dif_pos, 2) += (0.5/cd*(2+2.)*(1.-temp)) * (fabs(temp)<1.) ;	
			
			#endif	
		}

		#ifdef _OPENMP

		for(v = 0; v < tot_ev_per_dim[1]; v++)
			#pragma omp atomic
			pdf->PDF[densPosition[v]] += densValues[v];
	
		#endif
	}
}

void compute2DBox_3D(PDFPtr pdf, double * PC, double * lower, int * tot_ev_per_dim, double * gridpoint,	size_t * dif_pos, 
	double * x0,double * dx, double h2, double cd, MAT * eigenvectors, double * restrict densValues, int * restrict densPosition)
{
	int u,v,l; //Loop variables
	double temp; //Will contain the absolute distance value from gridpoint to sample.
	double PCdot[3] __attribute__((aligned(64)));
	
	//Compiler flag to inform about structure alignment		
	__assume_aligned(densValues,64);
	__assume_aligned(densPosition,64);				

	for(gridpoint[0] = lower[0], u = 0; u <= tot_ev_per_dim[0]; gridpoint[0] += dx[0], u++)
	{  
		int HalfPosition = (((gridpoint[0] - x0[0])/ dx[0]) * pdf->pdfcumsize[0]) + (((gridpoint[2] - x0[2])/ dx[2]) * pdf->pdfcumsize[2]);
		
		#pragma simd private(PCdot) assert	
		for(v = 0; v <= tot_ev_per_dim[1]; v++) 
		{  	
			//Conversion to PC space			
		    	PCdot[0] = (eigenvectors->me[0][0] * gridpoint[0]) + (eigenvectors->me[0][1] * (lower[1] + (dx[1] * v))) + (eigenvectors->me[0][2] * gridpoint[2]);
		    	PCdot[1] = (eigenvectors->me[1][0] * gridpoint[0]) + (eigenvectors->me[1][1] * (lower[1] + (dx[1] * v))) + (eigenvectors->me[1][2] * gridpoint[2]);		    
		    	PCdot[2] = (eigenvectors->me[2][0] * gridpoint[0]) + (eigenvectors->me[2][1] * (lower[1] + (dx[1] * v))) + (eigenvectors->me[2][2] * gridpoint[2]);		
		
			//Absolute distance calculation
			temp = (((PC[0] - PCdot[0]) * (PC[0] - PCdot[0])) + ((PC[1] - PCdot[1]) * (PC[1] - PCdot[1])) + ((PC[2] - PCdot[2]) * (PC[2] - PCdot[2]))) / h2;
					
			//If OpenMP version, store the density value in an auxiliar vector densValues, previous to storing in the final PDF structure
			//Vector densPosition will contain the position of the gridpoint in the final PDF structure
			#ifdef _OPENMP

			//PDFposition			
			densPosition[v] = HalfPosition + ((((lower[1] + (dx[1] * v)) - x0[1])/ dx[1]) * pdf->pdfcumsize[1]);		  
			densValues[v] = (0.5/cd*(3+2.)*(1.-temp)) * (fabs(temp)<1.);

			//If serial version, store the density value of the sample over the gridpoint in the PDF structure
			#else

			gridpoint[1] = (lower[1] + (dx[1] * v));

		        dif_pos[0] = (gridpoint[0] - x0[0])/ dx[0];
			dif_pos[1] = (gridpoint[1] - x0[1])/ dx[1];
			dif_pos[2] = (gridpoint[2] - x0[2])/ dx[2];

			*PDFitem(pdf ,dif_pos, 3) += (0.5/cd*(3+2.)*(1.-temp)) * (fabs(temp)<1.) ;	
			
			#endif	

		}

		#ifdef _OPENMP

		for(v = 0; v <= tot_ev_per_dim[1]; v++)
			#pragma omp atomic
			pdf->PDF[densPosition[v]] += densValues[v];
	
		#endif
	}
}

//Compute the density in the bounding box of a sample - Generic function, used for spaces of dimensionality higher than 3
void compute2DBox_ND(PDFPtr pdf, double * PC, double * lower, int * tot_ev_per_dim, double * gridpoint, size_t * dif_pos, double * x0, 
	double * dx, int dim, double h2, double cd, MAT * eigenvectors, double * restrict densValues, int * restrict densPosition, 
	double * restrict PCdot_vec, double * restrict temp_vec, double * restrict gridpoint_vec)
{
	int u,v,d,l; //Loop variables
	
	int HalfPosition;
	int dimGreaterThanTwoPosition = 0;
	double HalfTemp = 0;

	#ifdef _OPENMP //Initializations for vector implementation

	#pragma simd reduction(+:dimGreaterThanTwoPosition) assert
	for(d = 2; d < dim; d++)
		dimGreaterThanTwoPosition += (dif_pos[d] * pdf->pdfcumsize[d]);
		
	for(v = 0; v < tot_ev_per_dim[1]; v++) 
		for(d = 2; d < dim; d++)
			gridpoint_vec[v * dim + d] = gridpoint[d];
	#endif
	
	for(gridpoint[0] = lower[0], u = 0; u < tot_ev_per_dim[0]; gridpoint[0] += dx[0], u++)
	{  		
		//Compiler flag to inform about structure alignment		
		__assume_aligned(densValues,64);
		__assume_aligned(densPosition,64);		
		
		#ifdef _OPENMP //Vector friendly implementation
		
		HalfPosition = (((gridpoint[0] - x0[0])/ dx[0]) * pdf->pdfcumsize[0]) + dimGreaterThanTwoPosition;		
				
		for(v = 0; v < tot_ev_per_dim[1]; v++) 
			gridpoint_vec[v * dim + 0] = gridpoint[0];
			
		for(v = 0; v < tot_ev_per_dim[1]; v++) 	
			temp_vec[v] = 0;
		
		for(v = 0; v < tot_ev_per_dim[1]; v++) 
			gridpoint_vec[v * dim + 1] = (lower[1] + (dx[1] * v));
		
		for(v = 0; v < tot_ev_per_dim[1] * dim; v++) 
			PCdot_vec[v] = 0;	
				
		for(v = 0; v < tot_ev_per_dim[1]; v++) 
			for(d = 0; d < dim; d++)	
				#pragma simd reduction(+:PCdot_vec[v * dim + d]) assert
				for(l = 0; l < dim; l++)
					PCdot_vec[v * dim + d] += (eigenvectors->me[d][l] * gridpoint_vec[v * dim + l]);		
						
		for(v = 0; v < tot_ev_per_dim[1]; v++) 			
			#pragma simd reduction(+:temp_vec[v]) assert
			for(d = 0; d < dim; d++)				
				temp_vec[v] += ((PC[d] - PCdot_vec[v * dim + d]) * (PC[d] - PCdot_vec[v * dim + d]));
		
		for(v = 0; v < tot_ev_per_dim[1]; v++) 		
			temp_vec[v] /= h2;
			
		for(v = 0; v < tot_ev_per_dim[1]; v++) 
			densPosition[v] = HalfPosition + ((((lower[1] + (dx[1] * v)) - x0[1])/ dx[1]) * pdf->pdfcumsize[1]);		  
		
		for(v = 0; v < tot_ev_per_dim[1]; v++) 
			densValues[v] = (0.5/cd*(dim + 2.)*(1.-temp_vec[v])) * (fabs(temp_vec[v])<1.);
			
		for(v = 0; v < tot_ev_per_dim[1]; v++)
			#pragma omp atomic
			pdf->PDF[densPosition[v]] += densValues[v];		
		
		#else	// Serial implementation
		
		double temp;		
		dif_pos[0] = (gridpoint[0] - x0[0])/ dx[0];		
		
		for(v = 0; v < tot_ev_per_dim[1]; v++) 
		{  						
			gridpoint[1] = (lower[1] + (dx[1] * v));
			
			//Conversion to PC space						
			for(d = 0; d < dim; d++)
				PCdot_vec[d] = 0;	
			
			for(d = 0; d < dim; d++)	
				#pragma simd reduction(+:PCdot_vec[d]) assert
				for(l = 0; l < dim; l++)
					PCdot_vec[d] += (eigenvectors->me[d][l] * gridpoint[l]);
				
			//Absolute distance calculation	
			temp = 0;
			
			#pragma simd reduction(+:temp) assert
			for(d = 0; d < dim; d++)
				temp += ((PC[d] - PCdot_vec[d]) * (PC[d] - PCdot_vec[d]));
				
			temp /= h2;	
			
 			dif_pos[1] = (gridpoint[1] - x0[1])/ dx[1];

			*PDFitem(pdf ,dif_pos, dim) += (0.5/cd*(dim + 2.)*(1.-temp)) * (fabs(temp)<1.) ;				
		}		
		
		#endif		
	}
}

/**** Functions to calculate PDF, called from main ****/

//Compute the PDF of a one-dimensional grid space
void computePDF1D(MPDFEstimatorPtr mpdf, PDFPtr pdf, MAT *Sm1 , double h , double detSm1 , double *x0, 
		double *x1, double *dx, double *bounds, MAT *eigenvectors )
{
  int i,j,u; //Loop variables
  int dim = 1;	//Dimensions of grid space
  double cd = volumeConstant(dim); //Volume constants to calculate kernel values    
  double h2=h*h;  //Squared bandwith value
  double *PC; // Current sample (PC space)
  double theintegral = 0.0;  
  double total_vol = 0.0;  
  double * sample; 
  double k=1./sqrt(detSm1)/mpdf->current/pow(h,mpdf->length);  //Constant to recover the volume in the X space from the volume in the PC space
  double PCdot;
   
 //Variables to calculate coordinates and number of gridpoints of bounding box
  int steps;	  
  double upper, lower, gridpoint;
  int tot_ev;
  size_t dif_pos[1];
  double abs_bound,temp;   
  
  //Auxiliary vectors for OpenMP version
  double * densValues;
  int * densPosition;    
        
  #pragma omp parallel default(none) \
  shared(stdout,mpdf,pdf,dim,x0,x1,dx,theintegral,total_vol,bounds,eigenvectors,cd,h2,k) \
  private(i,j,u,sample,PC,lower,upper,steps,abs_bound,tot_ev,dif_pos,gridpoint,PCdot,densValues,densPosition,temp) 
  {	
 
  #ifdef _OPENMP 

  int dim0_max_size = ((ceil(bounds[0] / dx[0]) * 2) + 3);

  densValues = (double *)_mm_malloc(sizeof(double) * dim0_max_size,64); //Vector to hold density values of each sample-gridpoint combination
  densPosition = (int *)_mm_malloc(sizeof(int) * dim0_max_size,64); //Vector to hold the positions of densValues values in the PDF structure

  #endif     
  
  //Initialize PDF structure to 0s
  #pragma omp for
  for(i = 0; i < pdf->total_size; i++)
	pdf->PDF[i] = 0.0f;  
  
  //Main calculation loop. For each sample calculate the PDF of its influence area and store in the PDF structure
  #pragma omp for
  for(i=0;i<mpdf->current;i++) 
  {	
	sample = MPDFPosition(mpdf,i); //Get current sample
	PC = MPDFPCPosition(mpdf,i); //Get current sample (scaled as PC)
	  	  
	//For each sample, calculate its boundaries
	
	//Lower corner
	abs_bound = sample[0] - bounds[0];
	if (x0[0] > abs_bound)
		lower = x0[0];
	else
	{
		steps = floor((abs_bound - x0[0]) / dx[0]);
		lower = x0[0] + (steps * dx[0]);
	}

	//Upper corner
	abs_bound = sample[0] + bounds[0];
	if (x1[0] < abs_bound)
		upper = x1[0];
	else
	{
		steps = ceil((abs_bound - x0[0]) / dx[0]);
		upper = x0[0] + (steps * dx[0]);
	}	
	
	//Calculate number of eval points per dimension	
	tot_ev = rint((upper - lower)/dx[0]) + 1;		   
  

	//Calculate the PDF of the defined 1D space
	#ifdef _OPENMP
	#pragma simd private(PCdot,temp) assert
	#endif
	for(u = 0; u < tot_ev; u++)
	{  		
	    PCdot = (eigenvectors->me[0][0] * (lower + (dx[0] * u)));
	
		//Absolute distance calculation
		temp = ((PC[0] - PCdot) * (PC[0] - PCdot)) / h2;

		//If OpenMP version, store the density value in an auxiliar vector densValues, previous to storing in the final PDF structure
		//Vector densPosition will contain the position of the gridpoint in the final PDF structure
		#ifdef _OPENMP

		//PDFposition			
		densPosition[u] = (((lower + (dx[0] * u)) - x0[0])/ dx[0]) * pdf->pdfcumsize[0];		  
		densValues[u] = (0.5/cd*(1+2.)*(1.-temp)) * (fabs(temp)<1.);

		//If serial version, store the density value of the sample over the gridpoint in the PDF structure
		#else

		dif_pos[0] = ((lower + (dx[0] * u)) - x0[0])/ dx[0];
		*PDFitem(pdf ,dif_pos, 1) += (0.5/cd*(1+2.)*(1.-temp)) * (fabs(temp)<1.) ;	
		
		#endif			
	}	
	
	#ifdef _OPENMP

	for(u = 0; u < tot_ev; u++)
		#pragma omp atomic
		pdf->PDF[densPosition[u]] += densValues[u];

	#endif	    	 
  } 
    
  #ifdef _OPENMP
  _mm_free(densValues);
  _mm_free(densPosition);      
  #endif

  //Apply k constant to PDF
  #pragma omp for
  for(i=0; i < pdf->total_size; i++)
  	  pdf->PDF[i] = pdf->PDF[i] * k;

  //Calculate integral of PDF
  #pragma omp for reduction(+:theintegral)  
  for(i=0; i < pdf->total_size; i++)
      theintegral += pdf->PDF[i];
      
  #pragma omp single
  theintegral = theintegral * dx[0];     

  //Renormalize PDF using integral
  #pragma omp for
  for(i=0; i < pdf->total_size; i++)
     pdf->PDF[i] = pdf->PDF[i]/theintegral;
     
  //Calculate total volume of renormalized PDF   
  #pragma omp for reduction(+:total_vol)  
  for(i=0; i < pdf->total_size; i++)
	total_vol += pdf->PDF[i];
  
  }//End of parallel OpenMP Region    
  
  printf("Total integrated PDF: %g. The integral: %f\n",total_vol*dx[0],theintegral);	   
} 

//Compute the PDF of a 2D grid space
void computePDF2D(MPDFEstimatorPtr mpdf, PDFPtr pdf, MAT *Sm1 , double h , double detSm1 , double *x0, 
		double *x1, double *dx, double *bounds, MAT *eigenvectors )
{
  int i,j; //Loop variables
  int dim = 2;	//Dimensions of grid space
  double cd = volumeConstant(dim); //Volume constants to calculate kernel values    
  double h2=h*h;  //Squared bandwith value
  double *PC; // Current sample (PC space)
  double theintegral = 0.0;  
  double total_vol = 0.0;  
  double total_dx = dx[0] * dx[1]; 
  double * sample; 
  double k=1./sqrt(detSm1)/mpdf->current/pow(h,mpdf->length);  //Constant to recover the volume in the X space from the volume in the PC space
  double * PCdot;
 
  //Variables to calculate coordinates and number of gridpoints of bounding box
  int steps;	  
  double upper, lower[2], gridpoint[2];
  int tot_ev_per_dim[2];
  size_t dif_pos[2];
  double abs_bound;

  //Auxiliary vectors for OpenMP version
  double * densValues;
  int * densPosition;    
              
  #pragma omp parallel default(none) \
  shared(mpdf,pdf,dim,x0,x1,dx,total_dx,theintegral,total_vol,bounds,eigenvectors,cd,h2,k) \
  private(i,j,sample,PC,lower,upper,steps,abs_bound,tot_ev_per_dim,dif_pos,gridpoint,PCdot,densValues,densPosition) 
  {	
 	
  #ifdef _OPENMP 

  int dim1_max_size = ((ceil(bounds[1] / dx[1]) * 2) + 3);

  densValues = (double *)_mm_malloc(sizeof(double) * dim1_max_size,64); //Vector to hold density values of each sample-gridpoint combination
  densPosition = (int *)_mm_malloc(sizeof(int) * dim1_max_size,64); //Vector to hold the positions of densValues values in the PDF structure

  #endif   
 
  //Initialize PDF structure to 0s  
  #pragma omp for
  for(i = 0; i < pdf->total_size; i++)
	pdf->PDF[i] = 0.0f;

  //Main calculation loop. For each sample calculate the PDF of its influence area and store in the PDF structure
  #pragma omp for
  for(i=0;i<mpdf->current;i++) 
  {	  
	sample = MPDFPosition(mpdf,i); //Get current sample
	PC = MPDFPCPosition(mpdf,i); //Get current sample (scaled as PC)

	//For each sample, calculate its bounding box, 
	//expressed as coordinates of lower corner and number of gridpoints per dimensions
	for(j = 0; j < 2; j++)
	{	
		//Lower corner
		abs_bound = sample[j] - bounds[j];
		if (x0[j] > abs_bound)
			lower[j] = x0[j];
		else
		{
			steps = floor((abs_bound - x0[j]) / dx[j]);
			lower[j] = x0[j] + (steps * dx[j]);
		}

		//Upper corner
		abs_bound = sample[j] + bounds[j];
		if (x1[j] < abs_bound)
			upper = x1[j];
		else
		{
			steps = ceil((abs_bound - x0[j]) / dx[j]);
			upper = x0[j] + (steps * dx[j]);
		}	
		
		//Calculate number of eval points per dimension	
		tot_ev_per_dim[j] = rint((upper - lower[j])/dx[j]) + 1;			
	}    

	//Calculate the PDF of the defined 2D box	
	compute2DBox_2D(pdf,PC,lower,tot_ev_per_dim,gridpoint,dif_pos,x0,dx,h2,cd,eigenvectors,densValues,densPosition);	
  }
       
  #ifdef _OPENMP
  _mm_free(densValues);
  _mm_free(densPosition);      
  #endif

  //Apply k constant to PDF
  #pragma omp for
  for(i=0; i < pdf->total_size; i++)
  	  pdf->PDF[i] = pdf->PDF[i] * k;

  //Calculate integral of PDF
  #pragma omp for reduction(+:theintegral)   
  for(i=0; i < pdf->total_size; i++)
      theintegral += pdf->PDF[i];

  #pragma omp single
  theintegral = theintegral * total_dx;

  //Renormalize PDF using integral
  #pragma omp for
  for(i=0; i < pdf->total_size; i++)
     pdf->PDF[i] = pdf->PDF[i]/theintegral;
     
  //Calculate total volume of renormalized PDF   
  #pragma omp for reduction(+:total_vol)  
  for(i=0; i < pdf->total_size; i++)
	total_vol += pdf->PDF[i];
  
  }//End of parallel OpenMP Region    
  
  printf("Total integrated PDF: %g. The integral: %f\n",total_vol*dx[0]*dx[1],theintegral);
	
}

#define DEBUG_TEMPS 1
#undef  DEBUG_TEMPS

//Compute the PDF of grid spaces of dimension 3 or higher
void computePDF3D(MPDFEstimatorPtr mpdf, PDFPtr pdf, MAT *Sm1 , double h , 
		  double detSm1 , double *x0, double *x1, double *dx, 
		  double *bounds, MAT *eigenvectors)
{
  int dim = 3; 
  int i,j,l,u,w; //Loop variables	
  double cd = volumeConstant(dim); //Volume constant
  double k=1./sqrt(detSm1)/mpdf->current/pow(h,mpdf->length); //Constant to recover the volume in the X space from the volume in the PC space  
  double h2=h*h; //Square of bandwith value
  double *PC; // Current sample (PC space)
  double total_vol=0.0;
  double theintegral=0.0;
  double * sample; //Current sample
  double * PCdot;  

  //Variables to calculate the bounding box of a sample
  double lower[3];
  double upper;
  double gridpoint[3];
  int tot_ev_per_dim[3];
  size_t dif_pos[3];
  int total_ev;	
  int steps;
  double abs_bound; //Absolute bound per sample and dimension, given by ellipsoid shape

  //Calculate acumulated volume for the grid space
  double total_dx = 1.0;
  for (i = 0; i < dim; i++)
	total_dx *= dx[i];   

  //Variables to perform the calculation of the 2D layering
  double A,B,C,F,Z,theta,cosTheta,sinTheta,X2,Y2,X,Y,XY,termY2,valor,termX2,upy,rightx,upx_rot,upy_rot,rightx_rot,righty_rot; 
  double bound[2],box_center[2],box_min[2],box_max[2],box_steps[2],box_upper[2];
      			
  //Calculate partial equations for the 2D layering    							
  A = Sm1->me[0][0];
  B = 2 * Sm1->me[0][1];
  C = Sm1->me[1][1];
  theta = atan(B/(A-C))/2;		 
  cosTheta = cos(theta);
  sinTheta = sin(theta);					
  X2 =    Sm1->me[0][0]*cosTheta*cosTheta + 2*Sm1->me[0][1]*cosTheta*sinTheta +   Sm1->me[1][1]*sinTheta*sinTheta;
  XY = -2*Sm1->me[0][0]*cosTheta*sinTheta + 2*Sm1->me[0][1]*cosTheta*cosTheta - 2*Sm1->me[0][1]*sinTheta*sinTheta + 2*Sm1->me[1][1]*cosTheta*sinTheta;
  Y2 =    Sm1->me[0][0]*sinTheta*sinTheta - 2*Sm1->me[0][1]*cosTheta*sinTheta +   Sm1->me[1][1]*cosTheta*cosTheta;		
  
  //Aux vector for OpenMP version
  double * densValues; 
  int * densPosition;
  double * temp_vec;
  double * PCdot_vec;
  double * gridpoint_vec;
 
  //Beginning of OpenMP parallel region								
  #pragma omp parallel default(none)\
  shared(stdout,theintegral,total_vol,total_dx,k,mpdf,pdf,cd,dim,bounds,x0,x1,dx,Sm1,cosTheta,sinTheta,eigenvectors,X2,XY,Y2,h2,h) \
  private(i,j,l,u,w,sample,PC,gridpoint,total_ev,abs_bound,lower,box_upper,tot_ev_per_dim,box_steps,F,X,Y,Z,termX2,termY2,upy,rightx,upx_rot,upy_rot, \
  valor,rightx_rot,righty_rot,bound,box_center,box_min,box_max,PCdot,dif_pos,steps,upper,densValues,densPosition,temp_vec,gridpoint_vec,PCdot_vec)
  {			
								
  #ifdef _OPENMP

  int dim1_max_size = ((ceil(bounds[1] / dx[1]) * 2) + 3);

  densValues = (double *)_mm_malloc(sizeof(double) * dim1_max_size,64);
  densPosition = (int *)_mm_malloc(sizeof(int) * dim1_max_size,64);
	 
  #endif

  //Initialize PDF structure to 0s
  #pragma omp for 
  for(i = 0; i < pdf->total_size; i++)
	pdf->PDF[i] = 0.0f;

  //Main calculation loop. For each sample calculate the PDF of its influence area and store in the PDF structure
  #pragma omp for  
  for(i=0;i<mpdf->current;i++) 
  {	  
	sample = MPDFPosition(mpdf,i); //Get current sample
	PC = MPDFPCPosition(mpdf,i); //X is the current sample (scaled as PC)

	//Calculate boundaries for Z axis
	//Lower corner
	abs_bound = sample[2] - bounds[2];
	if (x0[2] > abs_bound)
		lower[2] = x0[2];
	else
	{
		steps = floor((abs_bound - x0[2]) / dx[2]);
		lower[2] = x0[2] + (steps * dx[2]);
	}

	//Upper corner
	abs_bound = sample[2] + bounds[2];
	if (x1[2] < abs_bound)
		upper = x1[2];
	else
	{
		steps = ceil((abs_bound - x0[2]) / dx[2]);
		upper = x0[2] + (steps * dx[2]);
	}	
		
	//Calculate number of grid points per dimension	
	total_ev = rint((upper - lower[2])/dx[2]) + 1;	
																		
	//For each gridpoint in dimensions 3 to N			
	for(j = 0; j < total_ev; j++)
	{							
		//Calculate location of grid point
		gridpoint[2] = lower[2] + (dx[2] * j); 			
		dif_pos[2] = (gridpoint[2] - x0[2])/ dx[2];

		/* This code calculates, a 2D plane formed by the first two dimensions of the space, the optimal
		 * box inside the initial bounding box */

		Z = gridpoint[2] - sample[2];

		//X,Y, along with X2,XY,Y2 form the equation of the 2D rotated plane

		F = Sm1->me[2][2] * Z * Z - 1;

		X  =  2*Sm1->me[0][2]*Z*cosTheta + 2*Sm1->me[1][2]*Z*sinTheta;
		Y  = -2*Sm1->me[0][2]*Z*sinTheta + 2*Sm1->me[1][2]*Z*cosTheta;

		//Calculate displacements and obtain formula (x-xo)^2 / a^2 + % (y-yo)^2/b^2 = 1
		
		termX2 = (X/X2)/2;
		termY2 = (Y/Y2)/2; 
		valor = -F + termX2*termX2*X2 + termY2*termY2*Y2;

		//Calculate new rotated bounding box. UP and RIGHT are the corners of the new bounding box

		upy = sqrt(1/(Y2/valor)) * h;
		rightx = sqrt(1/(X2/valor)) * h;
		
		upx_rot    =  0 * cosTheta + upy * sinTheta;
		upy_rot    = -0 * sinTheta + upy * cosTheta;
		rightx_rot =  rightx * cosTheta + 0 * sinTheta;
		righty_rot = -rightx * sinTheta + 0 * cosTheta;
				
		//Calculate original displacement (rotated ellipse)	
				
		box_center[0] = termX2*cosTheta-termY2*sinTheta;
		box_center[1] = termX2*sinTheta+termY2*cosTheta;
		
		bound[0] = sqrt(upx_rot*upx_rot+rightx_rot*rightx_rot);
		bound[1] = sqrt(upy_rot*upy_rot+righty_rot*righty_rot);
			
		//Calculate lower and upper bound of new BoundingBox	
		for(u = 0; u < 2; u++)
		{
			box_min[u] = (sample[u] - box_center[u]) - bound[u]; 
			box_steps[u] = floor((box_min[u] - x0[u]) / dx[u]);
			lower[u] = (x0[u] > box_min[u])?(x0[u]):(x0[u] + (box_steps[u] * dx[u]));

			box_max[u] = (sample[u] - box_center[u]) + bound[u]; 
			box_steps[u] = ceil((box_max[u] - x0[u]) / dx[u]); 
			box_upper[u] = (x1[u] < box_max[u])?(x1[u]):(x0[u] + (box_steps[u] * dx[u])); 

			tot_ev_per_dim[u] = rint((box_upper[u] - lower[u])/dx[u]);			
		}
	
	    //Calculate the PDF of the defined 2D box
	    compute2DBox_3D(pdf,PC,lower,tot_ev_per_dim,gridpoint,dif_pos,x0,dx,h2,cd,eigenvectors,densValues,densPosition);	
		
	}//End of "per gridpoint" for

  } //End of "per sample" for
    
  #ifdef _OPENMP
  _mm_free(densValues);
  _mm_free(densPosition);
  #endif  

  //Apply k constant to PDF
  #pragma omp for
  for(i=0; i < pdf->total_size; i++)
  	  pdf->PDF[i] = pdf->PDF[i] * k;

  //Calculate integral of PDF  
  #pragma omp for reduction(+:theintegral) 
  for(i=0; i < pdf->total_size; i++)
      theintegral += pdf->PDF[i];

  #pragma omp single
  theintegral = theintegral * total_dx;
  
  //Renormalize PDF using integral
  #pragma omp for
  for(i=0; i < pdf->total_size; i++)
     pdf->PDF[i] = pdf->PDF[i]/theintegral;
  
  //Calculate total volume of renormalized PDF     
  #pragma omp for reduction(+:total_vol)  
  for(i=0; i < pdf->total_size; i++)
	total_vol += pdf->PDF[i];

  }//End of parallel OpenMP Region   
  
  printf("Total integrated PDF: %g. The integral: %f\n",total_vol*total_dx,theintegral);
}

//Compute the PDF of grid spaces of dimension 3 or higher
void computePDFND(MPDFEstimatorPtr mpdf, PDFPtr pdf, MAT *Sm1 , double h , 
		  double detSm1 , double *x0, double *x1, double *dx, 
		  double *bounds, MAT *eigenvectors, int dim)
{
  int i,j,l,u,w; //Loop variables	
  double cd = volumeConstant(dim); //Volume constant
  double k=1./sqrt(detSm1)/mpdf->current/pow(h,mpdf->length); //Constant to recover the volume in the X space from the volume in the PC space  
  double h2=h*h; //Square of bandwith value
  double *PC; // Current sample (PC space)
  double total_vol=0.0;
  double theintegral=0.0;
  double * sample; //Current sample
  double * PCdot;
  
  //Variables to calculate the bounding box of a sample
  double * lower;
  double upper;
  double * gridpoint;
  int * tot_ev_per_dim;
  size_t * dif_pos;
  int total_ev;	
  int steps;
  double abs_bound; //Absolute bound per sample and dimension, given by ellipsoid shape

  //Calculate acumulated volume for the grid space
  double total_dx = 1.0;
  for (i = 0; i < dim; i++)
	total_dx *= dx[i];   

  //Variables to perform the calculation of the 2D layering
  double A,B,C,F,Z,theta,cosTheta,sinTheta,X2,Y2,X,Y,XY,termY2,valor,termX2,upy,rightx,upx_rot,upy_rot,rightx_rot,righty_rot; 
  double bound[2],box_center[2],box_min[2],box_max[2],box_steps[2],box_upper[2];
      			
  //Calculate partial equations for the 2D layering    							
  A = Sm1->me[0][0];
  B = 2 * Sm1->me[0][1];
  C = Sm1->me[1][1];
  theta = atan(B/(A-C))/2;		 
  cosTheta = cos(theta);
  sinTheta = sin(theta);					
  X2 =    Sm1->me[0][0]*cosTheta*cosTheta + 2*Sm1->me[0][1]*cosTheta*sinTheta +   Sm1->me[1][1]*sinTheta*sinTheta;
  XY = -2*Sm1->me[0][0]*cosTheta*sinTheta + 2*Sm1->me[0][1]*cosTheta*cosTheta - 2*Sm1->me[0][1]*sinTheta*sinTheta + 2*Sm1->me[1][1]*cosTheta*sinTheta;
  Y2 =    Sm1->me[0][0]*sinTheta*sinTheta - 2*Sm1->me[0][1]*cosTheta*sinTheta +   Sm1->me[1][1]*cosTheta*cosTheta;		
  
  //Aux vector for OpenMP version
  double * densValues; 
  int * densPosition;
  double * temp_vec;
  double * PCdot_vec;
  double * gridpoint_vec;
 
  //Beginning of OpenMP parallel region								
  #pragma omp parallel default(none)\
  shared(stdout,theintegral,total_vol,total_dx,k,mpdf,pdf,cd,dim,bounds,x0,x1,dx,Sm1,cosTheta,sinTheta,eigenvectors,X2,XY,Y2,h2,h) \
  private(i,j,l,u,w,sample,PC,gridpoint,total_ev,abs_bound,lower,box_upper,tot_ev_per_dim,box_steps,F,X,Y,Z,termX2,termY2,upy,rightx,upx_rot,upy_rot, \
  valor,rightx_rot,righty_rot,bound,box_center,box_min,box_max,PCdot,dif_pos,steps,upper,densValues,densPosition,temp_vec,gridpoint_vec,PCdot_vec)
  {			
								
  //Allocate variables to calculate the bounding box of a sample			
  lower = (double *)malloc(sizeof(double) * dim);
  gridpoint = (double *)malloc(sizeof(double) * dim);
  tot_ev_per_dim = (int *)malloc(sizeof(int) * dim);
  dif_pos = (size_t *)malloc(sizeof(size_t) * dim);			

  #ifdef _OPENMP

  int dim1_max_size = ((ceil(bounds[1] / dx[1]) * 2) + 3);

  densValues = (double *)_mm_malloc(sizeof(double) * dim1_max_size,64);
  densPosition = (int *)_mm_malloc(sizeof(int) * dim1_max_size,64);
  
  temp_vec = (double *)_mm_malloc(sizeof(double) * dim1_max_size,64);
  gridpoint_vec = (double *)_mm_malloc(sizeof(double) * dim1_max_size * dim,64);
  PCdot_vec = (double *)_mm_malloc(sizeof(double) * dim1_max_size * dim,64);
  
  #else
  
  PCdot_vec = (double *)malloc(sizeof(double) * dim);

  #endif

  //Initialize PDF structure to 0s
  #pragma omp for 
  for(i = 0; i < pdf->total_size; i++)
	pdf->PDF[i] = 0.0f;

  //Main calculation loop. For each sample calculate the PDF of its influence area and store in the PDF structure
  #pragma omp for  
  for(i=0;i<mpdf->current;i++) 
  {	  
	sample = MPDFPosition(mpdf,i); //Get current sample
	PC = MPDFPCPosition(mpdf,i); //X is the current sample (scaled as PC)

	//For each sample, calculate its bounding box, 
	//expressed as coordinates of lower corner and number of gridpoints per dimensions
	total_ev = 1;			
	for(j = 2; j < dim; j++)
	{	
		//Lower corner
		abs_bound = sample[j] - bounds[j];
		if (x0[j] > abs_bound)
			lower[j] = x0[j];
		else
		{
			steps = floor((abs_bound - x0[j]) / dx[j]);
			lower[j] = x0[j] + (steps * dx[j]);
		}

		//Upper corner
		abs_bound = sample[j] + bounds[j];
		if (x1[j] < abs_bound)
			upper = x1[j];
		else
		{
			steps = ceil((abs_bound - x0[j]) / dx[j]);
			upper = x0[j] + (steps * dx[j]);
		}	
		
		//Calculate number of grid points per dimension	
		tot_ev_per_dim[j] = rint((upper - lower[j])/dx[j]) + 1;	
		total_ev *= tot_ev_per_dim[j] ;				
	}		
																				
	//For each gridpoint in dimensions 3 to N			
	for(j = 0; j < total_ev; j++)
	{							
		//Calculate location of grid point
		int divisor;
		int eval_point = j;
		for(u = 2; u < dim-1; u++)
		{			
			divisor = 1;
			for(w = u+1; w < dim; w++)
				divisor *= tot_ev_per_dim[w];
			
			gridpoint[u] = lower[u] + (dx[u] * (eval_point / divisor));			
			eval_point = eval_point % divisor;			
		}
		gridpoint[dim-1] = lower[dim-1] + (dx[dim-1] * eval_point); //Last case			
																						    
		//Fill structure with gridpoint position                        
		for(l = 2; l < dim; l++)
			dif_pos[l] = (gridpoint[l] - x0[l])/ dx[l];

		/* This code calculates, a 2D plane formed by the first two dimensions of the space, the optimal
		 * box inside the initial bounding box */

		Z = gridpoint[2] - sample[2];

		//X,Y, along with X2,XY,Y2 form the equation of the 2D rotated plane

		F = Sm1->me[2][2] * Z * Z - 1;

		X  =  2*Sm1->me[0][2]*Z*cosTheta + 2*Sm1->me[1][2]*Z*sinTheta;
		Y  = -2*Sm1->me[0][2]*Z*sinTheta + 2*Sm1->me[1][2]*Z*cosTheta;

		//Calculate displacements and obtain formula (x-xo)^2 / a^2 + % (y-yo)^2/b^2 = 1
		
		termX2 = (X/X2)/2;
		termY2 = (Y/Y2)/2; 
		valor = -F + termX2*termX2*X2 + termY2*termY2*Y2;

		//Calculate new rotated bounding box. UP and RIGHT are the corners of the new bounding box

		upy = sqrt(1/(Y2/valor)) * h;
		rightx = sqrt(1/(X2/valor)) * h;
		
		upx_rot    =  0 * cosTheta + upy * sinTheta;
		upy_rot    = -0 * sinTheta + upy * cosTheta;
		rightx_rot =  rightx * cosTheta + 0 * sinTheta;
		righty_rot = -rightx * sinTheta + 0 * cosTheta;
				
		//Calculate original displacement (rotated ellipse)	
				
		box_center[0] = termX2*cosTheta-termY2*sinTheta;
		box_center[1] = termX2*sinTheta+termY2*cosTheta;
		
		bound[0] = sqrt(upx_rot*upx_rot+rightx_rot*rightx_rot);
		bound[1] = sqrt(upy_rot*upy_rot+righty_rot*righty_rot);
			
		//Calculate lower and upper bound of new BoundingBox	
		for(u = 0; u < 2; u++)
		{
			box_min[u] = (sample[u] - box_center[u]) - bound[u]; 
			box_steps[u] = floor((box_min[u] - x0[u]) / dx[u]);
			lower[u] = (x0[u] > box_min[u])?(x0[u]):(x0[u] + (box_steps[u] * dx[u]));

			box_max[u] = (sample[u] - box_center[u]) + bound[u]; 
			box_steps[u] = ceil((box_max[u] - x0[u]) / dx[u]); 
			box_upper[u] = (x1[u] < box_max[u])?(x1[u]):(x0[u] + (box_steps[u] * dx[u])); 

			tot_ev_per_dim[u] = rint((box_upper[u] - lower[u])/dx[u]);			
		}
	
	    //Calculate the PDF of the defined 2D box	
	    compute2DBox_ND(pdf,PC,lower,tot_ev_per_dim,gridpoint,dif_pos,x0,dx,dim,h2,cd,eigenvectors,densValues,densPosition,PCdot_vec,temp_vec,gridpoint_vec);	
	
	}//End of "per gridpoint" for

  } //End of "per sample" for
    
  //Delete memory structures created by threads   
  free(lower);	
  free(tot_ev_per_dim);				
  free(dif_pos);  
  free(gridpoint);       
  

  #ifdef _OPENMP
  _mm_free(densValues);
  _mm_free(densPosition);
  _mm_free(PCdot_vec);
  _mm_free(temp_vec);
  _mm_free(gridpoint_vec);   
  #else
  free(PCdot_vec);
  #endif  

  //Apply k constant to PDF
  #pragma omp for
  for(i=0; i < pdf->total_size; i++)
  	  pdf->PDF[i] = pdf->PDF[i] * k;

  //Calculate integral of PDF  
  #pragma omp for reduction(+:theintegral) 
  for(i=0; i < pdf->total_size; i++)
      theintegral += pdf->PDF[i];

  #pragma omp single
  theintegral = theintegral * total_dx;
  
  //Renormalize PDF using integral
  #pragma omp for
  for(i=0; i < pdf->total_size; i++)
     pdf->PDF[i] = pdf->PDF[i]/theintegral;
  
  //Calculate total volume of renormalized PDF     
  #pragma omp for reduction(+:total_vol)  
  for(i=0; i < pdf->total_size; i++)
	total_vol += pdf->PDF[i];

  }//End of parallel OpenMP Region   
  
  printf("Total integrated PDF: %g. The integral: %f\n",total_vol*total_dx,theintegral);
}

