// 
//	kernel.cu - Kernels for fractal (Julia and Mandelbrot) set generation
//

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <GL/freeglut.h>
#include <stdio.h>
#include <stdlib.h> // including: #define EXIT_SUCCESS    0  #define EXIT_FAILURE    1

struct cuComplex 
{
    GLdouble   r;
    GLdouble   i;
    __device__ cuComplex( GLdouble a, GLdouble b ) : r(a) , i( b)  { }
    __device__ GLdouble magnitude2( void ) 
    {
        return r * r + i * i;
    }
    __device__ cuComplex operator*(const cuComplex& a) 
    {
        return cuComplex( r*a.r - i*a.i, i*a.r + r*a.i) ;
    }
    __device__ cuComplex operator+(const cuComplex& a) 
    {
        return cuComplex( r+a.r, i+a.i) ;
    }
};

__device__ float toColor (int i)
{
	float intensity [10] = {0, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f};

	return intensity[i%10];
}

//=================================================================================================
//
//	Mandelbrot and Julia set GPU device functions
//
//=================================================================================================

__device__ int mandelbrot(int row, int col, int width, int height, GLdouble Rmin, GLdouble Rmax, GLdouble Imin, GLdouble Imax, int nIterations)
{
	// Normalize (row, col) to {(R,I) | Rmin < R < Rmax, Imin < I < Imax }
	GLdouble R = ((Rmax - Rmin)/(float)width) * (float)col + Rmin;
    //float I = ((Imax - Imin)/(float)height) * (float)row + Imin;
	GLdouble I = Imax - ((Imax - Imin)/(float)height) * (float)row;
	//if (I < 0.00001f && I > -0.00001f) I = 0.0f;

    cuComplex c(R, I);
    cuComplex a(R, I) ;
    int i = 0;
    for (i = 0; i < nIterations; i++) 
    {
        a = a * a + c;
        if (a. magnitude2() > 4.0)
            return i;
    }

    return i;
}

/////////////////////////////////////////////////////////////////////////////////
//
//	Kernel : kernel1
//
//		Computes colors of width by height pixels representing Julia set in
//		{(u,v) | -scale < u < scale, -scale < v < scale }
//
//	Thread grid requirments:
//
//		1) 2D grid of 2D thread blocks covering width by height pixels
//		2) one pixel per thread computing
//
/////////////////////////////////////////////////////////////////////////////////

__global__ void Mandelbrot_kernel(int *ptr,int width, int height, GLdouble Rmin, GLdouble Rmax, GLdouble Imin, GLdouble Imax, int nIterations) 
{
    // map from threadIdx/BlockIdx to pixel position
    int col = threadIdx.x + blockIdx.x * blockDim.x;	// column index to the width X height pixels
    int row = threadIdx.y + blockIdx.y * blockDim.y;	// row index to the width X height pixels

    // Assuming the origin of the width X height pixels is at upper-left corner
    if (row < height && col < width) 
    {
		// Calculate Mandelbrot value at (x,hy) position
		int index = (col + (height-row-1) * width);
		ptr[index] = mandelbrot(row, col, width, height, Rmin, Rmax, Imin, Imax, nIterations);
/*
		if (*ptr[index] == nIterations)
                {}
            else {                 
                       if(cnt>=0&&cnt<=31)   {b=cnt*4; g=cnt*8; r=0;  }
                  else if(cnt>=32&&cnt<=63)  {b=200; g=500-cnt*8; r=0;  }                      
		  else if(cnt>=64&&cnt<=95)  {b=200; g=0; r=(cnt-64)*4;}
		  else if(cnt>=96&&cnt<=127) {r=200; g=0; b=1000-cnt*8;}
		  else if(cnt>=128&&cnt<=159){r=200; g=(cnt-128)*8; b=0;}
		  else if(cnt>=160&&cnt<=191){g=200; r=1500-cnt*8; b=0;}
		  else if(cnt>=192&&cnt<=223){g=200; r=0; b=(cnt-192)*8;}
          else if(cnt>=224&&cnt<=255){g=230; r=(cnt-224)*8; b=256;}


                 //to change color by prssing key 'c'
                   	
                                tr=r;tb=b;tg=g;
                                switch(c)
				{
				case 0: break;
				case 1: r=tb;b=tr;break;
				case 2: r=tg;g=tr;break;
				case 3: b=tg;g=tb;break;
				case 4: r=tg; g=tb; b=tr; break;
				case 5: r=tb; g=tr; b=tg; break;
                      
				}	
			
			   
                    
                   glColor3f(r/256,g/256,b/256);
                   
                glVertex3d(col - nx / 2, row - ny / 2, 0.0f);
*/
	}
}


 #define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPU assert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//=======================================================================================================
//
//	Compute Julia set using CUDA
//
//=======================================================================================================

#define block_size (16)

//=======================================================================================================
//
//	Compute Mandelbrot set using CUDA
//
//=======================================================================================================

extern "C" int cuComputeMandelbrotSet (int *ptr,int width, int height, GLdouble Rmin, GLdouble Rmax, GLdouble Imin, GLdouble Imax, int nIterations)
{
	
	printf("@cuComputeMandelbrotSet %d == %d\n",width, height);
	int *d_ptr = 0;
    
	cudaError_t cudaStatus;

	cudaDeviceReset();

	// Make sure CUDA device 0 is available
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) 
	{
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        return 1;
    }
 
    // Allocate device memory
    if (cudaMalloc((void **)&d_ptr, width * height * sizeof(int)) != cudaSuccess)
    {
		printf("cuda mem failed");
		
        fprintf(stderr, "!!!! device memory allocation error (allocate A)\n");
		return EXIT_FAILURE;
    }
   
	cudaDeviceSynchronize();
	// Copy host memory to device
	cudaStatus = cudaMemcpy (d_ptr, ptr, width*height*sizeof(int), cudaMemcpyHostToDevice);
	 if (cudaStatus != cudaSuccess)
    {
        printf("cudaMemcpy (d_ptr, ptr) returned error code\n", cudaStatus);
        exit(EXIT_FAILURE);
    }

	// Setup execution parameters and call kernel
    dim3 block(block_size, block_size);
    dim3 grid ((width+block_size-1)/block_size, (height+ block_size-1)/block_size);
	Mandelbrot_kernel<<< grid, block >>>(d_ptr, width, height, Rmin, Rmax, Imin, Imax, nIterations); 
	//gpuErrchk( cudaPeekAtLastError() ); 
	//gpuErrchk( cudaDeviceSynchronize() );

	// Copy result from device to host
    cudaStatus = cudaMemcpy(ptr, d_ptr, width*height*sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess)
    {
        printf("cudaMemcpy (ptr, d_ptr) returned error code %d\n", cudaStatus);
        exit(EXIT_FAILURE);
    }

	
    // Device memory clean up
    if (cudaFree(d_ptr) != cudaSuccess)
    {
        fprintf(stderr, "!!!! memory free error (d_ptr)\n");
		return EXIT_FAILURE;
    }
	return 0;
}
/*
int MandelbrotSet_GPU_GL (float *d_dst, int width, int height, float Rmin, float Rmax, float Imin, float Imax, int nIterations)
{
	// Setup execution parameters and call kernel
    dim3 block(block_size, block_size);
    dim3 grid ((width+block_size-1)/block_size, (height+ block_size-1)/block_size);
	Mandelbrot_kernel<<< grid, block >>>(d_dst, width, height, Rmin, Rmax, Imin, Imax, nIterations); 

	return 0;
}
*/

extern "C" bool resetCUDADevice()
{
	cudaError_t cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) 
	{
        fprintf(stderr, "cudaDeviceReset failed!");
        return false;
    }

	return true;
}