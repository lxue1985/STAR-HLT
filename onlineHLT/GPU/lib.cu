#include <cstdlib>
#include <cstdio>
#include <string.h>
#include "lib_kernel.cuh"


#  define CUDA_SAFE_CALL_NO_SYNC( call) do {                                 \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#  define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);                                            \


#include <cuda_gl_interop.h>

extern "C"
{



void allocateArray(void **devPtr, size_t size)
{
    CUDA_SAFE_CALL(cudaMalloc(devPtr, size));
}

void freeArray(void *devPtr)
{
    CUDA_SAFE_CALL(cudaFree(devPtr));
}

void threadSync()
{
    CUDA_SAFE_CALL(cudaThreadSynchronize());
}

void copyArrayFromDevice(void* host, const void* device,  int size)
{   

    CUDA_SAFE_CALL(cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost));

}

void copyArrayToDevice(void* device, const void* host, int offset, int size)
{
    CUDA_SAFE_CALL(cudaMemcpy((char *) device + offset, host, size, cudaMemcpyHostToDevice));
}


int iDivUp(int a, int b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

// compute grid and thread block size for a given number of elements
void computeGridSize(int n, int blockSize, int &numBlocks, int &numThreads)
{
    numThreads = min(blockSize, n);
    numBlocks = iDivUp(n, numThreads);
}


void Gtest(float* position, int*hash,int*map,int* manager, float *out, int *number,int cores){

    // int numThreads, numBlocks;
    // computeGridSize(number, 16, numBlocks, numThreads);
	
	//cout<<numBlocks<<"  "<<numThreads<<endl;
	// printf("\n numBlocks %d  numThreads%d number %d \n",numBlocks,numThreads,number);
	// CUDA_SAFE_CALL(cudaBindTexture(0, positionTex, position, 10000*sizeof(float4)));
    // CUDA_SAFE_CALL(cudaBindTexture(0, hashTex, hash, HASHSIZE*sizeof(int2)));

test<<<cores,24>>>((float4*)position,(int2*)hash,map,manager,out,number);

    // CUDA_SAFE_CALL(cudaUnbindTexture(positionTex));
    // CUDA_SAFE_CALL(cudaUnbindTexture(hashTex));
}

// void Gcalc_hash(int number,float *fpara,float *fhit,int * volumeC,int * rowC){
    // int numThreads, numBlocks;
    // computeGridSize(number, 256, numBlocks, numThreads);
	// calc_hash<<< numBlocks, numThreads >>>(number,fpara,fhit,volumeC,rowC);
	

// }


// void Gclear(int * dataC,int number){
    // int numThreads, numBlocks;
    // computeGridSize(number, 256, numBlocks, numThreads);
	// clear<<< numBlocks, numThreads >>>(dataC,number);

// }

// void GSclass_hash(int number,float *fpara,float *fhit,int * volumeC,int * rowC){

	// class_hash<<<1,1>>>(number,fpara,fhit,volumeC,rowC);

// }

// void Gtracking(int number,float *fpara,float *fhit,int * ivolumeC,int * irowC,int * itrackC,float *ftrack){
    // int numThreads, numBlocks;
    // computeGridSize(number, 16, numBlocks, numThreads);
	
	// tracking<<<numBlocks, numThreads>>>( number,fpara,fhit, ivolumeC, irowC, itrackC, ftrack);
// }





}