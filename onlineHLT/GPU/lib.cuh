extern "C"
{


void allocateArray(void **devPtr, int size);
void freeArray(void *devPtr);

void threadSync();

void copyArrayFromDevice(void* host, const void* device,  int size);
void copyArrayToDevice(void* device, const void* host, int offset, int size);

int iDivUp(int a, int b);

// compute grid and thread block size for a given number of elements
void computeGridSize(int n, int blockSize, int &numBlocks, int &numThreads);


void Gtest(float* position, int*hash, int*map, int* manager, float *out, int* number,int cores);
// void Gcalc_hash(int number,float *fpara,float *fhit,int * volumeC,int * rowC);
// void Gclear(int * dataC,int number);
// void GSclass_hash(int number,float *fpara,float *fhit, int * volume,int * row);
// void Gtracking(int number,float *fpara,float *fhit,int * ivolumeC,int * irowC,int * itrackC,float *ftrack);




}
