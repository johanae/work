
#define RNZ 16
#define TRUE 1
#define FALSE 0

#define W 0.002  //Smoothing parameter
#ifdef HASCUDA
	#include <cuda.h>
	#include <cuda_runtime.h>
#endif

#ifdef HASPHI
	#include <scif.h>
	#include <source/COIProcess_source.h>
	#include <source/COIEngine_source.h>
#endif


#define MPI_CHECK(call) \
if((call) != MPI_SUCCESS) { \
printf("MPI error calling \""#call"\"\n"); \
exit(-1); }

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct {

    int n;
    int numNodes;
    
} Nodemesh;

typedef struct {
	int numtet;
	int numstates;
	int* rownz;
	int* I;				//warning, numbers become negative during initialization
	double* D;
	double* A;
	double dt;
	double t;
	double* ODEdata;
	double* ODEcoordinates;
	double* V ;
	double* centroid;
	double* volume;
} ELLmatrixVector;

 typedef struct {
	int* elements;
	int* neighbours;
	double* Gx;
	double* Gy;
	double* Gz;
	
	double* nodes;
	double* centroid;
	double* area;
	double* Nx;
	double* Ny;
	double* Nz;
	double* Sx;
	double* Sy;
	double* Sz;
	
	double* tensor;
	double* volume;
	double totalVolume;
	int numtet;
	int numnodes;
} meshdata;


typedef struct {
	int totalsize;
	int parts;
	int mysize;
	int myrank;
	int sepcount;

	//int t;
	int timesteps;
	int* limits;
	//int* getpart;

//legacy
	int OMPs;
	int GPUs;
	int* GPUstart;
	int* GPUsep;
	int OMPsep;
	int devicecount;
	//int MPIseppartcount;	
//modern
	int* devicelimits;
	int* devicesepsize;
	int* deviceMPIsepsize;
//-------


	//mpi comm variables
	int remoteVcount;			//length of the external part of V
	int sendneighbourcount;		// #neighbours to send to
	int recvneighbourcount;		// #neighbours to receive from
	int* sendneighbours;		// list of neighbours to send to
	int* recvneighbours;		// list of neighbours to receive from
	int* sendcounts;			// #elements to sent to each neighbour 
	int* recvecounts;			// #elements to receive from each neighbour 
	int** sendlists;
	double** sendbuffer;
	
	int numV;
	int* recvstart;
	double* newV;
	double* newODEdata;

} localparttype;

typedef struct {
/*	
 	double*   A __attribute__ ((aligned (32)));
	double*   D __attribute__ ((aligned (32)));
	int*   I __attribute__ ((aligned (32)));
	
	double*   SA __attribute__ ((aligned (32)));
	double*   SD __attribute__ ((aligned (32)));
	int*   SI __attribute__ ((aligned (32)));

	double*   MSA __attribute__ ((aligned (32)));
	double*   MSD __attribute__ ((aligned (32)));
	int*   MSI __attribute__ ((aligned (32)));
*/

	double*   A ;
	double*   D ;
	int*   I ;

double*   V ;
	
	double*   SA ;
	double*   SD ;
	int*   SI;

	double*   MSA ;
	double*   MSD ;
	int*   MSI ;

	int msize;	
	int mstart;

	int ssize;	
	int sstart;

	int mssize;	
	int msstart;

	int mainrest;
	int mainblocks;
	int seprest;
	
	int mainsize;
	int separatorsize;
	int** separatormask;
	
} localOMPtype;


typedef struct {
	int my_tid;
	int size;	
	int gstart;
	int* perm;

	int my_commid;
	int sendsize;
	int sendstart;
	int recvsize;
	int recvstart;
	
} localthreadtype;


typedef struct {
	
 	double*  A;
	double*   D;
	double*   V;
	int*   I;
	double*   newV;
	size_t pitch;
	
	int cudanum;

	int mainrest;
	int mainblocks;
	int seprest;
	
	int size;
	int sepsize;
	int** separatormask;
	
	int mainremainderthreads;
	int mainblocksPerGrid;

	int sepremainderthreads;
	int sepblocksPerGrid;

#ifdef HASCUDA
	cudaStream_t mainstream;
	cudaStream_t sepstream;
	cudaStream_t cpystream;
#endif
	
} localCUDAtype;
