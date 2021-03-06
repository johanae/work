// This file is due for a renaming.

typedef struct {

int n, numnodes, numpart;
int *xadj;
int *adjncy;
int *ell;

int *partsizes, *limits, *nodesizes, *nodelimits;

int *separatorsizes;
int *separraysizes;
int *seplimits;
int *separators, *sepparts;

int *elements;
double *nodes;

} TetgenMesh;   // rename to something more politically correct

// Based on the localparttype struct in commondefs.h
typedef struct {

// GLOBAL mesh-related variables
int global_n,global_nodecount;
int *limits;

// LOCAL mesh-related arrays and variables
int n, nodecount, sepsize;
int *sep, *seppart;
int *ell, *elemts;
double *nodes;

// Communication-related variables
int myrank, numprocs;
int neighcount; // number of neighbours & neighbours of neighbours.
int *edges;  // corresponding array.
int *sendidx, *recvidx; // indices of sends/recvs. Partitions are implied by send/recv-disps
int *recvsizes, *recvdisps, *sendsizes, *senddisps;

double *centroids;  // If the time comes.

int ghostcount;   // ~ remoteVcount in commondefs.h. Size of ethereal appendage.


} Nodemeshdata;


void swap(int *a, int *b);
void quicksortPerm(int *arr, int *perm, int beg, int end);
int binarySearch(int key, int first, int last, int initial, int *list);
void deallocateTetgenMesh(TetgenMesh *M);
void deallocateNodemeshdata(Nodemeshdata *d);
void ell2csr(int n, int *ell, int **xadj_ptr, int **adjncy_ptr);
