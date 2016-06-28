typedef struct {

int n, numnodes;
int *xadj;
int *adjncy;

int *elements;
double *nodes;

} TetgenMesh;



void swap(int *a, int *b);
void quicksortPerm(int *arr, int *perm, int beg, int end);
void deallocateTetgenMesh(TetgenMesh *M);
//int sddMPI_getpart(int i,localparttype* L);
