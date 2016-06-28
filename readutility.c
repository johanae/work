#include <stdlib.h>
#include "readutility.h"

void swap(int *a, int *b) {

    int t = *a;
    *a = *b;
    *b = t;
}

/*
Utility thing. Consider switching sorting algorithm at some depth, something must obviously be done with the pivot element...
*/
void quicksortPerm(int *arr, int *perm, int beg, int end) {
    int l, r;
    if (end > beg + 1) {
        int piv = arr[perm[beg]], l = beg + 1, r = end;
    while (l < r) {
        if (arr[perm[l]] <= piv) {
            l++;
        } else {
            swap(&perm[l], &perm[--r]);
        }
    }
    swap(&perm[--l], &perm[beg]);
    quicksortPerm(arr, perm, beg, l);
    quicksortPerm(arr, perm, r, end);
  }
}

void deallocateTetgenMesh(TetgenMesh *M) {

    free(M->xadj);
    free(M->adjncy);
    free(M->elements);
    free(M->nodes);
    free(M);
}

/*
Identify process containing node #i.
*/
/*int sddMPI_getpart(int i, localparttype* L) {
    int k=0;
    while(i >= L->limits[k+1])
        k++;
    assert(k <= L->parts);
    return(k);
}*/