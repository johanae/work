#include <stdio.h>
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

int binarySearch(int key, int first, int last, int initial, int *list) {

    int middle = initial;

    while (first <= last) {

        if (list[middle] < key) {
            first = middle + 1;
        } else if (list[middle] > key) {
            last = middle - 1;
        } else {
            return middle;
        }

        middle = (first + last)/2;
    }

    if (first > last) {
        printf("Something has gone like catastrophically wrong.\n");
        return -1;
    }
}

void deallocateTetgenMesh(TetgenMesh *M) {

    free(M->xadj);
    free(M->adjncy);
    free(M->ell);
    free(M->elements);
    free(M->nodes);

    free(M->partsizes);
    free(M->nodesizes);
    free(M->separatorsizes);
    free(M->separraysizes);
    free(M->seplimits);
    free(M->separators);
    free(M->sepparts);
    free(M);
}


void deallocateNodemeshdata(Nodemeshdata *d) {

    free(d->limits);
    free(d->sep);
    free(d->seppart);
    free(d->ell);
    free(d->elemts);
    free(d->nodes);

    if (d->edges) free(d->edges);

    if (d->recvsizes) free(d->recvsizes);
    if (d->recvdisps) free(d->recvdisps);
    if (d->sendsizes) free(d->sendsizes);
    if (d->senddisps) free(d->senddisps);
}

void ell2csr(int n, int *ell, int **xadj_ptr, int **adjncy_ptr) {

    int *xadj, *adjncy;

    *xadj_ptr = xadj = malloc((n + 1) * sizeof(int));
    xadj[0] = 0;

    int sum = 0;
    for (int i = 0; i < n; i ++) {
        for (int j = 0; j < 4; j ++) {
            if (ell[4 * i + j] != -1) sum ++;
        }
        xadj[i + 1] = sum;
    }

    *adjncy_ptr = adjncy = malloc(sum * sizeof(int));

    int k = 0;
    for (int i = 0; i < n; i ++) {
        for (int j = 0; j < 4; j ++) {
            if (ell[4 * i + j] != -1) adjncy[k++] = ell[4 * i + j];
        }
    }
}