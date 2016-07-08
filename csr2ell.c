#include <stdio.h>
#include <readutility.h>
#include <commondef.h>

void traverseSeparator(TetgenMesh *M, int myrank, int numprocs, int *separator, int separatorsize) {

    // make an all-to-all call here
    int current, current_neighbour, current_neighbour_squared;

    for (int i = 0; i < separatorsize; i ++) {
        for (int j = M->xadj[i]; j < M->xadj[i + 1]; i ++) {

        }
    }
}

void csr2ell(TetgenMesh *M, int limits, int nodelimits, int myrank, int numprocs, int separatorsize) {


}