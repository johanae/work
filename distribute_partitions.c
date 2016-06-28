/*
Continuation of the serial partitioning program. In partition_serial.c, the data is partitioned, sorted and prepared.
Here, the root node distributes the data to the nodes.

The goal of this program is not really efficiency, but rather just to partition the data in such a way that the resulting format is identical to that
of the parallellised partitioning, on each partition.

Dividing and sending the xadj partitions is potentially a bit confusing. Since local_xadj[0] = 0 always, don't send that index. Then process j
will receive process_j.local_n values.

TODO: Make an assert function and some assert data if at all possible.
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include <mpi.h>
#include <omp.h>

#include "readutility.h"
#include "common.h"

void readAndPrepareData(char *filename, int numpart, TetgenMesh *M, int **partsizes_ptr, int **separatorsizes_ptr);

void writeToFile(TetgenMesh *M, char *filename);

/*
* Goes through the elements array to obtain the indices of required elements.

* Subject to improvements.

* Pending said improvements, idea is currently to make a sort-perm, from which you make a giftlist, deleting duplicates
* as you go.

sortnodes is probably a misleading name.
*/
void sortnodes(int n, int numnodes, int *elements, double *nodes, int *nodelimits, int numpart, int myrank) {

    int i, j, uniques, count, min, max;
    int *perm, *invperm;
    int *wishlist_sizes = (int *)calloc(numpart * sizeof(int), sizeof(int));
    int *giftlist_sizes;
    double *sendnodes, *recvnodes;
    double *nodes_result;

    perm = malloc(4 * n * sizeof(int));
    for (i = 0; i < 4 * n; i ++) perm[i] = i;

    quicksortPerm(elements, perm, 0, 4 * n);

    invperm = malloc(4 * n * sizeof(int));
    for (i = 0; i < 4 * n; i ++) invperm[perm[i]] = i;

    int currentpart = uniques = 0;
    i = -1;
    while (++i < 4 * n) {
        if (i != 0 && (elements[perm[i]] == elements[perm[i - 1]])) continue;

        if (elements[perm[i]] >= nodelimits[currentpart + 1]) currentpart ++;

        uniques ++;
        wishlist_sizes[currentpart] ++;

    }

    int *global2local_elemts = malloc(uniques * sizeof(int));    // local relabelling essentially
    int *local2global_elemts = malloc(uniques * sizeof(int));
    i = -1; count = 0;
    while (++i < 4 * n) {

        if (i != 0 && (elements[perm[i]] == elements[perm[i - 1]])) continue;

        global2local_elemts[count] = elements[perm[i]];
        local2global_elemts[count] = invperm[global2local_elemts[count ++]];

    }

    if (myrank == 0) for (i = 0; i < uniques; i ++) printf("global2local: %d, local2global_elemts: %d.\n", global2local_elemts[i], local2global_elemts[i]);
    // Step 1: Sort the elts, [012345...].
    // Step 2: Give the elts local names.
    // Step 3: Compute the global names.
    // The permutation is bijective on the partition, but not globally.
    // Consider what you have. global2local_elemts[0, 1, ...] gives the SORTED

    giftlist_sizes = (int *)malloc(numpart * sizeof(int));
    MPI_Alltoall(wishlist_sizes, 1, MPI_INT, giftlist_sizes, 1, MPI_INT, MPI_COMM_WORLD);

    if (myrank == 0) for (i = 0; i < nodelimits[1]; i ++) printf("%d\n", global2local_elemts[i]);

    /*// Verify symmetry
    for (i = 0; i < numpart; i ++) {
        printf("I am process %d and %d requests this many of my elements: %d.\n", myrank, i, giftlist_sizes[i]);
        printf("I am process %d and from %d I request this many elements: %d.\n", myrank, i, wishlist_sizes[i]);
    }*/


    // Reminder: nodes currently require their UNSORTED LOCAL names for proper access, but require UNSORTED GLOBAL names for access.
    // SORTED GLOBAL is a nonsense term. Don't use it. The wishlist consists of UNSORTED GLOBAL names.
    //nodes_result = (double *)3 * uniques * sizeof(double);

}

// TODO Cleanup on aisle 3
int main(int narg, char **args) {

    // MPI Variables
    int myrank, numprocs;

    // Root node-only variables
    TetgenMesh *M;
    int numnodes;
    int *partsizes, *separatorsizes;
    int *vertixlimits, *elemlimits;
    int *xadjs, *xadjsizes, *xadjdisps, *adjncy, *adjncysizes, *adjncydisps, *elem, *elemsizes, *elemdisps, *nodesizes, *nodedisps;
    double *nodes;

    // All node variables
    int i, j;
    int recvsize, local_n, local_numnodes, *nodeindices; // nodeindices[j] = value of 1st node on proc j.
    int local_separatorsize, *local_xadj, *local_adjncy, *local_elemts, *limits, *nodelimits;
    double *local_nodes;

    MPI_Init(&narg, &args);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    limits = (int *)malloc((numprocs + 1) * sizeof(int));
    nodelimits = (int *)malloc((numprocs + 1) * sizeof(int));

    if (myrank == 0) {
        char *filename = args[1];
        int numpart = numprocs;
        M = (TetgenMesh *)malloc(sizeof(TetgenMesh));

        readAndPrepareData(filename, numpart, M, &partsizes, &separatorsizes);

        /* For use in MPI_Scatterv, here are all the size arrays and disp arrays required by process 0. Redundancies ahoy! */
        xadjsizes = (int *) malloc(numprocs * sizeof(int));
        xadjdisps = (int *) malloc((numprocs + 1) * sizeof(int));

        xadjsizes[0] = partsizes[0];
        xadjdisps[0] = 0;
        for (i = 1; i < numprocs; i ++) {
            xadjsizes[i] = partsizes[i];
            xadjdisps[i] = xadjdisps[i - 1] + xadjsizes[i - 1];
        }
        xadjdisps[numprocs] = M->n; // so, yeah. I don't remember what this was for. Find this out.

        xadjs = M->xadj;
        for (i = 1; i < numprocs; i ++) {

            int shift = M->xadj[xadjdisps[i]];
            for (j = xadjdisps[i] + 1; j < xadjdisps[i + 1]; j ++) {
                xadjs[j] -= shift;
            }
        }

        adjncysizes = (int *)malloc(numprocs * sizeof(int));
        adjncydisps = (int *)malloc(numprocs * sizeof(int));

        adjncysizes[0] = xadjdisps[partsizes[0]];
        adjncydisps[0] = 0;
        for (i = 1; i < numprocs; i ++) {
            adjncysizes[i] = xadjs[partsizes[i]];
            adjncydisps[i] = adjncydisps[i - 1] + adjncysizes[i - 1];
        }

        elemsizes = (int *)malloc(numprocs * sizeof(int));
        elemdisps = (int *)malloc(numprocs * sizeof(int));

        elemsizes[0] = 4 * partsizes[0];
        elemdisps[0] = 0;
        for (i = 1; i < numprocs; i++) {
            elemsizes[i] = 4 * partsizes[i];
            elemdisps[i] = elemdisps[i - 1] + elemsizes[i - 1];
        }

        numnodes = M->numnodes;

        nodesizes = (int *)malloc(numprocs * sizeof(int));
        nodedisps = (int *)malloc(numprocs * sizeof(int));

        int incr;
        int q = numnodes / numprocs;
        int r = numnodes % numprocs;
        nodesizes[0] = q;
        if (r > 0) nodesizes[0] ++;
        nodedisps[0] = 0;
        nodelimits[0] = 0;
        for (i = 1; i <= numprocs; i ++) {
            incr = q;
            if (i < r) incr ++;
            //nodesizes[i] = 3 * q;
            //if (i < r) nodesizes[i] += 3;
            nodesizes[i] = 3 * incr;
            nodedisps[i] = nodedisps[i - 1] + nodesizes[i - 1];
            printf("%d\n", incr);
            nodelimits[i] = nodelimits[i - 1] + incr;

        }

        adjncy = M->adjncy;
        elem = M->elements;
        nodes = M->nodes;

        limits[0] = 0; nodelimits[0] = 0;
        for (i = 1; i <= numprocs; i ++) {
            limits[i] = limits[i - 1] + partsizes[i - 1];
        }

        // Phew.
    }

    // This can probably be made much faster, but that's not what's important right now.

    MPI_Scatter(partsizes, 1, MPI_INT, &local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(separatorsizes, 1, MPI_INT, &local_separatorsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numnodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(limits, numprocs + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nodelimits, numprocs + 1, MPI_INT, 0, MPI_COMM_WORLD);

    printf("Process %d local_n %d\n", myrank, local_n);

    // No packaging or anything is performed here, I just want to get the sending out of the way.

    local_xadj = (int *)malloc((local_n + 1) * sizeof(int));
    local_xadj[0] = 0;
    recvsize = local_n;

    MPI_Scatterv(&xadjs[1], xadjsizes, xadjdisps, MPI_INT, &local_xadj[1], recvsize, MPI_INT, 0, MPI_COMM_WORLD);

    local_adjncy = (int *)malloc(local_xadj[local_n] * sizeof(int));
    recvsize = local_xadj[local_n];
    MPI_Scatterv(adjncy, adjncysizes, adjncydisps, MPI_INT, local_adjncy, recvsize, MPI_INT, 0, MPI_COMM_WORLD);

    local_elemts = (int *)malloc(4 * local_n * sizeof(int));
    recvsize = 4 * local_n;
    MPI_Scatterv(elem, elemsizes, elemdisps, MPI_INT, local_elemts, recvsize, MPI_INT, 0, MPI_COMM_WORLD);
    printf("Process %d has received local_elemts array of size %d.\n", myrank, recvsize);
    printf("Limits[%d+1]: %d\n", myrank, limits[myrank+1]);
    printf("Nodelimits[%d+1]: %d\n", myrank, nodelimits[myrank+1]);

    local_numnodes = numnodes / numprocs;
    if (myrank < numnodes % numprocs) local_numnodes ++;
    recvsize = 3 * local_numnodes;
    local_nodes = (double *)malloc(3 * local_numnodes * sizeof(double));
    MPI_Scatterv(nodes, nodesizes, nodedisps, MPI_DOUBLE, local_nodes, recvsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    sortnodes(local_n, local_numnodes, local_elemts, nodes, nodelimits, numprocs, myrank);

    if (myrank == 0) {
        deallocateTetgenMesh(M);
        free(partsizes);
        free(separatorsizes);
    }

    MPI_Finalize();
}