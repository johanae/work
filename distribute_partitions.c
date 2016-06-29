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

void readAndPrepareData(char *filename, int numprocs, TetgenMesh *M, int **partsizes_ptr, int **separatorsizes_ptr);

void writeToFile(TetgenMesh *M, char *filename);

/*
* Goes through the elements array to obtain the indices of required elements.

* Subject to improvements.

* Pending said improvements, idea is currently to make a sort-perm, from which you make a giftlist, deleting duplicates
* as you go.

* The sending is done via three allToAlls--first send/recv sizes, then the send/recvlists, then the actual nodes.

sortnodes is probably a misleading name.
*/
void sortnodes(int n, int numnodes, int *elements, double **nodes_ptr, int *nodelimits, int numprocs, int myrank) {

    int i, j, uniques, count, min, max;
    int *perm, *invperm;
    int *wishlist, *wishlist_sizes, *wishlist_disps, *giftlist, *giftlist_sizes, *giftlist_disps;
    double *sendnodes, *recvnodes;
    double *nodes, *nodes_result;

    nodes = *nodes_ptr;

    perm = malloc(4 * n * sizeof(int));
    for (i = 0; i < 4 * n; i ++) perm[i] = i;

    quicksortPerm(elements, perm, 0, 4 * n);


    wishlist_sizes = (int *)calloc((numprocs + 1) * sizeof(int), sizeof(int));
    wishlist_disps = (int *)malloc((numprocs  + 1)* sizeof(int));

    int currentpart = uniques = 0;
    wishlist_disps[0] = 0;
    for (i = 0; i < 4 * n; i ++) {
        if (i != 0 && (elements[perm[i]] == elements[perm[i - 1]])) continue;

        if (elements[perm[i]] >= nodelimits[currentpart + 1]) {
            wishlist_disps[++currentpart] = wishlist_disps[currentpart - 1] + wishlist_sizes[currentpart - 1];
        }

        uniques ++;
        wishlist_sizes[currentpart] ++;

    }


    //printf("%d\n", wishlist_sizes[myrank]);

    for (i = currentpart + 1; i <= numprocs; i ++) {
        wishlist_sizes[i] = 0;
        wishlist_disps[i] = wishlist_disps[i - 1] + wishlist_sizes[i - 1];
    }

    wishlist = (int *)malloc(wishlist_disps[numprocs] * sizeof(int));
    giftlist_sizes = (int *)malloc(numprocs * sizeof(int));
    giftlist_disps = (int *)malloc((numprocs + 1) * sizeof(int));

    count = 0;
    for (i = 0; i < 4 * n; i ++) {

        if (i != 0 && (elements[perm[i]] == elements[perm[i - 1]])) continue;

        wishlist[count++] = elements[perm[i]];

    }

    // Send wishlists, receive giftlists. Keep in mind, then: for proc i wishlist[j] = giftlist[i] on proc j.
    MPI_Alltoall(wishlist_sizes, 1, MPI_INT, giftlist_sizes, 1, MPI_INT, MPI_COMM_WORLD);

    giftlist_disps[0] = 0;
    for (i = 1; i <= numprocs; i ++) {
        giftlist_disps[i] = giftlist_disps[i - 1] + giftlist_sizes[i - 1];
    }

    giftlist = (int *)malloc(giftlist_disps[numprocs] * sizeof(int));

    MPI_Alltoallv(wishlist, wishlist_sizes, wishlist_disps, MPI_INT, giftlist, giftlist_sizes, giftlist_disps, MPI_INT, MPI_COMM_WORLD);

    sendnodes = (double *)malloc(3 * giftlist_disps[numprocs] * sizeof(double));

    printf("proc %d sending %d nodes.\n", myrank, giftlist_disps[numprocs] - giftlist_disps[myrank]);
    printf("proc %d has %d out of %d nodes.\n", myrank, nodelimits[myrank + 1] - nodelimits[myrank], nodelimits[numprocs]);

    // retrieve nodes
    sendnodes = (double *)malloc(3 * giftlist_disps[numprocs] * sizeof(double));
    int shift =3 * nodelimits[myrank];
    printf("Proc %d shift: %d\n", myrank, shift);
    for (i = 0; i < numprocs; i ++) {

        for (j = giftlist_disps[i]; j < giftlist_disps[i + 1]; j ++) {

            sendnodes[3 * j] = nodes[3 * giftlist[j] - shift];
            sendnodes[3 * j + 1] = nodes[3 * giftlist[j]  - shift + 1];
            sendnodes[3 * j + 2] = nodes[3 * giftlist[j]  - shift + 2];
        }
    }

    for (i = 0; i <= numprocs; i ++) {
        giftlist_sizes[i] *= 3;
        giftlist_disps[i] *= 3;
        wishlist_sizes[i] *= 3;
        wishlist_disps[i] *= 3;
    }

    // rename the nodes
    for (i = 0; i < 4 * n; i ++) {

        //elements[i] = invperm[elements[i]]; // I mean, something like this?
    }

    recvnodes = (double *)malloc(wishlist_disps[numprocs] * sizeof(double));
    MPI_Alltoallv(sendnodes, giftlist_sizes, giftlist_disps, MPI_DOUBLE, recvnodes, wishlist_sizes, wishlist_disps, MPI_DOUBLE, MPI_COMM_WORLD);

    // TODO: rename nodes

    free(giftlist_disps);
    free(giftlist_sizes);
    free(wishlist_sizes);
    free(wishlist_disps);
    // Heavens, that felt good.

    *nodes_ptr = recvnodes;

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
        M = (TetgenMesh *)malloc(sizeof(TetgenMesh));

        readAndPrepareData(filename, numprocs, M, &partsizes, &separatorsizes);

        // For use in MPI_Scatterv, here are all the size arrays and disp arrays required by process 0. Redundancies ahoy!
        xadjsizes = (int *) malloc(numprocs * sizeof(int));
        xadjdisps = (int *) malloc((numprocs + 1) * sizeof(int));

        xadjsizes[0] = partsizes[0];
        xadjdisps[0] = 0;
        for (i = 1; i < numprocs; i ++) {
            xadjsizes[i] = partsizes[i];
            xadjdisps[i] = xadjdisps[i - 1] + xadjsizes[i - 1];
        }
        xadjdisps[numprocs] = M->n;

        xadjs = malloc((M->n + 1) * sizeof(int));
        memcpy(xadjs, M->xadj, (M->n + 1) * sizeof(int));
        //xadjs = M->xadj;
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

        nodesizes = (int *)malloc((numprocs + 1) * sizeof(int));
        nodedisps = (int *)malloc((numprocs + 1) * sizeof(int));

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

            nodesizes[i] = 3 * incr;
            nodedisps[i] = nodedisps[i - 1] + nodesizes[i - 1];
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
    if (myrank == 0) {printf("Basic data successfully transmitted.\n"); fflush(stdout);}

    // No packaging or anything is performed here, I just want to get the sending out of the way.
    local_xadj = (int *)malloc((local_n + 1) * sizeof(int));
    local_xadj[0] = 0;
    recvsize = local_n;

    MPI_Scatterv(&xadjs[1], xadjsizes, xadjdisps, MPI_INT, &local_xadj[1], recvsize, MPI_INT, 0, MPI_COMM_WORLD);
    if (myrank == 0) {printf("xadj data successfully transmitted.\n"); fflush(stdout);}

    local_adjncy = (int *)malloc(local_xadj[local_n] * sizeof(int));
    recvsize = local_xadj[local_n];
    MPI_Scatterv(adjncy, adjncysizes, adjncydisps, MPI_INT, local_adjncy, recvsize, MPI_INT, 0, MPI_COMM_WORLD);
    if (myrank == 0)  {printf("adjncy data successfully transmitted.\n"); fflush(stdout);}

    local_elemts = (int *)malloc(4 * local_n * sizeof(int));
    recvsize = 4 * local_n;
    MPI_Scatterv(elem, elemsizes, elemdisps, MPI_INT, local_elemts, recvsize, MPI_INT, 0, MPI_COMM_WORLD);
    if (myrank == 0)  {printf("elements data successfully transmitted.\n"); fflush(stdout);}
    //printf("Process %d has received local_elemts array of size %d.\n", myrank, recvsize);
    //printf("Limits[%d+1]: %d\n", myrank, limits[myrank+1]);
    //printf("Nodelimits[%d+1]: %d\n", myrank, nodelimits[myrank+1]);

    local_numnodes = numnodes / numprocs;
    if (myrank < numnodes % numprocs) local_numnodes ++;
    recvsize = 3 * local_numnodes;
    local_nodes = (double *)malloc(3 * local_numnodes * sizeof(double));
    MPI_Scatterv(nodes, nodesizes, nodedisps, MPI_DOUBLE, local_nodes, recvsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (myrank == 0)  {printf("nodes data successfully transmitted.\n"); fflush(stdout);}

    sortnodes(local_n, local_numnodes, local_elemts, &local_nodes, nodelimits, numprocs, myrank);
    if (myrank == 0) {printf("All nodes data successfully rearranged, transmitted.\n"); fflush(stdout);}

    if (myrank == 0) {
        free(xadjsizes);
        free(xadjdisps);
        free(adjncysizes);
        free(adjncydisps);
        free(elemsizes);
        free(elemdisps);
        free(nodesizes);
        free(nodedisps);
        free(partsizes);
        free(separatorsizes);

        deallocateTetgenMesh(M);
    }

    free(limits);
    free(nodelimits);

    MPI_Finalize();
}