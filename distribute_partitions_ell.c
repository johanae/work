/*
Continuation of the serial partitioning program. In partition_serial.c, the data is partitioned, sorted and prepared.
Here, the root node distributes the data to the nodes.

The goal of this program is not really efficiency, but rather just to partition the data in such a way that the resulting format is identical to that
of the parallellised partitioning, on each partition.

There is a risk of the xadj, adjncy arrays losing vital information from the mesh, so in this version, ELL format only is distributed.

Dividing and sending the xadj partitions is potentially a bit confusing. Since local_xadj[0] = 0 always, don't send that index. Then process j
will receive process_j.local_n values.

TODO: Make an assert function and some assert data if at all possible.
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include <mpi.h>
#include <omp.h>

#include "readutility.h"
#include "common.h"

void initialiseTopologies(Nodemeshdata *d, int *trimmedcount_ptr, int **sep_trimmed_ptr, int **seppart_trimmed_ptr, int **orig_trimmed_ptr);
void readAndPrepareData(char *filename, int numpart, TetgenMesh *M);


/*
Writes the partitioned data to binary files, I guess.
Rank 0 writes an init file, containing in order
n
limits[]
separatorsizes[]
The local data have format

*/
void storeDataBinary(char *fname_common, int myrank, int numprocs, int *limits, int *separatorsizes, int local_n, int local_numnodes, int *local_ell, int *local_elmts, double *local_nodes) {

    char fname[200];
    char suffix[20];
    FILE *file;
    int err;

    if (myrank == 0) {

        strcpy(fname, fname_common);
        strcpy(suffix,"_init.bin");
        strcat(fname, suffix);
        file = fopen(fname, "wb");    // possibly stick to nonbinary format for the init file

        err = fwrite(&limits[numprocs], sizeof(int), 1, file);
        err = fwrite(limits, sizeof(int), numprocs + 1, file);
        err = fwrite(separatorsizes, sizeof(int), numprocs, file);

        fclose(file);
    }


    sprintf(suffix, "%04d.bin", myrank);
    strcpy(fname, fname_common);
    strcat(fname, suffix);

    file = fopen(fname, "wb");

    err = fwrite(local_ell, sizeof(int), 4 * local_n, file);
    err = fwrite(local_elmts, sizeof(int), 4 * local_n, file);
    err = fwrite(local_nodes, sizeof(double), 3 * local_numnodes, file);

    fclose(file);

}

/*
I...guess this is just for me to verify that everything is as it should be.
*/
void storeData(char *fname_common, Nodemeshdata *d, int *separatorsizes) {
/*void storeData(char *fname_common, int d->myrank, int d->numprocs, int *limits, int *separatorsizes, int local_n, int local_numnodes, int *d->ell, int *local_elmts, double *local_nodes,
               int separatorsize, int *separator, int *seppart) { */

    char fname[200];
    char suffix[20];
    char buffer[200];
    FILE *file;
    int err;

    if (d->myrank == 0) {

        strcpy(fname, fname_common);
        strcpy(suffix,"_init.txt");
        strcat(fname, suffix);
        file = fopen(fname, "w");    // possibly stick to nonbinary format for the init file

        fprintf(file, "%d\n\n", d->limits[d->numprocs]);
        //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);

        for (int i = 0; i <= d->numprocs; i ++) {

            fprintf(file, "%d %d\n", i, d->limits[i]);
            //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);
        }

        fprintf(file, "\n");
        //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);

        for (int i = 0; i < d->numprocs; i ++) {

            fprintf(file, "%d %d\n", i, separatorsizes[i]);
            //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);
        }

        fclose(file);
    }

    sprintf(suffix, "%04d.txt", d->myrank);
    strcpy(fname, fname_common);
    strcat(fname, suffix);

    file = fopen(fname, "w");

    fprintf(file, "\nNeighbours:\n");

    for (int i = 0; i < d->n; i ++) {

        if (i == d->sepsize) fprintf(file, "----------------------------------------------------------------------------------------\n");
        fprintf(file, "%d %d %d %d %d\n", i, d->ell[4 * i], d->ell[4 * i + 1], d->ell[4 * i + 2], d->ell[4 * i + 3]);
        //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);
    }

    fprintf(file, "\nElements:\n");
    for (int i = 0; i < d->n; i ++) {

        fprintf(file, "%d %d %d %d %d\n", i, d->elemts[4 * i], d->elemts[4 * i + 1], d->elemts[4 * i + 2], d->elemts[4 * i + 3]);
        //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);
    }

    fprintf(file, "\nNodes:\n");
    for (int i = 0; i < d->nodecount; i ++) {

        fprintf(file, "%d %lf %lf %lf\n", i, d->nodes[3 *i ], d->nodes[3 * i + 1], d->nodes[3 * i + 2]);
    }

    fprintf(file, "\nSeparator:\n");
    for (int i = 0; i < d->sepsize; i ++) {
        for (int j = 0; j < 16; j ++) {
            fprintf(file, "%d ", d->sep[16 * i + j]);
        }
        fprintf(file, "   ");
        for (int j = 0; j < 16; j ++) {
            fprintf(file, "%d ", d->seppart[16 * i + j]);
        }
        fprintf(file, "\n");
    }


    fclose(file);

}

/*
* Goes through the elements array to obtain the indices of required elements.

* Subject to improvements.

* Pending said improvements, idea is currently to make a sort-perm, from which you make a giftlist, deleting duplicates
* as you go.

* The sending is done via three allToAlls--first send/recv sizes, then the send/recvlists, then the actual nodes.

sortnodes is probably a misleading name.
*/
void sortnodes(int n, int *numnodes_ptr, int **elements_ptr, double **nodes_ptr, int *nodelimits, int numprocs, int myrank) {

    int i, j, uniques, count, min, max, numnodes;
    int *perm, *uniquesperm, *invmap, *elements, *elements_new;
    int *wishlist, *wishlist_sizes, *wishlist_disps, *giftlist, *giftlist_sizes, *giftlist_disps;
    double *sendnodes, *recvnodes;
    double *nodes, *nodes_result;

    numnodes = *numnodes_ptr;

    elements = *elements_ptr;
    nodes = *nodes_ptr;

    perm = malloc(4 * n * sizeof(int));
    for (i = 0; i < 4 * n; i ++) perm[i] = i;

    quicksortPerm(elements, perm, 0, 4 * n);
    wishlist_sizes = (int *)calloc((numprocs + 1) * sizeof(int), sizeof(int));
    wishlist_disps = (int *)calloc((numprocs  + 1)* sizeof(int), sizeof(int));

    int currentpart = uniques = 0;
    wishlist_disps[0] = 0;
    for (i = 0; i < 4 * n; i ++) {
        if (i != 0 && (elements[perm[i]] == elements[perm[i - 1]])) continue;

        if (elements[perm[i]] >= nodelimits[currentpart + 1]) {
            currentpart ++;
            wishlist_disps[currentpart] = wishlist_disps[currentpart - 1] + wishlist_sizes[currentpart - 1];
        }

        uniques ++;
        wishlist_sizes[currentpart] ++;

    }

    for (i = currentpart + 1; i <= numprocs; i ++) {
        wishlist_sizes[i] = 0;
        wishlist_disps[i] = wishlist_disps[i - 1] + wishlist_sizes[i - 1];
    }


    assert(uniques == wishlist_disps[numprocs]);

    wishlist = (int *)malloc(wishlist_disps[numprocs] * sizeof(int));
    giftlist_sizes = (int *)calloc(numprocs * sizeof(int), sizeof(int));
    giftlist_disps = (int *)calloc((numprocs + 1) * sizeof(int), sizeof(int));

    count = 0;
    uniquesperm = (int *)malloc(uniques * sizeof(int));
    for (i = 0; i < 4 * n; i ++) {

        if (i != 0 && (elements[perm[i]] == elements[perm[i - 1]])) continue;

        uniquesperm[count] = perm[i];
        wishlist[count++] = elements[perm[i]];

    }

    // Send wishlists, receive giftlists.
    MPI_Alltoall(wishlist_sizes, 1, MPI_INT, giftlist_sizes, 1, MPI_INT, MPI_COMM_WORLD);

    giftlist_disps[0] = 0;
    for (i = 1; i <= numprocs; i ++) {
        giftlist_disps[i] = giftlist_disps[i - 1] + giftlist_sizes[i - 1];
    }
    giftlist = (int *)malloc(giftlist_disps[numprocs] * sizeof(int));

    MPI_Alltoallv(wishlist, wishlist_sizes, wishlist_disps, MPI_INT, giftlist, giftlist_sizes, giftlist_disps, MPI_INT, MPI_COMM_WORLD);

    // printf("proc %d sending %d nodes.\n", myrank, giftlist_disps[numprocs] - giftlist_disps[myrank]);
    // printf("proc %d has %d out of %d nodes.\n", myrank, nodelimits[myrank + 1] - nodelimits[myrank], nodelimits[numprocs]);

    sendnodes = (double *)malloc(3 * giftlist_disps[numprocs] * sizeof(double));

    int shift = 3 * nodelimits[myrank];
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

    recvnodes = (double *)malloc(wishlist_disps[numprocs] * sizeof(double));

    MPI_Alltoallv(sendnodes, giftlist_sizes, giftlist_disps, MPI_DOUBLE, recvnodes, wishlist_sizes, wishlist_disps, MPI_DOUBLE, MPI_COMM_WORLD);

    // rename the nodes
    //invmap = (int *)malloc(uniques * sizeof(int));
    float dilation = (float)numnodes / (float) nodelimits[numprocs];
    elements_new = (int *)malloc(n * 4 * sizeof(int));
    for (i = 0; i < 4 * n; i ++) {

        // find elements[i] in wishlist, right and take the index, right? Because the wishlist is already sorted. Also permuted is a misnomer.
        elements_new[i] = binarySearch(elements[i], 0, uniques, (int)(dilation * elements[i]), wishlist);
        // Oh god. I think that's it. Could this be it? Is it finally over?
        // modulo efficiency improvements OBVIOUSLY
        // It's so beautiful.
    }

    free(perm);
    free(uniquesperm);

    free(wishlist);
    free(giftlist);
    free(giftlist_disps);
    free(giftlist_sizes);
    free(wishlist_sizes);
    free(wishlist_disps);

    free(elements);
    free(nodes);
    // Heavens, that felt good.

    *nodes_ptr = recvnodes;
    *elements_ptr = elements_new;
    *numnodes_ptr = uniques;

}

// TODO Cleanup on aisle 3
int main(int narg, char **args) {

    // MPI Variables
    int myrank, numprocs;

    // Root node-only variables
    TetgenMesh *M;
    int *partsizes, *separatorsizes, *separraysize, *seplimits;
    int *ell, *ellsizes, *elldisps, *elem, *elemsizes, *elemdisps, *nodesizes, *nodedisps, *separators, *sepparts;
    double *nodes;

    // All node variables
    int n, numnodes;
    int recvsize, local_n, local_numnodes, *nodeindices; // nodeindices[j] = value of 1st node on proc j.
    int local_separatorsize, *local_ell, *local_elemts, *limits, *nodelimits, *separator, *seppart;
    double *local_nodes;

    MPI_Init(&narg, &args);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if (myrank != 0) limits = (int *)malloc((numprocs + 1) * sizeof(int));
    nodelimits = (int *)malloc((numprocs + 1) * sizeof(int));

    if (myrank == 0) {
        char *filename = args[1];
        M = (TetgenMesh *)malloc(sizeof(TetgenMesh));
        readAndPrepareData(filename, numprocs, M);
        printf("Partitioning and serial rearranging complete.\n");

        n = M->n;
        numnodes = M->numnodes;
        partsizes = M->partsizes;
        limits = M->limits;
        separatorsizes = M->separatorsizes;
        separraysize = M->separraysizes;
        seplimits = M->seplimits;
        separators = M->separators;
        sepparts = M->sepparts;

        ell = M->ell;
        elem = M->elements;
        nodes = M->nodes;

        // For use in MPI_Scatterv, here are all the size arrays and disp arrays required by process 0. Redundancies ahoy!
        assert(limits[numprocs] == M->n);

        elemsizes = (int *)malloc(numprocs * sizeof(int));
        elemdisps = (int *)malloc(numprocs * sizeof(int));

        elemsizes[0] = 4 * partsizes[0];
        elemdisps[0] = 0;
        for (int i = 1; i < numprocs; i++) {
            elemsizes[i] = 4 * partsizes[i];
            elemdisps[i] = elemdisps[i - 1] + elemsizes[i - 1];
        }

        nodesizes = (int *)malloc((numprocs + 1) * sizeof(int));
        nodedisps = (int *)malloc((numprocs + 1) * sizeof(int));

        int incr;
        int q = numnodes / numprocs;
        int r = numnodes % numprocs;

        incr = q; if (r > 0) incr ++;
        nodelimits[0] = 0;
        for (int i = 1; i <= numprocs; i ++) {
            nodelimits[i] = nodelimits[i - 1] + incr;
            incr = q; if (r > i) incr ++;
        }
        assert(nodelimits[numprocs] == numnodes);

        nodedisps[0] = 0;
        for (int i = 0; i < numprocs; i ++) {

            nodesizes[i] = 3 * (nodelimits[i + 1] - nodelimits[i]);
            nodedisps[i] = 3 * nodelimits[i];
        }

        ellsizes = elemsizes;
        elldisps = elemdisps;
        // Phew.

    }

    // This can probably be made much faster, but that's not what's important right now.

    if (myrank == 0) {printf("Starting data transfer.\n"); fflush(stdout);}
    MPI_Scatter(partsizes, 1, MPI_INT, &local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(separatorsizes, 1, MPI_INT, &local_separatorsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numnodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(limits, numprocs + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nodelimits, numprocs + 1, MPI_INT, 0, MPI_COMM_WORLD);    // <-- everything here could probably be packaged together into a single send.
    if (myrank == 0) {printf("Basic data successfully transmitted.\n"); fflush(stdout);}

    local_ell = (int *)malloc(4 * local_n * sizeof(int));
    recvsize = 4 * local_n;
    MPI_Scatterv(ell, ellsizes, elldisps, MPI_INT, local_ell, recvsize, MPI_INT, 0, MPI_COMM_WORLD);
    if (myrank == 0)  {printf("Neighbour data (ELL format) successfully transmitted.\n"); fflush(stdout);}

    local_elemts = (int *)malloc(4 * local_n * sizeof(int));
    recvsize = 4 * local_n;
    MPI_Scatterv(elem, elemsizes, elemdisps, MPI_INT, local_elemts, recvsize, MPI_INT, 0, MPI_COMM_WORLD);
    if (myrank == 0)  {printf("elements data successfully transmitted.\n"); fflush(stdout);}
    //printf("Process %d has received local_elemts array of size %d.\n", myrank, recvsize);
    //printf("Limits[%d+1]: %d\n", myrank, limits[myrank+1]);
    //printf("Nodelimits[%d+1]: %d\n", myrank, nodelimits[myrank+1]);

    local_numnodes = nodelimits[myrank + 1] - nodelimits[myrank];
    recvsize = 3 * local_numnodes;
    local_nodes = (double *)malloc(3 * local_numnodes * sizeof(double));

    MPI_Scatterv(nodes, nodesizes, nodedisps, MPI_DOUBLE, local_nodes, recvsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (myrank == 0)  {printf("nodes data successfully transmitted.\n"); fflush(stdout);}

    separator = (int *)malloc(16 * local_separatorsize * sizeof(int));
    seppart = (int *)malloc(16 * local_separatorsize * sizeof(int));
    recvsize = 16 * local_separatorsize;
    MPI_Scatterv(separators, separraysize, seplimits, MPI_INT, separator, recvsize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(sepparts, separraysize, seplimits, MPI_INT, seppart, recvsize, MPI_INT, 0, MPI_COMM_WORLD);
    if (myrank == 0)  {printf("Separator data successfully transmitted.\n"); fflush(stdout);}

    sortnodes(local_n, &local_numnodes, &local_elemts, &local_nodes, nodelimits, numprocs, myrank);
    if (myrank == 0) {printf("Node reordering and replacement successful. Hooray!\n"); fflush(stdout);}


    Nodemeshdata *local_data = (Nodemeshdata *)malloc(sizeof(Nodemeshdata));
    local_data->myrank = myrank;
    local_data->numprocs = numprocs;
    local_data->global_n = n;
    local_data->global_nodecount = local_numnodes;
    local_data->sepsize = local_separatorsize;
    local_data->limits = limits;
    local_data->n = local_n;
    local_data->nodecount = local_numnodes;
    local_data->ell = local_ell;
    local_data->elemts = local_elemts;
    local_data->nodes = local_nodes;
    local_data->sep = separator;
    local_data->seppart = seppart;

    int trimmedcount, *sep_trimmed, *seppart_trimmed, *ind_trimmed;
    initialiseTopologies(local_data, &trimmedcount, &sep_trimmed, &seppart_trimmed, &ind_trimmed);

    MPI_Barrier(MPI_COMM_WORLD);
    storeData("output_txt/helloworld", local_data, separatorsizes);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0) {
        deallocateTetgenMesh(M);
    }

    free(nodelimits);
    free(limits);

    MPI_Finalize();
}