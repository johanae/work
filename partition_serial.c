/*
* This file contains functions for processing and partitioning TETGEN data.
* That means that the .neigh, .ele and .node files are read to a single process.
* The single-process then partitions the vertices via the .neigh data, and then sorts the resulting .neigh and .ele data accordingly.
*
* Currently, the .node data is NOT sorted. Rather, it is partitioned evenly as-is to be sent to the various processes.
*
* The end goal is that upon successful distribution of data among the nodes, the distributed formats will be the same as if parMETIS and parallell initial
* data distribution was used, such that from this point, the same program can be used.
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include <omp.h>

#include "readutility.h"
#include "common.h"
#include "kaHIP_interface.h"
#include "metis.h"

/*
* This function is mostly for testing purposes--it merely reads part of the data rather than all of it
.*/
void readneigh(int *n, int **xadj_ptr, int **adjncy_ptr, char *infile) {

    int i, j, *xadj, *adjncy, *nodes;
    double *nodeVals;
    int edgesLocal;
    char fname[200], suffix[20];
    FILE *fpdata;

    int dump, read;

    int err;

    /* .neigh file*/
    strcpy(fname, infile);
    strcpy(suffix, ".neigh");
    strcat(fname, suffix);
    fpdata = fopen(fname, "r");
    if (fpdata == NULL) {
        printf("\nFailure to open input file: %s.\n", fname);
        exit(0);
    }

    err = fscanf(fpdata, "%d", n);


    *xadj_ptr = xadj = (int *)malloc((*n + 1) * sizeof(int));
    xadj[0] = 0;
    *adjncy_ptr = adjncy = (int *)malloc(4 * (*n) * sizeof(int));

    for(int j = 0; j < 1; ++j) {
        err = fscanf(fpdata, "%d", &dump);
    }

    for (int i = 0; i < *n; i++) {

        err = fscanf(fpdata, "%d", &dump);
        edgesLocal = 4;

        for (int j = 0; j < 4; j++) {
            err = fscanf(fpdata, "%d", &read);

            if (read != -1) {
                adjncy[xadj[i] + j] = read;
            } else {
                edgesLocal --;
            }
        }
        xadj[i + 1] = xadj[i] + edgesLocal;
    }

    fclose(fpdata);
}

/* Reads a .neigh input file to a CSR type format. This is a temporary function, subject to expansion and renaming. */
/* Arguments: (n, xadj, adjncy), filename. */
/* Input: .neigh file, obtained via Tetgen. */
/* Output: integer n, xadj, adjncy arrays.*/
void readmesh(TetgenMesh *M, char *infile) {
/*void readmesh(int *n, int **xadj_ptr, int **adjncy_ptr, int *numnodes, int **nodes_ptr, double **nodeVals_ptr, char *infile) {*/

    int i, j, n, numnodes;
    int *xadj, *adjncy, *elements;
    double *nodes;
    int edgesLocal;
    char fname[200], suffix[20];
    FILE *fpdata;

    // double in;
    int dump, read;

    int err;

    /* .neigh file*/
    strcpy(fname, infile);
    strcpy(suffix, ".neigh");
    strcat(fname, suffix);
    fpdata = fopen(fname, "r");
    if (fpdata == NULL) {
        printf("\nFailure to open input file: %s.\n", fname);
        exit(0);
    }

    err = fscanf(fpdata, "%d", &n);

    M->xadj = xadj = (int *)malloc((n + 1) * sizeof(int));
    M->n = n;
    M->adjncy = adjncy = (int *)malloc(4 * n * sizeof(int));

    xadj[0] = 0;
    for(int j = 0; j < 1; ++j) {
        err = fscanf(fpdata, "%d", &dump);
    }

    for (int i = 0; i < n; i++) {

        err = fscanf(fpdata, "%d", &dump);
        edgesLocal = 4;

        for (int j = 0; j < 4; j++) {
            err = fscanf(fpdata, "%d", &read);

            if (read != -1) {
                adjncy[xadj[i] + j] = read;
            } else {
                edgesLocal --;
            }
        }
        xadj[i + 1] = xadj[i] + edgesLocal;
    }

    fclose(fpdata);

    for (i = xadj[n]; i < 4 * n; i ++) adjncy[i] = -1;


    /* .ele file */
    strcpy(fname, infile);
    strcpy(suffix, ".ele");
    strcat(fname, suffix);
    fpdata = fopen(fname, "r");
    if (fpdata == NULL) {
        printf("\nFailure to open input file.\n");
        exit(0);
    }

    err = fscanf(fpdata, "%d", &dump);
    assert(n == dump);

    M->elements = elements = (int *)malloc(4 * n * sizeof(int));



    for (j = 0; j < 2; j ++) {
        err = fscanf(fpdata, "%d", &dump);
    }

    for (i = 0; i < n; i ++) {

        err = fscanf(fpdata, "%d", &dump);
        for (j = 0; j < 4; j ++) {
            fscanf(fpdata, "%d", &elements[4*i + j]);
        }
    }

    fclose(fpdata);


    /* .node file */
    strcpy(fname, infile);
    strcpy(suffix, ".node");
    strcat(fname, suffix);
    fpdata = fopen(fname, "r");
    if (fpdata == NULL) {
        printf("\nFailure to open input file.\n");
        exit(0);
    }

    err = fscanf(fpdata, "%d", &numnodes);

    M->numnodes = numnodes;
    M->nodes = nodes = (double *)malloc(3 * numnodes * sizeof(double));

    for (j = 0; j < 3; j ++) {
        err = fscanf(fpdata, "%d", &dump);
    }

    for (i = 0; i < numnodes; i ++) {

        err = fscanf(fpdata, "%d", &dump);
        for (j = 0; j < 3; j ++) {
            err = fscanf(fpdata, "%lf", &nodes[3 * i + j]);
        }
    }

    fclose(fpdata);

}

/*
Calls serial METIS partitioner.
ARGUMENTS:
(n, xadj, adjncy) CSR-type neighbourhood arguments.

Minor precaution: with numparts = 1, Metis names the 1st partition 1 rather than 0, leading to confusion down the line. For simplicity, ignore in this case.
*/
void partitionMETIS(TetgenMesh *M, int nparts, int *edgecut_ptr, int **part_ptr) {

    int constraints = 1;
    *part_ptr =(int*)calloc(M->n * sizeof(int), sizeof(int));

    if (nparts > 1) METIS_PartGraphRecursive(&M->n, &constraints, M->xadj, M->adjncy, NULL, NULL, NULL, &nparts, NULL, NULL, NULL, edgecut_ptr, *part_ptr);
}

/*
This function does three things, which turn out to be somewhat related. In order of importance:
-generates a permutation array based on the partitioning, such that printing part[perm[1:n]] outputs 00...111...etc...(n-1). Separator nodes are placed higher.
-calculates the size of each parttition.
-identifies the separator size of each partition.
A sensible global relabeling of the node vertices is given by vertix #perm[i] --> #i; the local relabeling is then obvious.
Quicksort as-is may not be the ideal algorithm choice. In fact, given the data it may be an embarrasingly poor choice.
Output: perm, counters (size of each partition), separatorSizes (size of separator on each partition)

By separator vertices we mean vertices whose neighbours or neighbours of neighbours do not belong to the same partition.
*/

void reindexPrep(TetgenMesh *M, int numpart, int *part, int **perm_ptr, int **invperm_ptr, int **partsizes_ptr, int **separatorSizes_ptr) {
//void reindexPartitions(int n, int *xadj, int *adjncy, int numpart, int *part, int **perm_ptr, int **partsizes_ptr, int **separatorSizes_ptr) {
//void reindexPartitions(int n, int numpart, int **perm_ptr, int **partsizes_ptr, int **separatorSizes_ptr, int *part, int *xadj, int *adjncy) {

    int i, j, k;
    int *xadj, *adjncy, *perm, *invperm, *partsizes, *sortweight, *separatorSizes;

    xadj = M->xadj; adjncy = M->adjncy;

    *partsizes_ptr = partsizes = (int *)calloc(numpart * sizeof(int), sizeof(int));
    *separatorSizes_ptr = separatorSizes = (int *)calloc(numpart * sizeof(int), sizeof(int));
    *perm_ptr = perm = (int *)malloc(M->n * sizeof(int));
    *invperm_ptr = invperm = (int *)malloc(M ->n * sizeof(int));

    sortweight = (int *)malloc(M->n * sizeof(int));

    for (i = 0; i < M->n; i ++) {
        perm[i] = i;
        partsizes[part[i]] ++;
        sortweight[i] = (1 + part[i])*17;
    }

    // Identify separator vertices.
    #pragma omp parallell for
    for (i = 0; i < M->n; i ++) {
        for (j = xadj[i]; j < xadj[i + 1]; j ++) {

            if (part[i] != part[adjncy[j]]) sortweight[i] --;
            for (k = xadj[adjncy[j]]; k < xadj[adjncy[j] + 1]; k ++) {

                if (part[i] != part[adjncy[k]]) sortweight[i] --;
            }
        }
    }

    // The sortweight array is the best tool on hand to identify the separator sizes, so that's done here too...
    for (i = 0; i < M->n; i ++) {
        if (sortweight[i] % 17 != 0) separatorSizes[part[i]] ++;
    }

    quicksortPerm(sortweight, perm, 0, M->n);

    for (i = 0; i < M->n; i ++) invperm[perm[i]] = i;

    free(sortweight);

}

/*
Rearranges xadj, adjncy and the elements array based on a permutation array. This is not done in-place. The result is that the
xadj, adjncy and elements array are permuted to receive new global names.
*/
void rearrange(TetgenMesh *M, int *perm, int *invperm) {

    int i, j, size;
    int *xadj, *xadj_new, *adjncy, *adjncy_new, *elemts, *elemts_new;

    xadj = M->xadj;
    adjncy = M->adjncy;
    elemts = M->elements;
    xadj_new = malloc((M->n + 1) * sizeof(int));
    adjncy_new = malloc(xadj[M->n] * sizeof(int));

    xadj_new[0] = 0;
    for (i = 0; i < M->n; i ++) {
        size = xadj[perm[i] + 1] - xadj[perm[i]];
        xadj_new[i + 1] = xadj_new[i] + size;
    }

    for (i = 0; i < M->n; i ++) {
        size = xadj_new[i + 1] - xadj_new[i];
        for (j = 0; j < size; j ++) {
            adjncy_new[xadj_new[i] + j] = invperm[adjncy[xadj[perm[i]] + j]];
        }
    }

    M->xadj = xadj_new;
    M->adjncy = adjncy_new;

    free(xadj);
    free(adjncy);

    elemts_new = malloc(4 * M->n * sizeof(int));

    for (i = 0; i < M->n; i ++) {
        for (j = 0; j < 4; j ++)  {
            elemts_new[4 * i + j] = elemts[4 * perm[i] + j];
        }
    }

    M->elements = elemts_new;

    free(elemts);

}

/*
* Renames the nodes in the adjncy array to their appropriate local names.
*/
/*void renameNeighbours(TetgenMesh *M, int numpart, int *part, int *partsizes) {

    int i, j;

    for (i = 0; i < numpart; i ++) {

        for (j = 0; j < partsizes; j ++) {

        }
    }

}*/

/*
* Based on the sorted, globally namedTetgenMesh, this function names the neighbours of the separator vertices.

* The output is in ELL format, such that the 1st (16 * separatorsizes[0]) corresponds to the entire
* separator belonging to partition 0, etc, and each separator vertex has 16 neighbours^2--4 neighbours and 12 neighbours^2
* (because vertix i is a neighbour^2 of itself, and so isn't counted). All indices not used in this array are set to -1.

* Quadruple for loops and redundancies ahoy.

* At this point, consider sorting the part array, or just find appropriate sortedpart[j] via partsizes or something.

* Holy shit this is a monstrous abomination though.
*/
void identifyNeighbours(TetgenMesh *M, int **allseparators_ptr, int *separatorsizes, int *partsizes, int numpart, int *part, int *invperm) {

    int i, j, k, l, ll, sum;
    int currentidx_base, currentidx, sepidx, currentneigh, currentneighsq;
    int *allseparators;

    int *xadj = M->xadj;
    int *adjncy = M->adjncy;

    sum = 0;
    for (i = 0; i < numpart; i ++) {
        sum += separatorsizes[i];
    }
    *allseparators_ptr = allseparators = (int *)malloc(16 * sum * sizeof(int));
    for (i = 0; i < 16 * sum; i ++) allseparators[i] = -1;

    currentidx_base = currentidx =  sepidx = 0;
    for (i = 0; i < numpart; i ++) {

        for (j = 0; j < separatorsizes[i]; j ++) {

            for (k = 0; k < xadj[currentidx + 1] - xadj[currentidx]; k ++) { // <-- NEIGHBOURS loop

                currentneigh = adjncy[xadj[currentidx] + k];

                allseparators[sepidx + k] = part[invperm[currentneigh]];

                ll = 0;
                for (l = 0; l < xadj[currentneigh + 1] - xadj[currentneigh]; l ++) {// NEIGHBOURS OF NEIGHBOURS loop

                    currentneighsq = adjncy[xadj[currentneigh] + l];

                    if (currentneighsq != currentidx) allseparators[sepidx + 4 + 3*k + ll++] = part[invperm[currentneighsq]];
                }
            }

            currentidx ++;
            sepidx += 16;
        }
        currentidx_base += 16 * partsizes[i];
        currentidx = currentidx_base;
    }

}


/*
Using the above functions, obtains all the Tetgen data and sorts it for distribution.
*/
void readAndPrepareData(char *filename, int numpart, TetgenMesh *M, int **partsizes_ptr, int **separatorsizes_ptr) {

    int edgecut;
    int *part, *perm, *invperm, *allseparators;

    readmesh(M, filename);
    partitionMETIS(M, numpart, &edgecut, &part);
    reindexPrep(M, numpart, part, &perm, &invperm, partsizes_ptr, separatorsizes_ptr);
    rearrange(M, perm, invperm);
    //identifyNeighbours(M, &allseparators, *separatorsizes_ptr, *partsizes_ptr, numpart, part, invperm);

    free(part);
    free(perm);
    free(invperm);
}