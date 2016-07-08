/*
* The eternal question: how many l's in parallell?
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

void readfiles(char *filename, int myrank, int numprocs, TetgenMesh *M_local) {

    int i, j, n, numnodes, local_n, local_numnodes;
    int *xadj, *adjncy, *elements;
    double *nodes;
    int edgesLocal;
    char fname[200], suffix[20];
    FILE initdata;
    MPI_FILE *fpdata;

    int dump, read;
    int err, abort;

    int *scan;

    // Seems best to demand an extra parallell file containing node counts, element counts

    if (myrank == 0) {

        strcpy(fname, infile);
        strcpy(suffix, "_init.para");
        strcat(fname, suffix);

        initdata = fopen(fname, "r");

        if (fpdata == NULL) {
            printf("\nFailure to open input file: %s.\n", fname);
            abort = 1;
        }


        err = fscanf(fpdata, "%d", &n);

    }

    // MPI_Bcast(blah blah blah)


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

    M->n = n;

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