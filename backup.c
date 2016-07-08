void storeData(char *fname_common, Nodemeshdata *d) {
/*void storeData(char *fname_common, int myrank, int numprocs, int *limits, int *separatorsizes, int local_n, int local_numnodes, int *local_ell, int *local_elmts, double *local_nodes,
               int separatorsize, int *separator, int *seppart) { */

    char fname[200];
    char suffix[20];
    char buffer[200];
    FILE *file;
    int err;

    if (myrank == 0) {

        strcpy(fname, fname_common);
        strcpy(suffix,"_init.txt");
        strcat(fname, suffix);
        file = fopen(fname, "w");    // possibly stick to nonbinary format for the init file

        fprintf(file, "%d\n\n", limits[numprocs]);
        //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);

        for (int i = 0; i <= numprocs; i ++) {

            fprintf(file, "%d %d\n", i, limits[i]);
            //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);
        }

        fprintf(file, "\n");
        //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);

        for (int i = 0; i < numprocs; i ++) {

            fprintf(file, "%d %d\n", i, separatorsizes[i]);
            //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);
        }

        fclose(file);
    }

    sprintf(suffix, "%04d.txt", myrank);
    strcpy(fname, fname_common);
    strcat(fname, suffix);

    file = fopen(fname, "w");

    fprintf(file, "\nNeighbours:\n");

    for (int i = 0; i < local_n; i ++) {

        if (i == separatorsize) fprintf(file, "----------------------------------------------------------------------------------------\n");
        fprintf(file, "%d %d %d %d %d\n", i, local_ell[4 * i], local_ell[4 * i + 1], local_ell[4 * i + 2], local_ell[4 * i + 3]);
        //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);
    }

    fprintf(file, "\nElements:\n");
    for (int i = 0; i < local_n; i ++) {

        fprintf(file, "%d %d %d %d %d\n", i, local_elmts[4 * i], local_elmts[4 * i + 1], local_elmts[4 * i + 2], local_elmts[4 * i + 3]);
        //err = fwrite(buffer, sizeof(char), sizeof(buffer), file);
    }

    fprintf(file, "\nNodes:\n");
    for (int i = 0; i < local_numnodes; i ++) {

        fprintf(file, "%d %lf %lf %lf\n", i, local_nodes[3 *i ], local_nodes[3 * i + 1], local_nodes[3 * i + 2]);
    }

    fprintf(file, "\nSeparator:\n");
    for (int i = 0; i < separatorsize; i ++) {
        for (int j = 0; j < 16; j ++) {
            fprintf(file, "%d ", separator[16 * i + j]);
        }
        fprintf(file, "   ");
        for (int j = 0; j < 16; j ++) {
            fprintf(file, "%d ", seppart[16 * i + j]);
        }
        fprintf(file, "\n");
    }


    fclose(file);

}

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

#include <limits.h>
#include "readutility.h"
#include "common.h"
#include "kaHIP_interface.h"
#include "metis.h"


/* Reads a .neigh input file to a CSR type format. This is a temporary function, subject to expansion and renaming. */
/* Arguments: (n, xadj, adjncy), filename. */
/* Input: .neigh file, obtained via Tetgen. */
/* Output: integer n, xadj, adjncy arrays.*/
void readmesh(TetgenMesh *M, char *infile) {
/*void readmesh(int *n, int **xadj_ptr, int **adjncy_ptr, int *numnodes, int **nodes_ptr, double **nodeVals_ptr, char *infile) {*/

    int i, j, n, numnodes;
    int *xadj, *adjncy, *ell, *elements;
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

    M->n = n;
    ell = (int *)calloc(4 * n * sizeof(int), sizeof(int));

    for(int j = 0; j < 1; ++j) {
        err = fscanf(fpdata, "%d", &dump);
    }

    for (int i = 0; i < n; i++) {

        err = fscanf(fpdata, "%d", &dump);

        for (int j = 0; j < 4; j ++) {
            fscanf(fpdata, "%d", &ell[4*i + j]);
        }
    }

    fclose(fpdata);

    ell2csr(n, ell, &xadj, &adjncy);
    M->ell = ell;
    M->xadj = xadj;
    M->adjncy = adjncy;

    printf("ELL data read!\n");

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
void partitionMETIS(TetgenMesh *M, int *edgecut_ptr, int **part_ptr) {

    int constraints = 1;
    *part_ptr =(int*)calloc(M->n * sizeof(int), sizeof(int));

    if (M->numpart > 1) METIS_PartGraphRecursive(&M->n, &constraints, M->xadj, M->adjncy, NULL, NULL, NULL, &M->numpart, NULL, NULL, NULL, edgecut_ptr, *part_ptr);
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

void reindexPrep(TetgenMesh *M, int numpart, int *part, int **perm_ptr, int **invperm_ptr) {
/*void reindexPrep(TetgenMesh *M, int numpart, int *part, int **perm_ptr, int **invperm_ptr, int **partsizes_ptr, int **limits_ptr,
                 int **separatorsizes_ptr, int **separraysizes_ptr, int **seplimits_ptr) {*/

    int *xadj, *adjncy, *perm, *invperm, *partsizes, *limits, *sortweight, *separatorsizes, *separraysizes, *seplimits;

    xadj = M->xadj; adjncy = M->adjncy;

    limits = (int *)malloc((numpart + 1) * sizeof(int));
    partsizes = (int *)calloc(numpart * sizeof(int), sizeof(int));
    separatorsizes = (int *)calloc(numpart * sizeof(int), sizeof(int));
    separraysizes= (int *)malloc(numpart * sizeof(int));
    seplimits = (int *)malloc((numpart + 1) * sizeof(int));

    sortweight = (int *)malloc(M->n * sizeof(int));
    *perm_ptr = perm = (int *)malloc(M->n * sizeof(int));
    *invperm_ptr = invperm = (int *)malloc(M ->n * sizeof(int));


    for (int i = 0; i < M->n; i ++) {
        perm[i] = i;
        partsizes[part[i]] ++;
        sortweight[i] = (1 + part[i])*17;
    }

    limits[0] = 0;
    for (int i = 1; i <= numpart; i ++) limits[i] = limits[i - 1] + partsizes[i - 1];

    // Identify separator vertices.
    // #pragma omp parallell for
    int n, nn;
    for (int i  = 0; i < M->n; i ++) {

        for (int j = 0; j < 4; j ++) {

            n = M->ell[4 * i + j];
            if (n == -1) continue;
            if (part[i] != part[n]) sortweight[i] --;

            for (int k = 0; k < 4; k ++) {

                nn = M->ell[4 * n + k];
                if (nn == -1) continue;
                if (part[i] != part[nn]) sortweight[i] --;

            }
        }
    }


    // The sortweight array is the best tool on hand to identify the separator sizes, so that's done here too.
    for (int i = 0; i < M->n; i ++) {
        if (sortweight[i] % 17 != 0) separatorsizes[part[i]] ++;
    }

    seplimits[0] = 0;
    for (int i = 1; i <= numpart; i ++) {
        separraysizes[i - 1] = 16 * separatorsizes[i - 1];
        seplimits[i] = seplimits[i - 1] + separraysizes[i - 1];
    }

    quicksortPerm(sortweight, perm, 0, M->n);

    for (int i = 0; i < numpart; i ++) {

    }

    for (int i = 0; i < M->n; i ++) invperm[perm[i]] = i;

    M->limits = limits;
    M->partsizes = partsizes;
    M->separatorsizes = separatorsizes;
    M->separraysizes = separraysizes;
    M->seplimits = seplimits;

    free(sortweight);

}

/*
Rearranges xadj, adjncy and the elements array based on a permutation array.
*/
void rearrange(TetgenMesh *M, int *perm, int *invperm) {

    int *ell, *ell_new, *elemts, *elemts_new;

    ell = M->ell;
    elemts = M->elements;
    ell_new = malloc(M->n * 4 * sizeof(int));

    for (int i = 0; i < M->n; i ++) {

        for (int j = 0; j < 4; j ++) {

            ell_new[4 * i + j] = ell[4 * perm[i] + j];
        }

        //memcpy(&ell_new[4 * perm[i]], &ell[4 * i], 4 * sizeof(int));
    }

    M->ell = ell_new;

    free(ell);

    elemts_new = malloc(4 * M->n * sizeof(int));

    for (int i = 0; i < M->n; i ++) {
        for (int j = 0; j < 4; j ++)  {
            elemts_new[4 * i + j] = elemts[4 * perm[i] + j];    // switch to memcpy
        }
    }

    M->elements = elemts_new;

    free(elemts);

}

/*
* Renames the nodes in the ell array to their appropriate local names, and sets up the separator.

* This is done in the serial part of the program because, and this is the clever part, I mean just really brilliant, can I get back to you on this?
*/
void renameNeighbours(TetgenMesh *M, int numpart, int *part, int *invperm) {
//void renameNeighbours(TetgenMesh *M, int numpart, int *part, int *partsizes, int *limits, int **separators_ptr, int **sepparts_ptr, int *separatorsizes, int *seplimits, int *invperm) {

    int *ell = M->ell;
    int oldid, newid, oldid2, newid2, currentpart;
    int *separators, *sepparts;

    separators = (int *)malloc(16 * M->seplimits[M->numpart] * sizeof(int));
    sepparts = (int *)malloc(16 * M->seplimits[M->numpart] * sizeof(int));

    for (int i = 0; i < 16 * M->seplimits[M->numpart]; i ++) {

        separators[i] = -1;
        sepparts[i] = -1;
    }

    for (int i = 0; i < numpart; i ++) {

        for (int j = M->limits[i]; j < M->limits[i + 1]; j ++) {

            int l = j - M->limits[i];

            for (int k = 0; k < 4; k ++) {

                oldid = ell[4 * j + k];
                newid = invperm[oldid];

                if (oldid == -1) continue;

                if (l < M->separatorsizes[i]) {    // we're inside the separator belonging to partition i. This is as far down the rabbit hole as we go, for now.

                    newid -= M->limits[part[oldid]];
                    separators[M->seplimits[i] + 16 * l + 4 * k] = newid;
                    sepparts[M->seplimits[i] +16 * l + 4 * k] = part[oldid];

                } else {
                    newid -= M->limits[i];
                }

                ell[4 * j + k] = newid;
            }
        }
    }

    M->separators = separators;
    M->sepparts = sepparts;

}

/*
Based on sorted and renamed data, expand the separators to include neighbours of neighbours. The original separators are traversed, and the neighbours are
queried for information.
*/
void fillSeparator(TetgenMesh *M) {
//void fillSeparator(int *separators, int *sepparts, int *seplimits, int *separatorsizes, int *partsizes, int numpart) {

int current, currentpart, nid, npart, nnid, nnpart;

current = -1; currentpart = 0;
for (int i = 0; i < M->seplimits[M->numpart]; i += 16) {

    current ++;
    if (current >= M->separatorsizes[currentpart]) {
        current = 0;
        currentpart ++;
    }

    for (int j = 0; j < 16; j += 4) {

        nid = M->separators[i + j];  // 16 * i + 4 * j or w/e
        npart = M->sepparts[i + j];

        if (nid == -1) continue;

        if (npart == currentpart && nid > M->separatorsizes[currentpart]) {

            for (int k = 1; k < 4; k ++) {

                M->separators[i + j + k] = -2;
                M->sepparts[i + j + k] = currentpart;
            }
        } else {
            int l = 1;
            for (int k = 0; k < 4; k ++) {

                nnid = M->separators[M->seplimits[npart] + 16 * nid + 4 * k]; // I think this should be correct.
                nnpart = M->sepparts[M->seplimits[npart] + 16 * nid + 4 * k];

                if (nnid == -1 || (nnid == current && nnpart == currentpart)) continue;

                M->separators[i + j + l] = nnid;
                M->sepparts[i + j + l++] = nnpart;

                if (l == 5) printf("Disaster has struck. Current: %d of %d, n: %d of %d. nn: %d of house %d \n", current, currentpart, nid, npart, nnid, nnpart);
            }
        }

    }

}

}


/*
Using the above functions, obtains all the Tetgen data and sorts it for distribution.
*/
void readAndPrepareData(char *filename, int numpart, TetgenMesh *M) {
//void readAndPrepareData(char *filename, int numpart, TetgenMesh *M, int **partsizes_ptr, int **limits_ptr,
                        //int **separatorsizes_ptr, int **separraysizes_ptr, int **seplimits_ptr, int **separators_ptr, int **sepparts_ptr) {

    M->numpart = numpart;

    int edgecut;
    int *part, *perm, *invperm;

    readmesh(M, filename);
    partitionMETIS(M, &edgecut, &part);
    reindexPrep(M, numpart, part, &perm, &invperm);
    rearrange(M, perm, invperm);
    //void renameNeighbours(TetgenMesh *M, int numpart, int *part, int *invperm)
    renameNeighbours(M, numpart, part, invperm);
    //renameNeighbours(M, numpart, part, *partsizes_ptr, *limits_ptr, separators_ptr, sepparts_ptr, *separatorsizes_ptr, *seplimits_ptr, invperm);
    fillSeparator(M);


    free(part);
    free(perm);
    free(invperm);
}

void fillSeparator(TetgenMesh *M) {
//void fillSeparator(int *separators, int *sepparts, int *seplimits, int *separatorsizes, int *partsizes, int numpart) {

    int current, currentpart, nid, npart, nnid, nnpart;

    current = -1; currentpart = 0;
    for (int i = 0; i < M->seplimits[M->numpart]; i += 16) {

        current ++;
        if (current >= M->separatorsizes[currentpart]) {
            current = 0;
            currentpart ++;
        }

        for (int j = 0; j < 16; j += 4) {

            nid = M->separators[i + j];  // 16 * i + 4 * j or w/e
            npart = M->sepparts[i + j];

            if (nid == -1) continue;

            if (npart == currentpart && nid > M->separatorsizes[currentpart]) {

                for (int k = 1; k < 4; k ++) {

                    M->separators[i + j + k] = -2;
                    M->sepparts[i + j + k] = currentpart;
                }
            } else {
                int l = 1;
                for (int k = 0; k < 4; k ++) {

                    nnid = M->separators[M->seplimits[npart] + 16 * nid + 4 * k]; // I think this should be correct.
                    nnpart = M->sepparts[M->seplimits[npart] + 16 * nid + 4 * k];

                    if (nnid == -1 || (nnid == current && nnpart == currentpart)) continue;

                    M->separators[i + j + l] = nnid;
                    M->sepparts[i + j + l++] = nnpart;

                    if (l == 5) printf("Something's wrong with the ordering. Current: %d of %d, n: %d of %d. nn: %d of house %d \n", current, currentpart, nid, npart, nnid, nnpart);
                }
            }

        }

    }

}