#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <limits.h>

#include "common.h"
#include "readutility.h"

// Pending cleanup.
/*
Compute initial data for creating the MPI graph topology communicators.
Currently, don't do any renaming. Hopefully MPI_Alltoall_Neighbour respects the original edges array.

General idea: Permsort seppart array, maybe trim neighcount. Make a sort of CSR matrix out of the data?

Names:
sep, seppart ~ separator, seppart arrays
d->sepsize ~ length of separator (not of the array)
edges ~ list of neighbours in increasing order, used for generating topology

All-to-all communication variables:
*sendsizes
*senddisps
*recvsizes
*recvdisps
*/
void initialiseTopologies(Nodemeshdata *d, int *trimmedcount_ptr, int **sep_trimmed_ptr, int **seppart_trimmed_ptr, int **orig_trimmed_ptr) {

    int trimmedcount, counter, nn, nneighcount, *nneighs;
    int *sep_trimmed, *seppart_trimmed, *ind_trimmed, *owner;

    nneighs = (int *)calloc(d->numprocs * sizeof(int), sizeof(int));

    // Traverse to find out How many neighbours a node has, and that other thing.
    nneighcount = trimmedcount = 0;
    for (int i = 0; i < 20 * d->sepsize; i ++) {

        nn = d->seppart[i];
        if (nn == -1 || nn == d->myrank) continue;

        trimmedcount ++;
        if (nneighs[nn]++ < 1) nneighcount ++;

    }

    d->neighcount = nneighcount;
    d->senddisps = (int *)malloc((nneighcount + 1) * sizeof(int));
    d->edges = (int *)malloc(nneighcount * sizeof(int));
    counter = d->senddisps[0] = 0;
    for (int i = 0; i < d->numprocs; i ++) {
        if(i != d->myrank && nneighs[i] > 0) {
            d->edges[counter++] = i;
            d->senddisps[counter] = d->senddisps[counter - 1] + nneighs[i];
        }
    }

    ind_trimmed = (int *)malloc(trimmedcount * sizeof(int));
    *orig_trimmed_ptr = owner = (int *)malloc(trimmedcount * sizeof(int));
    // Traverse slightly more carefully, make lists.
    counter = 0;
    for (int i = 0; i < 20 * d->sepsize; i ++) {

        nn = d->seppart[i];
        if (nn == -1 || nn == d->myrank) continue;
        owner[counter] = i / 20;
        ind_trimmed[counter++] = i;
    }

    seppart_trimmed = (int *)malloc(trimmedcount * sizeof(int));
    int *perm = (int *)malloc(trimmedcount * sizeof(int));
    for (int i = 0; i < trimmedcount; i ++) {
        seppart_trimmed[i] = d->seppart[ind_trimmed[i]];
        perm[i] = i;
    }

    quicksortPerm(seppart_trimmed, perm, 0, trimmedcount);

    int *perm2 = (int *)malloc(trimmedcount * sizeof(int));
    memcpy(perm2, perm, trimmedcount * sizeof(int));

    sep_trimmed = (int *)malloc(trimmedcount * sizeof(int));
    for (int i = 0; i < trimmedcount; i ++) {

        sep_trimmed[i] = d->sep[ind_trimmed[i]];
    }

    for (int i = 0; i < nneighcount; i ++) {
        quicksortPerm(sep_trimmed, perm, d->senddisps[i], d->senddisps[i + 1]);
        quicksortPerm(owner, perm2, d->senddisps[i], d->senddisps[i + 1]);
    }

    int *trimmed_sorted = (int *)malloc(trimmedcount * sizeof(int));
    int *part_ts = (int *)malloc(trimmedcount * sizeof(int));
    int *owner_sorted = (int *)malloc(trimmedcount * sizeof(int));
    for (int i = 0; i < trimmedcount; i ++) {
        trimmed_sorted[i] = sep_trimmed[perm[i]];
        part_ts[i] = seppart_trimmed[perm[i]];  // ditch, I think
        owner_sorted[i] = owner[perm2[i]];
    }

    //if (d->myrank == 0) for (int i = 0; i < trimmedcount; i ++) printf("Receive: %d from %d\n", trimmed_sorted[i], seppart_trimmed[perm[i]]);
    //if (d->myrank == 0) for (int i = 0; i < trimmedcount; i ++) printf("Send: %d to %d\n", owner_sorted[i], seppart_trimmed[perm[i]]);    // Holy shit, I can't believe this actually worked.

    d->sendsizes = (int *)malloc(d->neighcount * sizeof(int));
    d->senddisps = (int *)malloc(d->neighcount * sizeof(int));
    // Trim further to obtain recvsizes & disps! Those are pretty fucking important.
    int trimmedcount2 = trimmedcount;
    for (int i = 1; i < trimmedcount; i ++) {

        if (trimmed_sorted[i] == trimmed_sorted[i - 1]) trimmedcount2 --;
    }

    int *sep_short = (int *)malloc(trimmedcount2 * sizeof(int));
    int *seppart_short = (int *)malloc(trimmedcount2 * sizeof(int));    // probably drop this. it's implied by sorting + recvsizes/disps.
    int *trimmed2counts = (int *)calloc(trimmedcount2 * sizeof(int), sizeof(int));  // probably drop this

    counter = -1;
    for (int i = 0; i < trimmedcount; i ++) {

        if (i == 0 || trimmed_sorted[i] != trimmed_sorted[i - 1]) {
            sep_short[++counter] = trimmed_sorted[i];
            seppart_short[counter] = seppart_trimmed[perm[i]];
        }
        trimmed2counts[counter] ++;

    }

    d->recvsizes = (int *)calloc(d->neighcount * sizeof(int), sizeof(int));
    counter = 0;
    for (int i = 0; i < trimmedcount2; i ++) {

        if (i != 0 && seppart_short[i] != seppart_short[i - 1]) counter ++;
        d->recvsizes[counter]++;
    }

    d->recvdisps = (int *)malloc((d->neighcount + 1)* sizeof(int));
    d->recvdisps[0] = 0;
    for (int i = 1; i <= d->neighcount; i ++) d->recvdisps[i] = d->recvdisps[i - 1] + d->recvsizes[i - 1];

    if (d->myrank == 0) for (int i = 0; i <= d->neighcount; i ++) printf("%d vs %d\n", d->recvdisps[i], trimmedcount2   );

    // Rename the appropriate indices in the separator.
    int part, shift;
    for (int i = 0; i < trimmedcount; i++) {

        part = binarySearch(seppart_trimmed[i], 0, d->neighcount, d->neighcount>>1, d->edges);
        shift = binarySearch(sep_trimmed[i], d->recvdisps[part], d->recvdisps[part + 1], (d->recvdisps[part] + d->recvdisps[part + 1])/2, sep_short);
        d->sep[ind_trimmed[i]] = d->n + shift;
    }

    free(perm); free(perm2);

}