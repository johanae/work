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
        sortweight[i] = (1 + part[i]) * (numpart + 3);
    }

    limits[0] = 0;
    for (int i = 1; i <= numpart; i ++) limits[i] = limits[i - 1] + partsizes[i - 1];

    // Identify separator vertices.
    // #pragma omp parallell for
    int n, nn;
    int leading, weight, doublebreak;

    for (int i  = 0; i < M->n; i ++) {

        leading = -1; weight = 0; doublebreak = 0;

        for (int j = 0; j < 4; j ++) {
            n = M->ell[4 * i + j];
            if (n == -1) continue;
            if (part[i] != part[n]) {
                if (leading == -1) {    // 1st external node encountered
                    leading = part[n];
                    weight = numpart - part[n] + 2;
                } else if (leading != part[n]) {    // multiple neighbours. Put it at the bottom.
                        weight = 1;
                        break;
                }
            }

            for (int k = 0; k < 4; k ++) {

                nn = M->ell[4 * n + k];
                if (nn == -1) continue;
                if (part[i] != part[nn]) {
                    if (leading == -1) {    // 1st external node encountered
                        leading = part[nn];
                        weight = numpart - part[nn] + 2;
                    } else if (leading != part[nn]) {
                            weight = 1;
                            doublebreak = 1;
                            break;
                    }
                }

            }
            if (doublebreak) break;
        }

        sortweight[i] -= weight;
    }


    // The sortweight array is the best tool on hand to identify the separator sizes, so that's done here too.
    for (int i = 0; i < M->n; i ++) {
        if (sortweight[i] % (numpart + 3) != 0) separatorsizes[part[i]] ++;
    }