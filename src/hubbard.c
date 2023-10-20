#include "hubbard.h"
#include "readgraphmllib.h"

int get_matelem(size_t deti, size_t detj) {
  exc_number_t exij;
  determinant_t d1[1];
  determinant_t d2[1];
  d1[0] = deti;
  d2[0] = detj;
  exij = exc_degree(1, d1, d2);
  return exij;
}

// Function to calculate binomial coefficient using lgamma function
long long binomialCoeff(int n, int k) {
    return round(exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)));
}

void printBits(int num, int len) {
    for (int bit = len - 1; bit >= 0; --bit) {
        printf("%d", (num >> bit) & 1);
    }
    printf("\n");
}

void generateConfigurations(int norb, int nelec, int* configAll, int* size) {
    *size = 0;

    for (int i = 0; i < (1 << norb); ++i) {
        if (__builtin_popcount(i) == nelec) {
            configAll[(*size)++] = i;
        }
    }
}

// Function to compare two configurations for qsort and bsearch
int compare(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

// Function to find the positions of a list of configurations in a sorted list
void findPositions(int* configList, int sizeList, int* configs, int sizeConfigs, int* positions) {
    for (int i = 0; i < sizeConfigs; ++i) {
        int* item = (int*) bsearch(&configs[i], configList, sizeList, sizeof(int), compare);
        if (item != NULL) {
            positions[i] = item - configList;
        } else {
            positions[i] = -1;
        }
    }
}

void printPositions(int* positions, int size) {
    for (int i = 0; i < size; ++i) {
        printf("Position: %d\n", positions[i]);
    }
}

// Function to generate all possible alpha determinants
void generateAlphaDeterminants(int* configAlpha, int sizeAlpha, const igraph_t* graph, size_t alphaConfig, igraph_vector_t* alphaDeterminants) {
    // Get the number of orbitals
    int norb = igraph_vcount(graph);

    // Loop over each orbital
    for (int i = 0; i < norb; ++i) {
        // Check if the orbital is occupied
        if ((alphaConfig >> i) & 1) {
            // Get the connected vertices
            igraph_vector_t orbital_id_allowed;
            igraph_vector_init(&orbital_id_allowed, 0);
            getConnectedVertices(graph, i, &orbital_id_allowed);

            // Loop over each connected vertex
            for (int j = 0; j < igraph_vector_size(&orbital_id_allowed); ++j) {
                int orbital_id = VECTOR(orbital_id_allowed)[j];

                // Check if the connected vertex is unoccupied
                if (!((alphaConfig >> orbital_id) & 1)) {
                    // Create a new alpha determinant by moving the electron
                    size_t newAlphaConfig = alphaConfig ^ ((1 << i) | (1 << orbital_id));

                    // Find the position of the new alpha determinant in the list and add it to alphaDeterminants
                    int pos;
                    findPositions(configAlpha, sizeAlpha, &newAlphaConfig, 1, &pos);

                    // Add the position of the new alpha determinant to the list
                    igraph_vector_push_back(alphaDeterminants, pos);
                }
            }

            igraph_vector_destroy(&orbital_id_allowed);
        }
    }
}

int* igraphVectorToIntArray(const igraph_vector_t* igraph_vector) {
    int* int_array = (int*)malloc(igraph_vector_size(igraph_vector) * sizeof(int));

    if (int_array == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }

    for (igraph_integer_t i = 0; i < igraph_vector_size(igraph_vector); i++) {
        int_array[i] = (int)VECTOR(*igraph_vector)[i];
    }

    return int_array;
}

int getPhase(size_t alphaConfig, size_t newAlphaConfig, int h, int p) {

    // Phase
    unsigned int nperm;

    determinant_t d1[1];
    determinant_t d2[1];
    d1[0] = alphaConfig;
    d2[0] = newAlphaConfig;
    orbital_t h1[1];
    orbital_t p2[1];
    h1[0] = h;
    p2[0] = p;
    nperm = get_nperm_single((unsigned int) 1, d1, d2, h1, p2);
    int phase = ((unsigned int) 1) & nperm;
    //printf(" %llu %llu (%d, %d) nperm = %d phase=%d \n",d1[0], d2[0], i,orbital_id,nperm,phase);
    return phase;
}
