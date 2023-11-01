#include "hubbard.h"
#include "readgraphmllib.h"

size_t get_matelem(size_t deti, size_t detj) {
  exc_number_t exij;
  determinant_t d1[1];
  determinant_t d2[1];
  d1[0] = deti;
  d2[0] = detj;
  exij = exc_degree(1, d1, d2);
  return exij;
}

// Function to find a global ID in the hash table given an alpha and beta ID
size_t findGlobalID(size_t alphaID, size_t betaID, size_t nalpha) {
    size_t id = (alphaID-1)*nalpha + betaID;
    return id;
}

// Function to calculate binomial coefficient using lgamma function
long long binomialCoeff(size_t n, size_t k) {
    return round(exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)));
}

void printBits(size_t num, size_t len) {
    for (int bit = len - 1; bit >= 0; --bit) {
        printf("%ld", (num >> bit) & 1);
    }
    printf("\n");
}

void generateConfigurations(size_t norb, size_t nelec, size_t* configAll, size_t* size) {
    *size = 0;

    for (size_t i = 0; i < (1 << norb); ++i) {
        if (__builtin_popcount(i) == nelec) {
            configAll[(*size)++] = i;
        }
    }
}

// Function to compare two configurations for qsort and bsearch
int compare(const void* a, const void* b) {
    return (*(size_t*)a - *(size_t*)b);
}

// Function to find the positions of a list of configurations in a sorted list
void findPositions(size_t* configList, size_t sizeList, size_t* configs, size_t sizeConfigs, size_t* positions) {
    for (size_t i = 0; i < sizeConfigs; ++i) {
        size_t* item = (size_t*) bsearch(&configs[i], configList, sizeList, sizeof(size_t), compare);
        if (item != NULL) {
            positions[i] = item - configList;
        } else {
            positions[i] = -1;
        }
    }
}

void printPositions(size_t* positions, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        printf("Position: %ld\n", positions[i]);
    }
}

// Function to generate all possible alpha determinants
void generateDeterminants(size_t* configAlpha, size_t sizeAlpha, const igraph_t* graph, size_t alphaConfig, igraph_vector_t* alphaDeterminants) {
    // Get the number of orbitals
    size_t norb = igraph_vcount(graph);

    // Loop over each orbital
    for (size_t i = 0; i < norb; ++i) {
        // Check if the orbital is occupied
        if ((alphaConfig >> i) & 1) {
            // Get the connected vertices
            igraph_vector_t orbital_id_allowed;
            igraph_vector_init(&orbital_id_allowed, 0);
            getConnectedVertices(graph, i, &orbital_id_allowed);

            // Loop over each connected vertex
            for (size_t j = 0; j < igraph_vector_size(&orbital_id_allowed); ++j) {
                size_t orbital_id = VECTOR(orbital_id_allowed)[j];

                // Check if the connected vertex is unoccupied
                if (!((alphaConfig >> orbital_id) & 1)) {
                    // Create a new alpha determinant by moving the electron
                    size_t newAlphaConfig = alphaConfig ^ ((1 << i) | (1 << orbital_id));

                    // Find the position of the new alpha determinant in the list and add it to alphaDeterminants
                    size_t pos;
                    findPositions(configAlpha, sizeAlpha, &newAlphaConfig, 1, &pos);

                    // Add the position of the new alpha determinant to the list
                    igraph_vector_push_back(alphaDeterminants, pos);
                }
            }

            igraph_vector_destroy(&orbital_id_allowed);
        }
    }
}

size_t* igraphVectorToIntArray(const igraph_vector_t* igraph_vector) {
    size_t* int_array = (size_t*)malloc(igraph_vector_size(igraph_vector) * sizeof(size_t));

    if (int_array == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }

    for (igraph_integer_t i = 0; i < igraph_vector_size(igraph_vector); i++) {
        int_array[i] = (size_t)VECTOR(*igraph_vector)[i];
    }

    return int_array;
}

size_t getPhase(size_t alphaConfig, size_t newAlphaConfig, size_t h, size_t p) {

    // Phase
    size_t nperm;

    determinant_t d1[1];
    determinant_t d2[1];
    d1[0] = alphaConfig;
    d2[0] = newAlphaConfig;
    orbital_t h1[1];
    orbital_t p2[1];
    h1[0] = h;
    p2[0] = p;
    nperm = get_nperm_single((size_t) 1, d1, d2, h1, p2);
    size_t phase = ((size_t) 1) & nperm;
    //printf(" %llu %llu (%d, %d) nperm = %d phase=%d \n",d1[0], d2[0], i,orbital_id,nperm,phase);
    return phase;
}

// Function to generate all possible alpha determinants for a list of given alpha configurations
void generateAllDeterminants(size_t *configAlpha, size_t sizeAlpha, const igraph_t* graph, size_t* alphaConfigs, size_t numConfigs, igraph_vector_t* allAlphaDeterminants) {
    for (size_t i = 0; i < numConfigs; ++i) {
        igraph_vector_t alphaDeterminants;
        igraph_vector_init(&alphaDeterminants, 0);

        generateDeterminants(configAlpha, sizeAlpha, graph, alphaConfigs[i], &alphaDeterminants);

        for (size_t j = 0; j < igraph_vector_size(&alphaDeterminants); ++j) {
            igraph_vector_push_back(allAlphaDeterminants, VECTOR(alphaDeterminants)[j]);
        }

        igraph_vector_destroy(&alphaDeterminants);
    }
}
