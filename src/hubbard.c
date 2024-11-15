#include "hubbard.h"
#include "readgraphmllib.h"
#include "get_s2.h"

double solveQuad(double a, double b, double c) {
  double root1 = (-b + sqrt(b*b-4.*a*c) ) / (2.*a);
  return(root1);
}

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
    size_t id = alphaID*nalpha + betaID;
    return id;
}

// Functions to find alpha and beta IDs from globalID
size_t findAlphaID(size_t globalID, size_t nalpha, size_t nbeta) {
    size_t alphaid = (globalID / nalpha);
    return alphaid;
}

size_t findBetaID(size_t globalID, size_t nalpha, size_t nbeta) {
    size_t alphaid = (globalID / nalpha);
    size_t betaid = globalID - alphaid*nalpha;
    return betaid;
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
        if (popcnt (i) == nelec) {
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

int getPhase(size_t alphaConfig, size_t newAlphaConfig, size_t h, size_t p) {

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
    // size_t phase = ((size_t) 1) & nperm;
    //printf(" %llu %llu (%d, %d) nperm = %d phase=%d \n",d1[0], d2[0], i,orbital_id,nperm,phase);
    return nperm;
}

int getExecDegree(size_t detI, size_t detJ) {
    // Result
    size_t result;
    determinant_t d1[1];
    determinant_t d2[1];
    d1[0] = detI;
    d2[0] = detJ;
    result = exc_degree((size_t) 1, d1, d2);
    return result;
}

int getHoles_1ex(size_t detI, size_t detJ, size_t *holesOut) {
    // Result
    size_t nholes;

    determinant_t d1[1];
    determinant_t d2[1];
    d1[0] = detI;
    d2[0] = detJ;
    orbital_t holes[2];
    nholes = get_holes((size_t) 1, d1, d2, holes);
    holesOut[0] = holes[0];
    return nholes;
}

int getPart_1ex(size_t detI, size_t detJ, size_t *particlesOut) {
    // Result
    size_t nparticles;

    determinant_t d1[1];
    determinant_t d2[1];
    d1[0] = detI;
    d2[0] = detJ;
    orbital_t particles[2];
    nparticles = get_particles((size_t) 1, d1, d2, particles);
    particlesOut[0] = particles[0];
    return nparticles;
}

void getElecList(size_t detI, size_t *holesOut, size_t nelec) {
    
    // Initialize determinant
    determinant_t d1[1];
    d1[0] = detI;
    orbital_t particles[nelec];
    orbital_t res;
    // Get the orbital (or holes) list
    res = to_orbital_list((size_t) 1, d1, particles);
    for( size_t i=0;i<nelec; ++i ) {
      holesOut[i] = particles[i];
    }
}

// Function to generate all possible alpha determinants
void generateDeterminants(size_t* configAlpha, size_t sizeAlpha, const igraph_t* graph, size_t alphaConfig, size_t betaConfig, igraph_vector_t* alphaDeterminants, igraph_vector_t* alphaMEs, int alphaBeta, double** wmat) {
    // Get the number of orbitals
    size_t norb = igraph_vcount(graph);
    int phase = 1;
    double MEt = 0.0;

    // Loop over each orbital
    for (size_t i = 0; i < norb; ++i) {
        // Check if the orbital is occupied
        if ((alphaConfig >> i) & 1) {
            // Get the connected vertices
            igraph_vector_int_t orbital_id_allowed;
            igraph_vector_int_init(&orbital_id_allowed, 0);
            getConnectedVertices(graph, (igraph_integer_t)i, &orbital_id_allowed);

            // Loop over each connected vertex
            for (size_t j = 0; j < igraph_vector_int_size(&orbital_id_allowed); ++j) {
                size_t orbital_id = VECTOR(orbital_id_allowed)[j];
                //printf("weight = %f \n",wmat[i][orbital_id]);
                MEt = wmat[i][orbital_id];

                // Check if the connected vertex is unoccupied
                if (!((alphaConfig >> orbital_id) & 1)) {
                    // Create a new alpha determinant by moving the electron
                    size_t newAlphaConfig = alphaConfig ^ ((1 << i) | (1 << orbital_id));

                    // Find the phase
                    size_t alphabetadet = alphaConfig ^ betaConfig;
                    alphabetadet = alphabetadet & ~((size_t)1 << (orbital_id));
                    size_t newalphabetadet = newAlphaConfig ^ betaConfig;
                    newalphabetadet = newalphabetadet & ~((size_t)1 << (i));
                    phase = getPhase(alphabetadet, newalphabetadet, i+1, orbital_id+1);
                    phase = phase & 1 == 1 ? -1 : 1;
                    if(alphaBeta == 0) {
                        if ( orbital_id > i ) {
                            if( ((betaConfig >> i) & 1) == 1) phase *= -1;
                        }
                        else {
                            if( ((betaConfig >> orbital_id) & 1) == 1) phase *= -1;
                        }
                    }
                    else {
                        if ( orbital_id > i ) {
                            if( ((betaConfig >> orbital_id) & 1) == 1) phase *= -1;
                        }
                        else {
                            if( ((betaConfig >> i) & 1) == 1) phase *= -1;
                        }
                    }

                    // Find the position of the new alpha determinant in the list and add it to alphaDeterminants
                    size_t pos;
                    findPositions(configAlpha, sizeAlpha, &newAlphaConfig, 1, &pos);

                    // Add the position of the new alpha determinant to the list
                    igraph_vector_push_back(alphaDeterminants, pos);

                    // Add the position of the new alpha determinant to the list
                    igraph_vector_push_back(alphaMEs, phase*MEt);
                    //printf(" MEt = %10.5f\n",MEt);
                }
            }

            igraph_vector_int_destroy(&orbital_id_allowed);
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

// Function to generate all possible alpha determinants for a list of given alpha configurations
void generateAllDeterminants(size_t *configAlpha, size_t sizeAlpha, const igraph_t* graph, size_t* alphaConfigs, size_t numConfigs, igraph_vector_t* allAlphaDeterminants, double** wmat) {
    for (size_t i = 0; i < numConfigs; ++i) {
        igraph_vector_t alphaDeterminants;
        igraph_vector_init(&alphaDeterminants, 0);
        igraph_vector_t alphaMEs;
        igraph_vector_init(&alphaMEs, 0);

        generateDeterminants(configAlpha, sizeAlpha, graph, alphaConfigs[i], 0, &alphaDeterminants, &alphaMEs, 1, wmat);

        for (size_t j = 0; j < igraph_vector_size(&alphaDeterminants); ++j) {
            igraph_vector_push_back(allAlphaDeterminants, VECTOR(alphaDeterminants)[j]);
        }

        igraph_vector_destroy(&alphaDeterminants);
    }
}

// Main function that calculates MEs
void getAllHubbardMEs(size_t Idet, igraph_vector_t* MElist, igraph_vector_t* Jdetlist, size_t *configAlpha, size_t sizeAlpha, size_t *configBeta, size_t sizeBeta, const igraph_t* graph, double** wmat) {
    double phaseAlpha;
    double phaseBeta;
    //Find alpha and beta ids
    size_t alphaID = findAlphaID(Idet, sizeAlpha, sizeBeta);
    size_t betaID  = findBetaID(Idet, sizeAlpha, sizeBeta);

    // Find allowed excitations
    igraph_vector_t alphaDeterminants;
    igraph_vector_init(&alphaDeterminants, 0);
    igraph_vector_t alphaMEs;
    igraph_vector_init(&alphaMEs, 0);
    generateDeterminants(configAlpha, sizeAlpha, graph, configAlpha[alphaID], configBeta[betaID], &alphaDeterminants, &alphaMEs, 1, wmat);
    igraph_vector_t betaDeterminants;
    igraph_vector_init(&betaDeterminants, 0);
    igraph_vector_t betaMEs;
    igraph_vector_init(&betaMEs, 0);
    generateDeterminants(configBeta, sizeBeta, graph, configBeta[betaID], configAlpha[alphaID], &betaDeterminants, &betaMEs, 0, wmat);

    for (size_t j = 0; j < igraph_vector_size(&alphaDeterminants); ++j) {
        size_t alphaJ = VECTOR(alphaDeterminants)[j];
        phaseAlpha = VECTOR(alphaMEs)[j];

        size_t foundGlobalID = findGlobalID(alphaJ, betaID, sizeAlpha);

        igraph_vector_push_back(Jdetlist, foundGlobalID);
        igraph_vector_push_back(MElist, phaseAlpha);
    }
    for (size_t k = 0; k < igraph_vector_size(&betaDeterminants); ++k) {

        size_t betaK = VECTOR(betaDeterminants)[k];
        phaseBeta = VECTOR(betaMEs)[k];

        size_t foundGlobalID = findGlobalID(alphaID, betaK, sizeAlpha);

        igraph_vector_push_back(Jdetlist, foundGlobalID);
        igraph_vector_push_back(MElist, phaseBeta);
    }
    igraph_vector_destroy(&alphaDeterminants);
    igraph_vector_destroy(&alphaMEs);
    igraph_vector_destroy(&betaDeterminants);
    igraph_vector_destroy(&betaMEs);
}

// Main function that generates S2 operator
void getS2Operator(size_t Idet, igraph_vector_t* MElist, igraph_vector_t* Jdetlist, size_t *configAlpha, size_t sizeAlpha, size_t *configBeta, size_t sizeBeta, const igraph_t* graph, int natom, int natomax) {
    int phaseAlpha;
    int phaseBeta;
    int ideter[natomax];
    int ideter2[natomax];
    int kko, kok, kkio;
    int iiii, iii, ii;
    size_t iaa2;
    double xmat = 0.0;
    //Find alpha and beta ids
    size_t alphaID = findAlphaID(Idet, sizeAlpha, sizeBeta);
    size_t betaID  = findBetaID(Idet, sizeAlpha, sizeBeta);
    size_t alphadet = configAlpha[alphaID];
    size_t betadet = configBeta[betaID];

    getdet (Idet, ideter, configAlpha, sizeAlpha, configBeta,  sizeBeta, natom);
    //printf(" %d : ",Idet);
    //for(kko=0;kko<natom;++kko) printf(" %d ",ideter[kko]);
    //printf("\n");
    for( kko = 0; kko < natom; ++kko ) {
        if(ideter[kko] != 4 && ideter[kko] != 3) xmat=xmat+(3.0/4.0);
        for( kok = kko + 1; kok < natom; ++kok ) {
            if(ideter[kko] == 1 && ideter[kok] == 1){
                xmat=xmat+(1.0/2.0);
            }
            if(ideter[kko] == 2 && ideter[kok] == 2){
                xmat=xmat+(1.0/2.0);
            }
            if(ideter[kko] == 1 && ideter[kok] == 2){
                xmat=xmat-(1.0/2.0);
                for(kkio=0;kkio<=natom-1;kkio++){
                    ideter2[kkio]=ideter[kkio];
                }
                ideter2[kko]=2;
                ideter2[kok]=1;
                adr (ideter2, &iaa2,
                     configAlpha, sizeAlpha,
                     configBeta,  sizeBeta, natom);
                igraph_vector_push_back(Jdetlist, iaa2);
                igraph_vector_push_back(MElist, 1.0);
            }
            if(ideter[kko] == 2 && ideter[kok] == 1){
                xmat=xmat-(1.0/2.0);
                for(kkio=0;kkio<=natom-1;kkio++){
                    ideter2[kkio]=ideter[kkio];
                }
                ideter2[kko]=1;
                ideter2[kok]=2;
                adr (ideter2, &iaa2,
                     configAlpha, sizeAlpha,
                     configBeta,  sizeBeta, natom);
                igraph_vector_push_back(Jdetlist, iaa2);
                igraph_vector_push_back(MElist, 1.0);
            }
        }
    }
    igraph_vector_push_back(Jdetlist, Idet);
    igraph_vector_push_back(MElist, xmat);
}

// Get the diagonal part of the hubbard Hamiltonian
int getHubbardDiag(size_t Idet, size_t *configAlpha, size_t sizeAlpha, size_t *configBeta, size_t sizeBeta) {
    //Find alpha and beta ids
    size_t alphaID = findAlphaID(Idet, sizeAlpha, sizeBeta);
    size_t betaID  = findBetaID(Idet, sizeAlpha, sizeBeta);

    size_t alphadet = configAlpha[alphaID];
    size_t betadet  = configBeta[betaID];
    size_t doublyOcc = popcnt(alphadet & betadet);
    return doublyOcc;
}

// Get the diagonal part of the hubbard Hamiltonian
void getHubbardDiagVij(size_t Idet, size_t *configAlpha, size_t sizeAlpha, size_t *configBeta, size_t sizeBeta, const igraph_t* graph, double* alphaMEs, double** vmat) {
    // Get the number of orbitals
    size_t norb = igraph_vcount(graph);
    //Find alpha and beta ids
    size_t alphaID = findAlphaID(Idet, sizeAlpha, sizeBeta);
    size_t betaID  = findBetaID(Idet, sizeAlpha, sizeBeta);

    size_t alphaConfig = configAlpha[alphaID];
    size_t betaConfig  = configBeta[betaID];
    size_t alphabeta   = alphaConfig | betaConfig;

    double MEv = 0.0;
    alphaMEs[0] = 0.0;

    // Loop over each orbital
    for (size_t i = 0; i < norb; ++i) {
        // Check if the orbital is occupied
        //if ((alphabeta >> i) & 1) {
            // Get the connected vertices
            igraph_vector_int_t orbital_id_allowed;
            igraph_vector_int_init(&orbital_id_allowed, 0);
            getConnectedVertices(graph, (igraph_integer_t)i, &orbital_id_allowed);

            // Loop over each connected vertex
            for (size_t j = 0; j < igraph_vector_int_size(&orbital_id_allowed); ++j) {
                size_t orbital_id = VECTOR(orbital_id_allowed)[j];
                //printf("weight = %f \n",wmat[i][orbital_id]);
                MEv = vmat[i][orbital_id];

                // Check if the connected vertex is unoccupied
                if (i < orbital_id) {
                    // Create a new alpha determinant by moving the electron
                    size_t ni = (( alphaConfig >> i) & 1) + ((betaConfig >> i) & 1);
                    size_t nj = (( alphaConfig >> orbital_id) & 1) + ((betaConfig >> orbital_id) & 1);

                    // Find the position of the new alpha determinant in the list and add it to alphaDeterminants
                    size_t pos = Idet;

                    // Add the position of the new alpha determinant to the list
                    //igraph_vector_push_back(alphaDeterminants, pos);

                    // Add the position of the new alpha determinant to the list
                    alphaMEs[0] += MEv * (1.0-ni)*(1.0-nj);
                    //printf(" ME=%10.5f\n",alphaMEs[0]);
                }
            }

            igraph_vector_int_destroy(&orbital_id_allowed);
        //}
    }
}

// Get maximum neighbors
int getMaxNeighbors(const igraph_t* graph, size_t nsites) {
    int max_nbrs = 0;
    for(int i=0; i<nsites; ++i) {
        int nbrs = getNumberOfConnectedVertices(graph, (igraph_integer_t)i);
        if(nbrs > max_nbrs) max_nbrs = nbrs;
    }
    return(max_nbrs);
}

// A function to declare a matrix of given size and initialize it to 0
double** declare_matrix(int rows, int cols) {
    // Allocate memory for the matrix
    double** matrix = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
    }

    // Initialize the matrix to 0
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = 0;
        }
    }

    // Return the matrix
    return matrix;
}

// A function to fill up the non zero elements of a matrix
void fill_matrix(int** matrix, int rows, int cols) {
    // Loop through the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // If the element is 0, assign a random value between 1 and 10
            if (matrix[i][j] == 0) {
                matrix[i][j] = rand() % 10 + 1;
            }
        }
    }
}

// A function to print a matrix
void print_matrix(int** matrix, int rows, int cols) {
    // Loop through the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Print the element with a space
            printf("%d ", matrix[i][j]);
        }
        // Print a new line
        printf("\n");
    }
}

// A function to print a matrix
void print_matrix_d(double** matrix, int rows, int cols) {
    // Loop through the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Print the element with a space
            printf("%5f ", matrix[i][j]);
        }
        // Print a new line
        printf("\n");
    }
}

// A function to save a matrix in a file in CSV format
void save_matrix(double** matrix, int rows, int cols, char* filename) {
    // Open the file in write mode
    FILE* file = fopen(filename, "w");

    // Check if the file is opened successfully
    if (file == NULL) {
        printf("Error: could not open the file %s\n", filename);
        return;
    }

    // Loop through the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Write the element to the file with a comma
            fprintf(file, "%10.5f,", matrix[i][j]);
        }
        // Write a new line to the file
        fprintf(file, "\n");
    }

    // Close the file
    fclose(file);
}

// A function to calculate the number of set bits in a before each bit of b
int calculate (int a, int b) {
    // Initialize the count of set bits in a
    int count = popcnt (a) - 1;

    // Initialize the result array
    int* result = (int*) malloc (32 * sizeof (int));

    // Initialize the mask variable
    int m = 1;
    int phase = 0;

    // Initialize the index variable
    int i = 0;

    // Loop through each bit of b from right to left
    while (b) {
        // Append the current count to the result array
        phase += count;
        i++;
        printf(" %d (%d) ",b, m);

        // If the current bit of b is 1, then decrement the count by 1 if the corresponding bit of a is also 1
        if (b & 1) {
            if (a & m) {
                count--;
            }
        }

        // Left-shift the mask variable by 1
        m = m | (m << 1);

        // Right-shift b by 1
        b >>= 1;
    }
    printf("\n");

    // Return the result array
    return phase;
}

// Function to get the num  operator
void getNumOperator(size_t detIa, size_t detIb, double *numval, double *numvala, size_t* configAlpha, size_t sizeAlpha, size_t* configBeta, size_t sizeBeta, int nblk, size_t* NumBlock) {

    double n=0.0;
    size_t pos;
    // Loop over positions
    for(size_t p=0;p<nblk; ++p) {
      pos = NumBlock[p];
      // Check if alpha electron
      if( detIa & (1 << (pos)) ) {
        n += 1.0;
      }
      numvala[p] = n;
      
      // Check if beta electron
      if( detIb & (1 << (pos)) ) {
        n += 1.0;
      }
      numval[p] = n;
      n = 0.0;
    }
}

// Function to get the Sz  operator
void getSzOperator(size_t detIa, size_t detIb, double *szval, size_t* configAlpha, size_t sizeAlpha, size_t* configBeta, size_t sizeBeta, int nblk, size_t* SzBlock) {

    double sz=0.0;
    size_t pos;
    // Loop over positions
    for(size_t p=0;p<nblk; ++p) {
      pos = SzBlock[p];
      // Check if alpha electron
      if( detIa & (1 << (pos)) ) {
        sz += 0.5;
      }
      
      // Check if beta electron
      if( detIb & (1 << (pos)) ) {
        sz -= 0.5;
      }
      szval[p] = sz;
      sz = 0.0;
    }
}

