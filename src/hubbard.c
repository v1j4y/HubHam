#include "hubbard.h"

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

