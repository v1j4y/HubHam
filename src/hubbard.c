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

// Function to calculate the binomial coefficient (n choose k) using gamma and lgamma
unsigned long long binomialCoeff(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }

    double result = exp(lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1));
    return (unsigned long long)result;
}

// Function to calculate the total number of electrons in a configuration
int calculateTotalElectrons(const char *config) {
    int totalElectrons = 0;
    for (int i = 0; config[i] != '\0'; i++) {
        if (config[i] == '1') {
            totalElectrons++;
        } else if (config[i] == '2') {
            totalElectrons += 2; // Double occupancy
        }
    }
    return totalElectrons;
}

// Function to generate all configurations
void generateConfigurations(int norb, int nalpha, int nbeta, int currentOrbital, int currentAlpha, int currentBeta, char *config, char **configAll, int *count) {
    if (currentOrbital == norb || (currentAlpha + currentBeta) == (nalpha + nbeta)) {
        if (currentOrbital < norb) {
          for ( int i=currentOrbital; i<norb; ++i) {
            config[i] = '0';
          }
        }
        config[norb] = '\0'; // Null-terminate the configuration string
        int totalElectrons = calculateTotalElectrons(config);
        // Check if the total number of electrons in the configuration is within the limit
        if (totalElectrons == (nalpha + nbeta)) {
            strcpy(configAll[(*count)++], config);
        }
    } else {
        // Place a beta electron in the current orbital
        config[currentOrbital] = '0';
        generateConfigurations(norb, nalpha, nbeta, currentOrbital + 1, currentAlpha, currentBeta, config, configAll, count);

        if (currentAlpha + currentBeta < nalpha + nbeta) {
            // Allow for single occupancy in the current orbital
            config[currentOrbital] = '1';
            generateConfigurations(norb, nalpha, nbeta, currentOrbital + 1, currentAlpha + 1, currentBeta + 0, config, configAll, count);
        }
        config[currentOrbital] = '\0';

        if (currentAlpha + currentBeta < nalpha + nbeta) {
            // Allow for double occupancy in the current orbital
            config[currentOrbital] = '2';
            generateConfigurations(norb, nalpha, nbeta, currentOrbital + 1, currentAlpha + 1, currentBeta + 1, config, configAll, count);
        }
        config[currentOrbital] = '\0';
    }
}

int generateConfigurationsDriver(int norb, int nalpha, int nbeta) {

    unsigned long long numConfigurations = binomialCoeff(norb, nalpha) * binomialCoeff(norb, nbeta);

    char **configAll = (char **)malloc(numConfigurations * sizeof(char *));
    for (int i = 0; i < numConfigurations; i++) {
        configAll[i] = (char *)malloc(norb + 1);
    }

    char config[norb + 1]; // +1 for the null-terminator
    int count = 0;

    generateConfigurations(norb, nalpha, nbeta, 0, 0, 0, config, configAll, &count);

    printf("All configurations: %d \n", count);
    for (int i = 0; i < count; i++) {
        int totalElectrons = calculateTotalElectrons(configAll[i]);
        printf("%s (Total Electrons: %d)\n", configAll[i], totalElectrons);
        free(configAll[i]);
    }
    free(configAll);

    return 0;
}

