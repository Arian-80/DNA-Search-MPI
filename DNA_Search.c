#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

int getSequenceData(FILE* file, int* sequenceCount, int* maxSequenceLength) {
    fscanf(file, "%d", sequenceCount);
    fscanf(file, "%d", maxSequenceLength);

    if (!sequenceCount || sequenceCount < 0 ||
        !maxSequenceLength || maxSequenceLength < 0) {
        printf("Invalid input.\n");
        return 0;
    }
    return 1;
}

void getDistributions(int *start, int *portion, int *remainder, int size,
                      int processorCount, int rank) {
    (*portion) = size/processorCount;
    (*remainder) = size % processorCount;
    if (rank < (*remainder)) { // Spread remainder evenly across processors
        (*portion)++; // Processors get a max 1 remainder bucket to deal with
        *start = (*portion)*rank; // All processors before have had remainder buckets
    }
    else {
        *start = ((*portion)+1)*(*remainder) + (*portion)*(rank-(*remainder));
    }
}

int malloc_pattern_arrays(char**** patterns, size_t** patternLengths,
                          int elementCount) {
    *patterns = (char***) malloc(elementCount * sizeof(char**));
    if (*patterns) {
        *patternLengths = (size_t *) malloc(elementCount * sizeof(size_t));
        if (!*patternLengths) {
            free(*patterns);
            return 0;
        }
        return 1;
    }
    return 0;
}

void free_inner_2DArray(void** array, size_t length) {
    for (int i = 0; i < length; i++) {
        free(array[i]);
    }
}

void free_full_2DArray(void** array, int length) {
    free_inner_2DArray(array, length);
    free(array);
}


int malloc_pattern_array(char*** array1, size_t elementCount, int repeat) {
    // The function parameters are specific for less complexity as this is the only purpose they are used for.
    for (int i = 0; i < repeat; i++) {
        array1[i] = (char**) malloc(elementCount * sizeof(char *));
        if (!array1[i]) {
            free_full_2DArray((void **) array1, i);
            return 0;
        }
    }
    return 1;
}

void free_pattern_arrays(char*** patterns, size_t* patternLengths,
                         int repeats, int size, int* patternCount) {
    for (int i = 0; i < repeats; i++) {
        free_inner_2DArray((void**) patterns[i], patternCount[i]);
    }
    free_full_2DArray((void**) patterns, size);
    free(patternLengths);
}

FILE* getFile(char* fileName) {
    FILE* file = fopen(fileName, "r");
    return file;
}

int distinguishPatterns(FILE* file, char*** patterns, size_t* patternLengths,
                        int sequenceCount, size_t estimatedPatternCount,
                        int* patternCount, int rank) {
    char patternBuffer[255]; // Max length of a pattern
    int activeBracket = 0;
    int k;
    int currentVariations;
    size_t currentLength;
    char currentLetter;
    for (int i = 0; i < sequenceCount; i++){
        fscanf(file, "%s", patternBuffer);
        currentLength = strlen(patternBuffer);
        char withoutVariation[currentLength]; // Pattern before variable input
        currentVariations = 0;
        k = 0;
        for (int j = 0; j < currentLength; j++) {
            currentLetter = patternBuffer[j];
            if (currentLetter == '[') {
                // Embedded brackets or more than one variation in the same pattern
                if (activeBracket || currentVariations) {
                    if (!rank) {
                        printf("Pattern %d is invalid.\nEmbedded brackets are not"
                               " allowed and only one variation is permitted"
                               " within each pattern.\n", i);
                    }
                    free_pattern_arrays(patterns, patternLengths, i,
                                        sequenceCount, patternCount);
                    return 0;
                }
                withoutVariation[k] = '\0';
                activeBracket++;
                continue;
            }
            if (currentLetter == ']') {
                // No corresponding open bracket, or empty bracket
                if (!activeBracket || !currentVariations) {
                    if (!rank) {
                        printf("Pattern %d is invalid.\nClosing brackets "
                               "must have corresponding opening brackets and "
                               "brackets may not be empty.\n", i);
                    }
                    free_pattern_arrays(patterns, patternLengths, i,
                                        sequenceCount, patternCount);
                    return 0;
                }
                // Variation within current pattern already exists - limited to 1.
                activeBracket--;
                k = 0;
                continue;
            }
            if (activeBracket) {
                if (currentVariations == estimatedPatternCount) {
                    patterns[i] = (char**) realloc(
                            patterns[i],currentVariations*2 * sizeof(char*));

                    if (!patterns[i] || patterns[i][currentVariations]) {
                        free_pattern_arrays(patterns, patternLengths, i,
                                            sequenceCount, patternCount);
                        return 0;
                    }
                }
                // currentLength is either the whole length or includes 2 brackets.
                patterns[i][currentVariations] = (char*)
                        malloc((currentLength-2+1) * sizeof(char));
                if (!patterns[i][currentVariations]) {
                    free_pattern_arrays(patterns, patternLengths, i,
                                        sequenceCount, patternCount);
                    return 0;
                }
                // Store one variation of the pattern
                strncpy(patterns[i][currentVariations], withoutVariation, k);
                patterns[i][currentVariations][k] = '\0';
                strncat(patterns[i][currentVariations], &currentLetter, 1);
                patternCount[i]++;
                currentVariations++;
            }
            else {
                withoutVariation[k] = currentLetter;
                k++;
            }
        }
        withoutVariation[k] = '\0';
        if (!currentVariations) {
            patterns[i] = (char**) realloc(patterns[i], 1 * sizeof(char*));
            if (patterns[i]) {
                patterns[i][0] = (char*) malloc((k+1) * sizeof(char));
            }
            if (!patterns[i] || !patterns[i][0]) {
                free_pattern_arrays(patterns, patternLengths, i,
                                    sequenceCount, patternCount);
                return 0;
            }
            strncpy(patterns[i][0], withoutVariation, k+1);
            patternLengths[i] = k;
            patternCount[i] = 1;
            continue;
        }
        currentLength = 0;
        for (int j = 0; j < currentVariations; j++) {
            strncat(patterns[i][j], withoutVariation, k+1);
            if (!currentLength) { // All patterns variations are of the same length
                currentLength = strlen(patterns[i][j]);
            }
            patternLengths[i] = currentLength;
            patterns[i][j] = (char*) realloc(patterns[i][j],
                                             (currentLength + 1) * sizeof(char));
            if (!patterns[i][j]) {
                free_pattern_arrays(patterns, patternLengths, i,
                                    sequenceCount, patternCount);
                return 0;
            }
        }
        patternCount[i] = currentVariations;
    }
    return 1;
}

int parseSequences(FILE* file, char** sequences, int sequenceCount, int maxSequenceLength) {
    char seqBuffer[maxSequenceLength+1];
    size_t currentLength;
    for (int i = 0; i < sequenceCount; i++) {
        fscanf(file, "%s", seqBuffer);
        currentLength = strlen(seqBuffer);
        sequences[i] = (char*) malloc(currentLength * sizeof(char));
        if (!sequences[i]) {
            free_full_2DArray((void**) sequences, i); // Free all up to i
            return 0;
        }
        seqBuffer[currentLength] = '\0';
        strncpy(sequences[i], seqBuffer, currentLength+1);
    }
    return 1;
}

int truncateSequences(char** sequences, int portion, int start) {
    size_t truncatedSize = portion*sizeof(char*);
    memcpy(&sequences[0], &sequences[start], truncatedSize);
    sequences = (char**) realloc(sequences, truncatedSize);
    if (!sequences) {
        return 0;
    }
    return 1;
}

void findCommonStart(char* commonStart, char** patterns, size_t patternLength,
                     int index, const int* patternCount) {
    for (int j = 0; j < patternLength; j++) {
        for (int k = 1; k < patternCount[index]; k++) {
            if (patterns[k][j] != patterns[k - 1][j]) {
                commonStart[j] = '\0';
                j = patternLength; // End outer loop - factor whole loop into function later and return instead
                break;
            }
            commonStart[j] = patterns[0][j]; // index 0 is allowed because all instances of the first letter are the same per the above loop
        }
    }
}

void calculateEstimatedOccurrence(char** sequences, int index,
                                  size_t patternLength, int* estimatedOccurrences) {
    size_t sequenceLength = strlen(sequences[index]);
    // Estimating (x/y)/2 occurrences will be found, where x is the sequence length and y is the pattern length.
    (*estimatedOccurrences) = sequenceLength / patternLength / 2;
    if (!(*estimatedOccurrences)) (*estimatedOccurrences) = 1; // Estimate at least one occurrence
}

int patternMatchCommonStart(int patternCount, size_t commonStartLength,
                            size_t patternLength, const char* matchedString,
                            char** patterns) {
    for (int j = 0; j < patternCount; j++) {
        for (int k = commonStartLength; k < patternLength; k++) {
            if (matchedString[k] != patterns[j][k]) {
                break;
            }
            if (k == patternLength-1) return 1; // Found
        }
    }
    return 0;
}

int findMatches(char** sequences, char** currentPatterns, char* commonStart,
                int** foundMatches, int* patternCount, int* currMatchCounter,
                int* totalMatches, int currIndex, int localIndex,
                size_t* currentSize, size_t patternLength) {
    size_t commonStartLength = strlen(commonStart);
    char* temp = sequences[currIndex];
    int found;
    while (temp[0] != '\0') { // End of sequence
        char* commonStartMatch = strstr(temp, commonStart); // Shortcut - initially only look for common start
        if (!commonStartMatch) break; // Not found
        found = patternMatchCommonStart(patternCount[localIndex], commonStartLength,
                                        patternLength, commonStartMatch,
                                        currentPatterns);
        if (found) {
            (*foundMatches)[*currMatchCounter] = (int) (commonStartMatch - sequences[currIndex]);
            temp = &(temp[commonStartMatch - temp + 1]);
            (*currMatchCounter)++;
        }
        else temp = &temp[1]; // Skip one character to remove the already-found common start
        if (*currMatchCounter >= *currentSize) { // Enlarge array
            (*currentSize)*=2;
            *foundMatches = (int*) realloc(*foundMatches, (*currentSize)*sizeof(int));
            if (!*foundMatches) {
                return 0;
            }
        }
    }
    (*totalMatches) += (*currMatchCounter);
    return 1;
}

int manageMatches(int portion, int start, int sequenceCount, int* totalMatches,
                  int* patternCount, int* matchCounter,  int** foundMatches,
                  const size_t* patternLengths, char*** patterns,
                  char** sequences) {
    char** currentPatterns;
    size_t patternLength;
    int iLocal;

    // Estimate there will be 3 matches in each sequence.
    size_t currentSize = sequenceCount * 3;

    (*foundMatches) = (int*) malloc(currentSize * sizeof(int));
    if (!*foundMatches) {
        return 0;
    }

    for (int i = 0; i < portion; i++) {
        iLocal = i+start;
        patternLength = patternLengths[iLocal];
        currentPatterns = patterns[iLocal];
        char commonStart[patternLength];
        if (patternCount[iLocal] == 1) strcpy(commonStart, currentPatterns[0]);
        else {
            findCommonStart(commonStart, currentPatterns, patternLength,
                            iLocal, patternCount);
        }
        matchCounter[i] = 0;
        if (!findMatches(sequences, currentPatterns, commonStart, foundMatches,
                    patternCount, &matchCounter[i], totalMatches, i, iLocal,
                    &currentSize, patternLength)) {
            free(*foundMatches);
            return 0;
        }
    }
    return 1;
}

void getDispls(int** displs, const int* recvcounts, int processorCount) {
    // Work out displacements
    (*displs)[0] = 0;
    for (int i = 1; i < processorCount; i++) {
        (*displs)[i] = recvcounts[i - 1] + (*displs)[i - 1];
    }
}

int getCombinedCounts(int** combinedMatchCount, int** individualCounts,
                      int* matchCounter, int rank, int portion,
                      int remainder, int sequenceCount,
                      int processorCount) {
    if (!rank) {
        *combinedMatchCount = (int*) calloc(processorCount, sizeof(int));
        if (*combinedMatchCount) {
            *individualCounts = (int*) calloc(processorCount, sizeof(int));
            if (!*individualCounts) {
                free(*combinedMatchCount);
                return 0;
            }
        }
        else {
            return 0;
        }

        int recvcounts[processorCount];
        int noRemainder = remainder ? portion - 1 : portion;

        if (remainder) {
            for (int i = 0; i < remainder; i++) {
                recvcounts[i] = portion; // Rank 0 always accounts for remainder
            }
        }
        for (int i = remainder; i < sequenceCount; i++) {
            recvcounts[i] = noRemainder;
        }
        for (int i = sequenceCount; i < processorCount; i++) {
            // Processors with rank > sequenceCount are auxiliary
            recvcounts[i] = 0;
        }

        int* displs = (int*) malloc(processorCount * sizeof(int));
        getDispls(&displs, recvcounts, processorCount);
        MPI_Gatherv(matchCounter, portion, MPI_INT, (*combinedMatchCount),
                    recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
        free(displs);
        int k = 0;
        for (int i = 0; i < processorCount; i++) {
            for (int j = 0; j < recvcounts[i]; j++) {
                (*individualCounts)[i] += (*combinedMatchCount)[k];
                k++;
            }
        }
        return 1;
    }
    else {
        int sendCount;
        if (rank >= sequenceCount) { // Auxiliary processors don't send
            sendCount = 0;
        }
        else {
            sendCount = portion;
        }
        MPI_Gatherv(matchCounter, sendCount, MPI_INT, NULL, NULL, NULL,
                    MPI_INT, 0, MPI_COMM_WORLD);
        return 1;
    }
}


int main(int argc, char** argv) {
    FILE *file = getFile("sequences.txt");
    if (!file) {
        printf("File not found.\n");
        return -1;
    }

    MPI_Init(&argc, &argv);

    int sequenceCount = 0;
    int maxSequenceLength = 0;

    if (!getSequenceData(file, &sequenceCount, &maxSequenceLength)) {
        MPI_Finalize();
        return -1;
    }


    int rank, processorCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processorCount);

    // Processors work on different or the same sequence based on their ranks
    int groupRank = rank % sequenceCount;

    int start, portion, remainder;
    getDistributions(&start, &portion, &remainder, sequenceCount,
                     processorCount, groupRank);


    char*** patterns;
    size_t* patternLengths;
    if (!malloc_pattern_arrays(&patterns, &patternLengths, sequenceCount)) {
        MPI_Finalize();
        return -1;
    }


    // Rough estimation that all patterns include 3 variations.
    size_t estimatedPatternCount = 3;
    if (!malloc_pattern_array(patterns,estimatedPatternCount, sequenceCount)) {
        MPI_Finalize();
        return -1;
    }

    int patternCount[sequenceCount]; // Number of patterns for each sequence
    // Separate out pattern variations
    if (!distinguishPatterns(file, patterns, patternLengths, sequenceCount,
                             estimatedPatternCount, patternCount, rank)) {
        MPI_Finalize();
        return -1;
    }

    char** sequences = (char**) malloc(sequenceCount * sizeof(char*));
    if (!sequences) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        MPI_Finalize();
        return -1;
    }

    // Fill "sequences" array with sequences found in the file.
    if (!parseSequences(file, sequences, sequenceCount, maxSequenceLength)) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        MPI_Finalize();
        return -1;
    }

    if (!truncateSequences(sequences, portion, start)) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        MPI_Finalize();
        return -1;
    }


    int matchCounter[portion];
    int* foundMatches;
    int totalMatches = 0;

    if (!manageMatches(portion, start, sequenceCount, &totalMatches,
                       patternCount, matchCounter,&foundMatches,
                       patternLengths, patterns, sequences)) {
        free_pattern_arrays(patterns, patternLengths,sequenceCount,
                            sequenceCount, patternCount);
        free_full_2DArray((void**) sequences, portion);
        MPI_Finalize();
        return -1;
    }

    int *combinedCount;
    int *combinedIndividualMatches;
    if (!getCombinedCounts(&combinedCount, &combinedIndividualMatches, matchCounter,
                           rank, portion, remainder, sequenceCount,
                           processorCount)) {
        free_pattern_arrays(patterns, patternLengths,sequenceCount,
                            sequenceCount, patternCount);
        free_full_2DArray((void**) sequences, portion);
        free(foundMatches);
        MPI_Finalize();
        return -1;
    }

    if (!rank) {
        int totalMatchesFound = 0;
        for (int i = 0; i < sequenceCount; i++) {
            totalMatchesFound += combinedCount[i];
        }
        int combinedMatches[totalMatchesFound];
        int* displs = (int*) malloc(processorCount * sizeof(int));
        getDispls(&displs, combinedIndividualMatches, processorCount);
        MPI_Gatherv(foundMatches, totalMatches, MPI_INT, combinedMatches,
                    combinedIndividualMatches, displs, MPI_INT, 0, MPI_COMM_WORLD);
        free(displs);
        int k = 0;
        for (int i = 0; i < sequenceCount; i++) {
            printf("Occurrences found in sequence %d: %d\n", i, combinedCount[i]);
            for (int j = 0; j < combinedCount[i]; j++) {
                printf("Occurrence %d: Index %d\n", j+1, combinedMatches[k]);
                k++;
            }
        }
    }
    else {
        MPI_Gatherv(foundMatches, totalMatches, MPI_INT, NULL,
                    NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
    }

//    int k = 0;
//    for (int i = 0; i < portion; i++) {
//        for (int j = 0; j < matchCounter[i]; j++) {
//            printf("Match[%d]: %d\t", j, foundMatches[k]);
//            k++;
//        }
//        printf("\n");
//    }

    free_pattern_arrays(patterns, patternLengths, sequenceCount,
                        sequenceCount, patternCount);
    free_full_2DArray((void**) sequences, portion);
    free(foundMatches);
    if (!rank) {
        free(combinedCount);
        free(combinedIndividualMatches);
    }
    MPI_Finalize();
    return 0;
}
