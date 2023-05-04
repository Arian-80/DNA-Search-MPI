#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

int getSequenceData(FILE* file, int* sequenceCount, int* maxSequenceLength) {
    // Get metadata about sequences
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
    /* This function has been taken from the MPI implementation of Bucketsort ..-
     * -.. within the same project. */
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
        // Finalize adding pattern
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
        sequences[i] = (char*) malloc((currentLength+1) * sizeof(char));
        if (!sequences[i]) {
            free_full_2DArray((void**) sequences, i); // Free all up to i
            return 0;
        }
        seqBuffer[currentLength] = '\0';
        strncpy(sequences[i], seqBuffer, (currentLength+1));
    }
    return 1;
}

int truncateSequences(char*** sequences, int portion, int start) {
    size_t truncatedSize = portion*sizeof(char*);
    memcpy(*sequences, &(*sequences)[start], truncatedSize);
    *sequences = (char**) realloc(*sequences, truncatedSize);
    if (!*sequences) {
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
                // End outer loop - factor whole loop into function later and return instead
                j = patternLength;
                break;
            }
            // index 0 is allowed because all instances of the first letter are the same
            commonStart[j] = patterns[0][j];
        }
    }
}

int patternMatchCommonStart(int patternCount, size_t patternLength,
                            const char* matchedString, char** patterns) {
    for (int j = 0; j < patternCount; j++) {
        for (int k = 0; k < patternLength; k++) {
            if (matchedString[k] != patterns[j][k]) {
                break;
            }
            if (k == patternLength-1) return 1; // Found
        }
    }
    return 0;
}

int findMatches(char* sequence, char** currentPatterns, char* commonStart,
                int** foundMatches, int* patternCount, int* currMatchCounter,
                int* totalMatches, int localIndex,
                size_t* currentSize, size_t patternLength) {
    char* temp = sequence;
    int found;
    while (temp[0] != '\0') { // End of sequence
        // Shortcut - initially only look for common start
        char* commonStartMatch = strstr(temp, commonStart);
        if (!commonStartMatch) break; // Not found
        found = patternMatchCommonStart(patternCount[localIndex],
                                        patternLength, commonStartMatch,
                                        currentPatterns);
        if (found) {
            (*foundMatches)[*totalMatches] = (int) (commonStartMatch - sequence);
            temp = &(temp[commonStartMatch - temp + 1]);
            (*totalMatches)++; (*currMatchCounter)++;
        }
        else temp = &temp[1]; // Skip one character to remove the already-found common start
        if (*totalMatches >= *currentSize) { // Enlarge array
            (*currentSize)*=2;
            *foundMatches = (int*) realloc(*foundMatches, (*currentSize)*sizeof(int));
            if (!*foundMatches) {
                return 0;
            }
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

int parallel_manageMatches(MPI_Comm communicator, char*** patterns,
                            char** sequence, int** foundMatches,
                            const size_t* patternLengths, int* totalMatches,
                            int* patternCount, int* matchCounter,
                            int start, size_t* currentSize) {
    // If parallel then portion is 1 sequence for each processor
    int rank, processorCount, localPortion, localStart, remainder;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &processorCount);
    int sequenceLength = (int) strlen(*sequence);
    getDistributions(&localStart, &localPortion, &remainder, sequenceLength,
                     processorCount, rank);

    size_t patternLength = patternLengths[start];
    // Start earlier in case split is at a matching case
    if (rank) {
        localStart -= (int) patternLength-1;
        localPortion += (int) patternLength-1;
    }
    memcpy(*sequence, &(*sequence)[localStart], localPortion*sizeof(char));
    *sequence = (char*) realloc(*sequence, (localPortion+1) * sizeof(char));
    if (!*sequence) return 0;
    (*sequence)[localPortion] = '\0';

    char** currentPatterns = patterns[start];
    char commonStart[patternLength];

    if (patternCount[start] == 1) strcpy(commonStart, currentPatterns[0]);
    else {
        findCommonStart(commonStart, currentPatterns, patternLength,
                        start, patternCount);
    }
    matchCounter[0] = 0;
    if (!findMatches(*sequence, currentPatterns, commonStart, foundMatches,
                     patternCount, &matchCounter[0], totalMatches, start,
                     currentSize, patternLength)) {
        free(*foundMatches);
        return 0;
    }
    // Index relative to full sequence.
    for (int i = 0; i < *totalMatches; i++) {
        (*foundMatches)[i] += localStart;
    }
    if (!rank) {
        int matchesFound[processorCount];
        MPI_Reduce(MPI_IN_PLACE, totalMatches, 1, MPI_INT, MPI_SUM, 0, communicator);
        // Need individual matches found, hence using gather instead of reduce.
        MPI_Gather(&matchCounter[0], 1, MPI_INT, matchesFound, 1, MPI_INT, 0,
                   communicator);
        matchCounter[0] = *totalMatches;
        *foundMatches = (int*) realloc(*foundMatches, matchCounter[0]*sizeof(int));
        if (!*foundMatches && matchCounter[0]) {
            return 0;
        }
        int* displs = (int*) malloc(processorCount * sizeof(int));
        if (!displs) {
            return 0;
        }
        getDispls(&displs, matchesFound, processorCount);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT, *foundMatches, matchesFound,
                    displs, MPI_INT, 0, communicator);
        free(displs);
    }
    else {
        MPI_Reduce(totalMatches, NULL, 1, MPI_INT, MPI_SUM, 0, communicator);
        MPI_Gather(&matchCounter[0], 1, MPI_INT, NULL, 1, MPI_INT, 0,
                   communicator);
        MPI_Gatherv(*foundMatches, matchCounter[0], MPI_INT, NULL, NULL,NULL,
                    MPI_INT, 0, communicator);
    }
    return 1;
}

int manageMatches(int parallel, int portion, int start, int sequenceCount,
                  int* totalMatches, int* patternCount, int* matchCounter,
                  int** foundMatches, const size_t* patternLengths,
                  char*** patterns, char** sequences,
                  MPI_Comm communicator) {
    char** currentPatterns;
    size_t patternLength;
    int iLocal;

    // Estimate there will be 5 matches in each sequence.
    size_t currentSize = sequenceCount * 5;

    (*foundMatches) = (int*) malloc(currentSize * sizeof(int));
    if (!*foundMatches) {
        return 0;
    }

    if (parallel) {
        if (!parallel_manageMatches(communicator, patterns, &sequences[0],
                                    foundMatches, patternLengths, totalMatches,
                                    patternCount, matchCounter, start,
                                    &currentSize)) {
            free(*foundMatches);
            return 0;
        }
        return 1;
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
        if (!findMatches(sequences[i], currentPatterns, commonStart, foundMatches,
                    patternCount, &matchCounter[i], totalMatches, iLocal,
                    &currentSize, patternLength)) {
            free(*foundMatches);
            return 0;
        }
    }
    return 1;
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
        if (processorCount > sequenceCount) remainder = 0;
        int noRemainder = remainder ? portion - 1 : portion;
        if (remainder) {
            for (int i = 0; i < remainder; i++) {
                recvcounts[i] = portion; // Rank 0 always accounts for remainder
            }
        }
        for (int i = remainder; i < processorCount; i++) {
            recvcounts[i] = noRemainder;
        }
        for (int i = sequenceCount; i < processorCount; i++) {
            // Processors with rank >= sequenceCount are auxiliary
            recvcounts[i] = 1;
        }

        int* displs = (int*) malloc(processorCount * sizeof(int));
        if (!displs) {
            free(*combinedMatchCount);
            free(*individualCounts);
            return 0;
        }
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
        for (int i = sequenceCount; i < processorCount; i++) {
            (*individualCounts)[i] = 0;
        }
        return 1;
    }
    else {
        int sendCount;
        if (rank >= sequenceCount) { // Auxiliary processors don't send
            sendCount = 1;
            matchCounter[0] = 0;
        }
        else {
            sendCount = portion;
        }
        MPI_Gatherv(matchCounter, sendCount, MPI_INT, NULL, NULL, NULL,
                    MPI_INT, 0, MPI_COMM_WORLD);
        return 1;
    }
}

int DNA_Search(FILE* file) {
    int sequenceCount = 0;
    int maxSequenceLength = 0;

    if (!getSequenceData(file, &sequenceCount, &maxSequenceLength)) {
        return 0;
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
        return 0;
    }


    // Rough estimation that all patterns include 3 variations.
    size_t estimatedPatternCount = 3;
    if (!malloc_pattern_array(patterns,estimatedPatternCount, sequenceCount)) {
        return 0;
    }

    int patternCount[sequenceCount]; // Number of patterns for each sequence
    // Separate out pattern variations
    if (!distinguishPatterns(file, patterns, patternLengths, sequenceCount,
                             estimatedPatternCount, patternCount, rank)) {
        return 0;
    }

    char** sequences = (char**) malloc(sequenceCount * sizeof(char*));
    if (!sequences) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        return 0;
    }
    // Fill "sequences" array with sequences found in the file.
    if (!parseSequences(file, sequences, sequenceCount, maxSequenceLength)) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        return 0;
    }

    // Load distribution
    if (!truncateSequences(&sequences, portion, start)) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        return 0;
    }

    int matchCounter[portion];
    int* foundMatches;
    int totalMatches = 0;

    int aux_processors = processorCount - sequenceCount;
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (aux_processors > 0) {
        MPI_Comm_split(MPI_COMM_WORLD, groupRank, rank, &communicator);
    }

    // Sort sequences in parallel
    if (!manageMatches(groupRank < aux_processors, portion, start, sequenceCount,
                       &totalMatches, patternCount, matchCounter,
                       &foundMatches, patternLengths, patterns,
                       sequences, communicator)) {
        free_pattern_arrays(patterns, patternLengths,sequenceCount,
                            sequenceCount, patternCount);
        free_full_2DArray((void**) sequences, portion);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OP);
        return 0;
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
        return 0;
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
//        int k = 0;
//        for (int i = 0; i < sequenceCount; i++) {
//            printf("Occurrences found in sequence %d: %d\n", i+1, combinedCount[i]);
//            for (int j = 0; j < combinedCount[i]; j++) {
//                printf("Occurrence %d: Index %d\n", j+1, combinedMatches[k]);
//                k++;
//            }
//        }
    }
    else {
        if (rank >= processorCount) totalMatches = 0;
        MPI_Gatherv(foundMatches, totalMatches, MPI_INT, NULL,
                    NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
    }
    free_pattern_arrays(patterns, patternLengths, sequenceCount,
                        sequenceCount, patternCount);
    free_full_2DArray((void**) sequences, portion);
    free(foundMatches);
    if (!rank) {
        free(combinedCount);
        free(combinedIndividualMatches);
        return 1;
    }
    else {
        return -1;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    FILE *file = fopen("test.txt", "r");
    if (!file) {
        printf("Input file not found.\n");
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_BUFFER);
        return 0;
    }
    int result;
    double start, end;
    start = MPI_Wtime();
    result = DNA_Search(file);
    end = MPI_Wtime();
    if (!result) {
        printf("An error has occurred.\n");
//        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OP);
        MPI_Finalize();
        return -1;
    } else if (result == -1) {
        MPI_Finalize();
        return 0;
    }
    printf("Time taken: %g\n", end - start);
    FILE *f = fopen("times.txt", "a");
    fprintf(f, "%g,", end - start);
    MPI_Finalize();
    return 0;
}
