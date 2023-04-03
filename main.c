#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

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

void free2DArray(void** array, size_t length) {
    for (int i = 0; i < length; i++) {
        free(array[i]);
    }
    free(array);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    FILE *file;
    file = fopen("sequences.txt", "r");
    if (!file) {
        printf("File not found.\n");
        return -1;
    }

    int sequenceCount = 0;
    int maxSequenceLength = 0;
    fscanf(file, "%d", &sequenceCount);
    fscanf(file, "%d", &maxSequenceLength);

    if (!sequenceCount || sequenceCount < 0 ||
    !maxSequenceLength || maxSequenceLength < 0) {
        printf("Invalid input.\n");
        return -1;
    }

    char*** patterns = (char***) malloc(sequenceCount * sizeof(char**));
    if (!patterns) return -1;

    size_t** patternLengths = (size_t**) malloc(sequenceCount * sizeof(size_t*));
    if (!patternLengths) {
        free(patterns);
        return -1;
    }

    // Rough estimation that all patterns include 3 variations.
    size_t estimatedPatternCount = 3;
    for (int i = 0; i < sequenceCount; i++) {
        patterns[i] = (char **) malloc(estimatedPatternCount * sizeof(char *));
        if (!patterns[i]) {
            free2DArray((void**) patterns, i);
            free2DArray((void**) patternLengths, i);
            return -1;
        }
        patternLengths[i] = (size_t*) malloc(estimatedPatternCount * sizeof(size_t));
        if (!patternLengths[i]) {
            free2DArray((void**) patterns, i+1);
            free2DArray((void**) patternLengths, i);
            return -1;
        }
    }

    int patternCount[sequenceCount];

    char patternBuffer[255];
    size_t currentLength;
    int activeBracket = 0;
    int k;
    int currentVariations;
    char currentLetter;
    for (int i = 0; i < sequenceCount; i++){
        fscanf(file, "%s", patternBuffer);
        currentLength = strlen(patternBuffer);
        if (!currentLength) {
            for (int a = 0; a < i; a++) {
                free2DArray((void**) patterns[a], patternCount[a]);
            }
            free2DArray((void**) patterns, sequenceCount);
            free2DArray((void**) patternLengths, sequenceCount);
            return -1;
        }
        char withoutVariation[currentLength]; // Pattern before variable input
        estimatedPatternCount = 3;
        currentVariations = 0;
        k = 0;
        for (int j = 0; j < currentLength; j++) {
            currentLetter = patternBuffer[j];
            if (currentLetter == '[') {
                // Embedded brackets or more than one variation in the same pattern
                if (activeBracket || currentVariations) {
                    printf("Pattern %d is invalid.\n Embedded brackets are not"
                           " allowed and only one variation is permitted"
                           " within each pattern.\n", i);
                    for (int a = 0; a < i; a++) {
                        free2DArray((void**) patterns[a], patternCount[a]);
                    }
                    free2DArray((void**) patterns, sequenceCount);
                    free2DArray((void**) patternLengths, sequenceCount);
                    return -1;
                }
                activeBracket++;
                continue;
            }
            if (currentLetter == ']') {
                // No corresponding open bracket, or empty bracket
                if (!activeBracket || !currentVariations) {
                    printf("Pattern %d is invalid.\nClosing brackets "
                           "must have corresponding opening brackets and "
                           "brackets may not be empty.\n", i);
                    for (int a = 0; a < i; a++) {
                        free2DArray((void**) patterns[a], patternCount[a]);
                    }
                    free2DArray((void**) patterns, sequenceCount);
                    free2DArray((void**) patternLengths, sequenceCount);
                    return -1;
                }
                // Variation within current pattern already exists - limited to 1.
                activeBracket--;
                k = 0;
                continue;
            }
            if (activeBracket) {
                if (currentVariations == estimatedPatternCount) {
                    estimatedPatternCount *= 2;
                    patterns[i] = (char**) realloc(
                            patterns[i],estimatedPatternCount * sizeof(char*));
                    if (!patterns[i]) {
                        for (int a = 0; a < i; a++) {
                            free2DArray((void**) patterns[a], patternCount[a]);
                        }
                        free2DArray((void**) patterns, sequenceCount);
                        free2DArray((void**) patternLengths, sequenceCount);
                        return -1;
                    }
                    patternLengths[i] = (size_t*) realloc(patternLengths[i],
                                                          estimatedPatternCount * sizeof(size_t));
                    if (!patternLengths[i]) {
                        for (int a = 0; a < i; a++) {
                            free2DArray((void**) patterns[a], patternCount[a]);
                        }
                        free2DArray((void**) patterns, sequenceCount);
                        free2DArray((void**) patternLengths, sequenceCount);
                        return -1;
                    }
                }
                // currentLength is either the whole length or includes 2 brackets.
                patterns[i][currentVariations] = (char*) malloc((currentLength-2+1) * sizeof(char));
                if (!patterns[i][currentVariations]) {
                    for (int a = 0; a < i; a++) {
                        free2DArray((void**) patterns[a], patternCount[a]);
                    }
                    free2DArray((void**) patterns, sequenceCount);
                    free2DArray((void**) patternLengths, sequenceCount);
                    return -1;
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
            // hello[123]he hello1he hello2he hello3he
            // hello[1]a
        }
        withoutVariation[k] = '\0';
        if (!currentVariations) {
            patterns[i] = (char**) realloc(patterns[i], 1 * sizeof(char*));
            if (!patterns[i]) {
                for (int j = 0; j < i; j++) {
                    free2DArray((void**) patterns[j], patternCount[j]);
                }
                free2DArray((void**) patterns, sequenceCount);
                free2DArray((void**) patternLengths, sequenceCount);
                return -1;
            }

            patterns[i][0] = (char*) malloc((k+1) * sizeof(char));
            if (!patterns[i][0]) {
                for (int j = 0; j < i; j++) {
                    free2DArray((void**) patterns[j], patternCount[j]);
                }
                free2DArray((void**) patterns, sequenceCount);
                free2DArray((void**) patternLengths, sequenceCount);
                return -1;
            }
            strncpy(patterns[i][0], withoutVariation, k+1);
            patternLengths[i][0] = k;
            patternCount[i] = 1;
            continue;
        }
        currentLength = 0;
        for (int j = 0; j < currentVariations; j++) {
            strncat(patterns[i][j], withoutVariation, k+1);
            if (!currentLength) { // All patterns variations are of the same length
                currentLength = strlen(patterns[i][j]);
            }
            patternLengths[i][j] = currentLength;
            patterns[i][j] = (char*) realloc(patterns[i][j],
                                  (currentLength + 1) * sizeof(char));
            if (!patterns[i][j]) {
                for (k = 0; k < i; k++) {
                    free2DArray((void**) patterns[k], patternCount[k]);
                }
                free2DArray((void**) patterns, sequenceCount);
                free2DArray((void**) patternLengths, sequenceCount);
                return -1;
            }
        }
        patternCount[i] = currentVariations;
    }

    char** sequences = (char**) malloc(sequenceCount * sizeof(char*));
    if (!sequences) {
        for (int i = 0; i < sequenceCount; k++) {
            free2DArray((void**) patterns[i], patternCount[i]);
        }
        free2DArray((void**) patterns, sequenceCount);
        free2DArray((void**) patternLengths, sequenceCount);
        return -1;
    }

    size_t sequenceLengths[sequenceCount];
    char seqBuffer[maxSequenceLength+1];
    for (int i = 0; i < sequenceCount; i++) {
        fscanf(file, "%s", seqBuffer);
        currentLength = strlen(seqBuffer);
        sequenceLengths[i] = currentLength;
        sequences[i] = (char*) malloc(currentLength * sizeof(char));
        if (!sequences[i]) {
            for (int j = 0; j < sequenceCount; j++) {
                free2DArray((void**) patterns[j], patternCount[j]);
            }
            free2DArray((void**) patterns, sequenceCount);
            free2DArray((void**) patternLengths, sequenceCount);
            free2DArray((void**) sequences, i); // Free all up to i
            return -1;
        }
        seqBuffer[currentLength] = '\0';
        strncpy(sequences[i], seqBuffer, currentLength+1);
    }

    /* Testing
    for (int i = 0; i < sequenceCount; i++) {
        printf("i: %d\tSequence: %s\tPatternCount: %d\n", i, sequences[i], patternCount[i]);
        for (int j = 0; j < patternCount[i]; j++) {
            printf("Pattern: %s\tPatternLength: %d\n", patterns[i][j], patternLengths[i][j]);
        }
    }
    return 0;
     */

    int rank, processorCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processorCount);

    // Processors work on different or the same sequence based on their ranks
    int groupRank = rank % sequenceCount;

    int start, portion, remainder;
    getDistributions(&start, &portion, &remainder, sequenceCount,
                     processorCount, groupRank);

    size_t truncatedSize = portion*sizeof(char*);
    memcpy(&sequences[0], &sequences[start], truncatedSize);
    sequences = (char**) realloc(sequences, truncatedSize);
    if (!sequences) {
        for (int j = 0; j < sequenceCount; j++) {
            free2DArray((void**) patterns[j], patternCount[j]);
        }
        free2DArray((void**) patterns, sequenceCount);
        free2DArray((void**) patternLengths, sequenceCount);
        return -1;
    }

    size_t currPatternLength;
    size_t commonStartLength;
    size_t sequenceLength;
    int iLocal;
    char** currentPatterns;
    int** foundMatches = (int**) malloc(portion * sizeof(int*));
    int matchCounter[portion];
    for (int i = 0; i < portion; i++) {
        iLocal = i+start;
        sequenceLength = sequenceLengths[iLocal];
        currentPatterns = patterns[iLocal];
        currPatternLength = patternLengths[iLocal][0];
        // Estimating (x/y)/2 occurrences will be found, where x is the sequence length and y is the pattern length.
        size_t estimatedOccurrences = sequenceLength/currPatternLength/2;
        foundMatches[i] = (int*) malloc(estimatedOccurrences * sizeof(int));
//        CHECK TO MAKE SURE MALLOC ABOVE² AND REALLOC BELOW WORK, OTHERWISE FREE EVEYRTHING ABOVE... :)
        char* commonStart = (char*) malloc(currPatternLength * sizeof(char));
        if (patternCount[iLocal] == 1) commonStart = currentPatterns[0];
        else {
            for (int j = 0; j < currPatternLength; j++) {
                for (int a = 1; a < patternCount[iLocal]; a++) {
                    if (currentPatterns[a][j] != currentPatterns[a - 1][j]) {
                        commonStart[j] = '\0';
                        j = currPatternLength; // End outer loop - factor whole loop into function later and return instead
                        break;
                    }
                }
                commonStart[j] = currentPatterns[0][j]; // index 0 is allowed because all instances of the first letter are the same per the above loop
            }
        }
        matchCounter[i] = 0;
        commonStartLength = strlen(commonStart);
        char* temp = sequences[i];
        while (temp[0] != '\0') {
            char* currentMatch = strstr(temp, commonStart);
            if (!currentMatch) break; // Not found
            if (matchCounter[i] >= estimatedOccurrences) {
                estimatedOccurrences*=2;
                foundMatches[i] = realloc(foundMatches[i], estimatedOccurrences*sizeof(int));
            }
            foundMatches[i][matchCounter[i]] = (int) (currentMatch - sequences[i]); // Index of first found instance
            temp = &(sequences[i][foundMatches[i][matchCounter[i]]+commonStartLength]); // Skip the matched part
            matchCounter[i]++;
        }
    }


    // Consider getting rid of sequenceLength
    // FREE COMMONSTART
    for (int i = 0; i < portion; i++) {
        for (int j = 0; j < matchCounter[i]; j++) {
            printf("Rank: %d, match[%d]: %d\n", rank, j, foundMatches[i][j]);
        }
    }

    for (int j = 0; j < sequenceCount; j++) {
        free2DArray((void**) patterns[j], patternCount[j]);
    }
    free(patterns);
    free2DArray((void**) patternLengths, sequenceCount);
    free2DArray((void**) sequences, portion);

    MPI_Finalize();
    return 0;
}
