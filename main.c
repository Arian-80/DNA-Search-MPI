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

void free_inner_2DArray(void** array, size_t length) {
    for (int i = 0; i < length; i++) {
        free(array[i]);
    }
}

void free_full_2DArray(void** array, int length) {
    free_inner_2DArray(array, length);
    free(array);
}


int malloc_pattern_arrays(char*** array1, size_t** array2, size_t elementCount, int repeat) {
    // The function parameters are specific for less complexity as this is the only purpose they are used for.
    for (int i = 0; i < repeat; i++) {
        array1[i] = (char **) malloc(elementCount * sizeof(char *));
        if (!array1[i]) {
            free_full_2DArray((void**) array1, i);
            free_full_2DArray((void**) array2, i);
            return 0;
        }
        array2[i] = (size_t*) malloc(elementCount * sizeof(size_t));
        if (!array2[i]) {
            free_full_2DArray((void**) array1, i+1);
            free_full_2DArray((void**) array2, i);
            return 0;
        }
    }
    return 1;
}

void free_pattern_arrays(char*** patterns, size_t** patternLengths,
                         int repeats, int size, int* patternCount) {
    for (int i = 0; i < repeats; i++) {
        free_inner_2DArray((void**) patterns[i], patternCount[i]);
    }
    free_full_2DArray((void**) patterns, size);
    free_full_2DArray((void**) patternLengths, size);
}

FILE* getFile(char* fileName) {
    FILE* file = fopen(fileName, "r");
    return file;
}

int distinguishPatterns(FILE* file, char*** patterns, size_t** patternLengths,
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

                    patternLengths[i] = (size_t*) realloc(
                            patternLengths[i],currentVariations*2 * sizeof(size_t));

                    if (!patterns[i] || !patternLengths[i] ||
                        patterns[i][currentVariations]) {
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


int main(int argc, char** argv) {
    FILE *file = getFile("sequences.txt");
    if (!file) {
        printf("File not found.\n");
        return -1;
    }

    MPI_Init(&argc, &argv);

    int sequenceCount = 0;
    int maxSequenceLength = 0;
    fscanf(file, "%d", &sequenceCount);
    fscanf(file, "%d", &maxSequenceLength);

    if (!sequenceCount || sequenceCount < 0 ||
    !maxSequenceLength || maxSequenceLength < 0) {
        printf("Invalid input.\n");
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


    char*** patterns = (char***) malloc(sequenceCount * sizeof(char**));
    if (!patterns) return -1;

    size_t** patternLengths = (size_t**) malloc(sequenceCount * sizeof(size_t*));
    if (!patternLengths) {
        free(patterns);
        return -1;
    }

    // Rough estimation that all patterns include 3 variations.
    size_t estimatedPatternCount = 3;
    if (!malloc_pattern_arrays(patterns, patternLengths,
                   estimatedPatternCount, sequenceCount)) {
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

    size_t patternLength;
    size_t commonStartLength;
    int iLocal;
    char** currentPatterns;
    int** foundMatches = (int**) malloc(portion * sizeof(int*));
    if (!foundMatches) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        MPI_Finalize();
        return -1;
    }
    int matchCounter[portion];
    for (int i = 0; i < portion; i++) {
        iLocal = i+start;
        patternLength = patternLengths[iLocal][0];
        currentPatterns = patterns[iLocal];

        int estimatedOccurrences;
        calculateEstimatedOccurrence(sequences, i,
                                     patternLength,&estimatedOccurrences);

        foundMatches[i] = (int*) malloc(estimatedOccurrences * sizeof(int));
        char commonStart[patternLength];
        if (!foundMatches[i]) {
            free_pattern_arrays(patterns, patternLengths,
                                sequenceCount, sequenceCount, patternCount);
            free_full_2DArray((void**) foundMatches, i);
            MPI_Finalize();
            return -1;
        }

        if (patternCount[iLocal] == 1) strcpy(commonStart, currentPatterns[0]);
        else {
            findCommonStart(commonStart, currentPatterns, patternLength, iLocal, patternCount);
        }

        matchCounter[i] = 0;
        commonStartLength = strlen(commonStart);
        char* temp = sequences[i];
        int found;
        while (temp[0] != '\0') { // End of sequence
            char* commonStartMatch = strstr(temp, commonStart); // Shortcut - initially only look for common start
            if (!commonStartMatch) break; // Not found
            found = patternMatchCommonStart(patternCount[iLocal], commonStartLength,
                                            patternLength, commonStartMatch,
                                            currentPatterns);
            if (found) {
                foundMatches[i][matchCounter[i]] = (int) (commonStartMatch - sequences[i]);
                temp = &(temp[commonStartMatch - temp + patternLength]);
                matchCounter[i]++;
            }
            else temp = &temp[1]; // Skip one character to remove the already-found common start
            if (matchCounter[i] >= estimatedOccurrences) { // Enlarge array
                estimatedOccurrences*=2;
                foundMatches[i] = realloc(foundMatches[i], estimatedOccurrences*sizeof(int));
                if (!foundMatches) {
                    free_pattern_arrays(patterns, patternLengths,sequenceCount,
                                        sequenceCount, patternCount);
                    free_full_2DArray((void**) sequences, sequenceCount);
                    free(foundMatches);
                    MPI_Finalize();
                    return -1;
                }
            }
        }
    }


    // Consider making patternLengths a 1d array instead of 2d.
    // Patterns for the same sequences are of the same length always, so 1d array of length 'portion' suffices.
    for (int i = 0; i < portion; i++) {
        for (int j = 0; j < matchCounter[i]; j++) {
            printf("Rank: %d, match[%d]: %d\n", rank, j, foundMatches[i][j]);
        }
    }
    free_pattern_arrays(patterns, patternLengths, sequenceCount,
                        sequenceCount, patternCount);
    free_full_2DArray((void**) sequences, portion);

    MPI_Finalize();
    return 0;
}
