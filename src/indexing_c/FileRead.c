#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Reads the content of a text file and counts its number of lines.
 * @param filename The name of the file to read from.
 * @param content A pointer to the buffer where the file content will be stored.
 * @return The number of lines in the file.
 */
long readFileAndCountLines(const char* filename, char** content) {
    long flength;
    long numLines = 0;
    FILE* file = fopen(filename, "rb"); // Open the file in binary mode.

    if (!file) {
        printf("Error: Cannot open file %s\n", filename);
        exit(1);
    }

    // Get the file length.
    fseek(file, 0, SEEK_END);
    flength = ftell(file);
    fseek(file, 0, SEEK_SET);

    // Allocate memory for the buffer.
    *content = (char*)malloc(flength + 1);
    if (!*content) {
        printf("Error allocating memory\n");
        fclose(file);
        exit(1);
    }

    // Read the file content into the buffer.
    fread(*content, 1, flength, file);
    (*content)[flength] = '\0'; // Null-terminate the string.

    fclose(file); // Close the file.

    // Count the number of lines in the buffer.
    for (long i = 0; i < flength; i++) {
        if ((*content)[i] == '\n') {
            numLines++;
        }
    }

    // Check for the last line without a newline character.
    if ((*content)[flength - 1] != '\n') {
        numLines++;
    }

    return numLines;
}
