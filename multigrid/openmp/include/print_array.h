#ifndef PRINT_ARRAY
#define PRINT_ARRAY

#include "array.h"
#include "parser.h"

static int
is_little_endian(void) {
    int num = 1;
    return (*((char *)&num) == 1);
}

void print_vtk(Array & array, Settings & settings, const char *fname) {

    FILE *f_ptr;
    uint_t written = 0;
    uint_t items = array.size;
    uint_t i,j,k;

    if ( (f_ptr = fopen(fname, "w")) == NULL ) {
       perror("No output! fopen()");
       return;
    }

    // Write VTK file header
    fprintf(f_ptr, "# vtk DataFile Version 3.0\n");
    fprintf(f_ptr, "saved from function print_vtk.\n");
    fprintf(f_ptr, "BINARY\n");
    fprintf(f_ptr, "DATASET STRUCTURED_POINTS\n");
    fprintf(f_ptr, "DIMENSIONS %d %d %d\n", array.shape[0], array.shape[1], array.shape[2]);
    fprintf(f_ptr, "ORIGIN %d %d %d\n", settings.origin[0],settings.origin[1],settings.origin[2]);
    fprintf(f_ptr, "SPACING %d %d %d\n", settings.h,settings.h,settings.h);
    fprintf(f_ptr, "POINT_DATA %lu\n", items);
    fprintf(f_ptr, "SCALARS %s %s 1\n", "gray", "double");
    fprintf(f_ptr, "LOOKUP_TABLE default\n");

    if ( is_little_endian() ) {
        // System is little endian, so we need to reverse the byte order.
        for (k = 0; k < array.shape[0]; ++k) {
            for (j = 0; j < array.shape[1]; ++j) {
                for (i = 0; i < array.shape[2]; ++i) {
            uint64_t crnt = *(uint64_t *)(array.at[array.idx(i,j,k)]); // Get double as int

            // Reverse byte order and write to file
            crnt = (crnt & 0x00000000FFFFFFFF) << 32 | (crnt & 0xFFFFFFFF00000000) >> 32;
            crnt = (crnt & 0x0000FFFF0000FFFF) << 16 | (crnt & 0xFFFF0000FFFF0000) >> 16;
            crnt = (crnt & 0x00FF00FF00FF00FF) << 8  | (crnt & 0xFF00FF00FF00FF00) >> 8;
            written += fwrite(&crnt, sizeof(uint64_t), 1, f_ptr);
                }
            }
        }
    } else {
        // System is big endian, so just dump the data.
        for (k = 0; k < array.shape[0]; ++k) {
            for (j = 0; j < array.shape[1]; ++j) {
                for (i = 0; i < array.shape[2]; ++i) {
                    written += fwrite(array.at[array.idx(i,j,k)], sizeof(double), 1, f_ptr);
                }
            }
        }
    }

    if ( written != items ) {
	    fprintf(stderr, "Writing failed:  only %lu of %lu items saved!\n",
		written, items);
    }

    fclose(f_ptr);
}

#endif