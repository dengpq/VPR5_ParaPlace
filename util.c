#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "globals.h"
#include "const.h"

/* This file contains utility functions widely used in *
 * my programs.  Many are simply versions of file and  *
 * memory grabbing routines that take the same         *
 * arguments as the standard library ones, but exit    *
 * the program if they find an error condition.        */

int linenum;  /* Line in file being parsed. */

/* Returns the min of cur and max. If cur > max, a warning
 * is emitted. */
int limit_value(int cur,
                int max,
                const char* name)
{
    if (cur > max) {
        printf(WARNTAG "%s is being limited from [%d] to [%d]\n",
               name, cur, max);
        return max;
    }

    return cur;
}

/* An alternate for strncpy since strncpy doesn't work as most
 * people would expect. This ensures null termination */
char* my_strncpy(char* dest,
                 const char* src,
                 size_t size)
{
    /* Find string's length */
    size_t len = strlen(src);

    /* Cap length at (num - 1) to leave room for \0 */
    if (size <= len) {
        len = (size - 1);
    }

    /* Copy as much of string as we can fit */
    memcpy(dest, src, len);

    /* explicit null termination */
    dest[len] = '\0';

    return dest;
}

/* Uses global var 'OutFilePrefix' */
FILE* my_fopen(const char* fname,
               const char* flag)
{
    FILE* fp;
    int Len;
    char* new_fname = NULL;

    /* Appends a prefix string for output files */
    if (OutFilePrefix) {
        if (strchr(flag, 'w')) {
            Len = 1;    /* NULL char */
            Len += strlen(OutFilePrefix);
            Len += strlen(fname);
            new_fname = (char*)my_malloc(Len * sizeof(char));
            strcpy(new_fname, OutFilePrefix);
            strcat(new_fname, fname);
            fname = new_fname;
        }
    }

    if (NULL == (fp = fopen(fname, flag))) {
        printf("Error opening file %s for %s access.\n", fname, flag);
        exit(1);
    }

    if (new_fname) {
        free(new_fname);
    }

    return (fp);
}

char* my_strdup(const char* str)
{
    int   Len = 1 + strlen(str);
    char* Dst = (char*)my_malloc(Len * sizeof(char));
    memcpy(Dst, str, Len);

    return Dst;
}

int my_atoi(const char* str)
{
    /* Returns the integer represented by the first part of the character       *
     * string.                                              */

    if (str[0] < '0' || str[0] > '9') {
        if (!(str[0] == '-' && str[1] >= '0' && str[1] <= '9')) {
            printf(ERRTAG "expected number instead of '%s'.\n", str);
            exit(1);
        }
        return (-1);
    }

    return (atoi(str));
}


void* my_calloc(size_t nelem,
                size_t size)
{
    void* ret;
    if ((ret = calloc(nelem, size)) == NULL) {
        fprintf(stderr, "Error:  Unable to calloc memory.  Aborting.\n");
        exit(1);
    }

    return (ret);
}


void* my_malloc(size_t size)
{
    void* ret;
    if ((ret = malloc(size)) == NULL) {
        fprintf(stderr, "Error:  Unable to malloc memory.  Aborting.\n");
        abort();
        exit(1);
    }

    return (ret);
}

void* my_realloc(void* ptr,
                 size_t size)
{
    void* ret = NULL;
    if (size <= 0) {
        printf("reallocating of size <= 0.\n");
    }
    if ((ret = realloc(ptr, size)) == NULL) {
        printf(ERRTAG "Unable to realloc memory. Aborting. "
               "ptr=%p, Size=%d.\n", ptr, size);

        if (ptr == NULL) {
            printf(ERRTAG "my_realloc: ptr == NULL. Aborting.\n");
        }

        exit(1);
    }
    return (ret);
}

/* This routine should be used for allocating fairly small data             *
     * structures where memory-efficiency is crucial.  This routine allocates   *
     * large "chunks" of data, and parcels them out as requested.  Whenever     *
     * it mallocs a new chunk it adds it to the linked list pointed to by       *
     * chunk_ptr_head.  This list can be used to free the chunked memory.       *
     * If chunk_ptr_head is NULL, no list of chunked memory blocks will be kept *
     * -- this is useful for data structures that you never intend to free as   *
     * it means you don't have to keep track of the linked lists.               *
     * Information about the currently open "chunk" must be stored by the       *
     * user program.  mem_avail_ptr points to an int storing how many bytes are *
     * left in the current chunk, while next_mem_loc_ptr is the address of a    *
     * pointer to the next free bytes in the chunk.  To start a new chunk,      *
     * simply set *mem_avail_ptr = 0.  Each independent set of data structures  *
     * should use a new chunk.                                                  */
void* my_chunk_malloc(size_t size,
                      linked_vptr_t** chunk_ptr_head,
                      int* mem_avail_ptr,
                      char** next_mem_loc_ptr)
{
    /* To make sure the memory passed back is properly aligned, I must *
     * only send back chunks in multiples of the worst-case alignment  *
     * restriction of the machine.  On most machines this should be    *
     * a long, but on 64-bit machines it might be a long long or a     *
     * double.  Change the typedef below if this is the case.          */
    typedef long Align;

#define CHUNK_SIZE 32768
#define FRAGMENT_THRESHOLD 100

    char* tmp_ptr = NULL;
    int aligned_size;

    assert(*mem_avail_ptr >= 0);
    if ((size_t)(*mem_avail_ptr) < size) {
        /* Need to malloc more memory. */
        if (size > CHUNK_SIZE) {
            /* Too big, use standard routine. */
            tmp_ptr = (char*)my_malloc(size);

            /* When debugging, uncomment the code below to see if memory allocation size */
            /* makes sense */
            /*#ifdef DEBUG
                   printf("NB:  my_chunk_malloc got a request for %d bytes.\n",
                      size);
                   printf("You should consider using my_malloc for such big requests.\n");
            #endif */

            if (chunk_ptr_head != NULL) {
                *chunk_ptr_head = insert_in_vptr_list(*chunk_ptr_head, tmp_ptr);
            }

            return tmp_ptr;
        }

        if (*mem_avail_ptr < FRAGMENT_THRESHOLD) {
            /* Only a small scrap left. */
            *next_mem_loc_ptr = (char*)my_malloc(CHUNK_SIZE);
            *mem_avail_ptr = CHUNK_SIZE;

            if (chunk_ptr_head != NULL)
                *chunk_ptr_head = insert_in_vptr_list(*chunk_ptr_head,
                                                      *next_mem_loc_ptr);
        } else {
        /* Execute else clause only when the chunk we want is pretty big,  *
         * and would leave too big an unused fragment.  Then we use malloc *
         * to allocate normally.                                           */
            tmp_ptr = (char*)my_malloc(size);
            if (chunk_ptr_head != NULL) {
                *chunk_ptr_head = insert_in_vptr_list(*chunk_ptr_head,
                                                      tmp_ptr);
            }

            return tmp_ptr;
        }
    }

    /* Find the smallest distance to advance the memory pointer and keep *
     * everything aligned.                                               */
    if (size % sizeof(Align) == 0) {
        aligned_size = size;
    } else {
        aligned_size = size + sizeof(Align) - size % sizeof(Align);
    }

    tmp_ptr = *next_mem_loc_ptr;
    *next_mem_loc_ptr += aligned_size;
    *mem_avail_ptr -= aligned_size;
    return (tmp_ptr);
}

void free_chunk_memory(linked_vptr_t* chunk_ptr_head)
{
    /* Frees the memory allocated by a sequence of calls to my_chunk_malloc. */
    linked_vptr_t* curr_ptr, *prev_ptr;
    curr_ptr = chunk_ptr_head;

    while (curr_ptr != NULL) {
        free(curr_ptr->data_vptr);  /* Free memory "chunk". */
        prev_ptr = curr_ptr;
        curr_ptr = curr_ptr->next;
        free(prev_ptr); /* Free memory used to track "chunk". */
    }
}

/* Inserts a new element at the head of a linked list of void pointers. *
 * Returns the new head of the list.                                    */
linked_vptr_t* insert_in_vptr_list(linked_vptr_t* head,
                                   void* vptr_to_add)
{
    linked_vptr_t* linked_vptr =
        (linked_vptr_t*)my_malloc(sizeof(struct s_linked_vptr));

    linked_vptr->data_vptr = vptr_to_add;
    linked_vptr->next = head;
    return  linked_vptr; /* New head of the list */
}

/* Inserts a new element at the head of a linked list of integers.  Returns  *
 * the new head of the list.  One argument is the address of the head of     *
 * a list of free ilist elements.  If there are any elements on this free    *
 * list, the new element is taken from it.  Otherwise a new one is malloced. */
t_linked_int* insert_in_int_list(t_linked_int* head,
                                 int data,
                                 t_linked_int** free_list_head_ptr)
{
    t_linked_int* linked_int = NULL;
    if (*free_list_head_ptr != NULL) {
        linked_int = *free_list_head_ptr;
        *free_list_head_ptr = linked_int->next;
    } else {
        linked_int = (t_linked_int*) my_malloc(sizeof(t_linked_int));
    }

    linked_int->data = data;
    linked_int->next = head;
    return linked_int;
}

/* This routine truly frees (calls free) all the integer list elements    *
 * on the linked list pointed to by *head, and sets head = NULL.          */
void free_int_list(t_linked_int** int_list_head_ptr)
{
    t_linked_int* linked_int, *next_linked_int;
    linked_int = *int_list_head_ptr;

    while (linked_int != NULL) {
        next_linked_int = linked_int->next;
        free(linked_int);
        linked_int = next_linked_int;
    }

    *int_list_head_ptr = NULL;
}

/* Allocates an integer vector with num_items elements and copies the       *
 * integers from the list pointed to by list_head (of which there must be   *
 * num_items) over to it.  The int_list is then put on the free list, and   *
 * the list_head_ptr is set to NULL.                                        */
void alloc_ivector_and_copy_int_list(t_linked_int** list_head_ptr,
                                     int num_items,
                                     vector_t* ivec,
                                     t_linked_int** free_list_head_ptr)
{
    int i, *list;
    t_linked_int* list_head = *list_head_ptr;
    t_linked_int* linked_int;
    if (num_items == 0) {
        /* Empty list. */
        ivec->nelem = 0;
        ivec->list = NULL;
        if (list_head != NULL) {
            printf(ERRTAG
                   "alloc_ivector_and_copy_int_list: Copied %d elements, "
                   "but list at %p contains more.\n",
                   num_items, (void*)list_head);
            exit(1);
        }

        return;
    }

    ivec->nelem = num_items;
    list = (int*)my_malloc(num_items * sizeof(int));
    ivec->list = list;
    linked_int = list_head;

    for (i = 0; i < num_items - 1; i++) {
        list[i] = linked_int->data;
        linked_int = linked_int->next;
    }

    list[num_items - 1] = linked_int->data;
    if (linked_int->next != NULL) {
        printf("Error in alloc_ivector_and_copy_int_list:\n Copied %d elements, "
               "but list at %p contains more.\n",
               num_items, (void*)list_head);
        exit(1);
    }

    linked_int->next = *free_list_head_ptr;
    *free_list_head_ptr = list_head;
    *list_head_ptr = NULL;
}

static int cont;        /* line continued? */

char* my_fgets(char* buf,
               int max_size,
               FILE* fp)
{
    /* Get an input line, update the line number and cut off *
     * any comment part.  A \ at the end of a line with no   *
     * comment part (#) means continue.                      */
    cont = 0;
    linenum++;

    char* val = fgets(buf, max_size, fp);
    if (val == NULL) {
        return (val);
    }

    /* Check that line completely fit into buffer.  (Flags long line   *
     * truncation).                                                    */
    int i;
    for (i = 0; i < max_size; i++) {
        if (buf[i] == '\n') {
            break;
        }

        if (buf[i] == '\0') {
            printf("Error on line %d -- line is too long for input buffer.\n",
                   linenum);
            printf("All lines must be at most %d characters long.\n",
                   BUFSIZE - 2);
            printf("The problem could also be caused by a missing newline.\n");
            exit(1);
        }
    }


    for (i = 0; i < max_size && buf[i] != '\0'; i++) {
        if (buf[i] == '#') {
            buf[i] = '\0';
            break;
        }
    }

    if (i < 2) {
        return (val);
    }

    if (buf[i - 1] == '\n' && buf[i - 2] == '\\') {
        cont = 1;   /* line continued */
        buf[i - 2] = '\n'; /* May need this for tokens */
        buf[i - 1] = '\0';
    }

    return (val);
}

/* Get next token, and wrap to next line if \ at end of line.    *
 * There is a bit of a "gotcha" in strtok.  It does not make a   *
 * copy of the character array which you pass by pointer on the  *
 * first call.  Thus, you must make sure this array exists for   *
 * as long as you are using strtok to parse that line.  Don't    *
 * use local buffers in a bunch of subroutines calling each      *
 * other; the local buffer may be overwritten when the stack is  *
 * restored after return from the subroutine.                    */
char* my_strtok(char* ptr,
                char* tokens,
                FILE* fp,
                char* buf)
{
    char* val = strtok(ptr, tokens);
    for (;;) {
        if (val != NULL || cont == 0) {
            return (val);
        }

        /* return unless we have a null value and a continuation line */
        if (my_fgets(buf, BUFSIZE, fp) == NULL) {
            return (NULL);
        }

        val = strtok(buf, tokens);
    }
}


void free_ivec_vector(vector_t* ivec_vector,
                      int min_row,
                      int max_row)
{
    /* Frees a 1D array of integer vectors.                              */
    int i;
    for (i = min_row; i <= max_row; i++)
        if (ivec_vector[i].nelem != 0) {
            free(ivec_vector[i].list);
        }

    free(ivec_vector + min_row);
}


void free_ivec_matrix(vector_t** ivec_matrix,
                      int min_row,
                      int max_row,
                      int min_col,
                      int max_col)
{
    /* Frees a 2D matrix of integer vectors (ivecs).                     */
    int i, j;

    for (i = min_row; i <= max_row; i++) {
        for (j = min_col; j <= max_col; j++) {
            if (ivec_matrix[i][j].nelem != 0) {
                free(ivec_matrix[i][j].list);
            }
        }
    }

    free_matrix(ivec_matrix, min_row, max_row, min_col, sizeof(vector_t));
}


void free_ivec_matrix3(vector_t** *ivec_matrix3,
                       int min_row,
                       int max_row,
                       int min_col,
                       int max_col,
                       int ndmin,
                       int ndmax)
{
    /* Frees a 3D matrix of integer vectors (ivecs).                     */
    int i, j, k;
    for (i = min_row; i <= max_row; i++) {
        for (j = min_col; j <= max_col; j++) {
            for (k = ndmin; k <= ndmax; k++) {
                if (ivec_matrix3[i][j][k].nelem != 0) {
                    free(ivec_matrix3[i][j][k].list);
                }
            }
        }
    }

    free_matrix3(ivec_matrix3, min_row, max_row, min_col, max_col, ndmin,
                 sizeof(vector_t));
}

/* allocates an generic matrix with (max_row-min_row + 1) rows and max_col - *
 * min_col + 1 columns, with each element of size elsize. i.e.           *
 * returns a pointer to a storage block [min_row..max_row][min_col..max_col].  *
 * Simply cast the returned array pointer to the proper type.          */
void** alloc_matrix(int min_row,
                    int max_row,
                    int min_col,
                    int max_col,
                    size_t elem_size)
{
    char** cptr = (char**)my_malloc((max_row - min_row + 1) * sizeof(char*));
    cptr -= min_row;

    /* Row-based Matrix */
    int i = 0;
    for (i = min_row; i <= max_row; ++i) {
        cptr[i] = (char*)my_malloc((max_col - min_col + 1) * elem_size);
        cptr[i] -= min_col * elem_size / sizeof(char);
    }

    return ((void**)cptr);
} /* end of alloc_matrix */


/* NB:  need to make the pointer type void * instead of void ** to allow   *
 * any pointer to be passed in without a cast.                             */

void free_matrix(void* vptr,
                 int min_row,
                 int max_row,
                 int min_col,
                 size_t elsize)
{
    char** cptr = (char**)vptr;

    int i;
    for (i = min_row; i <= max_row; i++) {
        free(cptr[i] + min_col * elsize / sizeof(char));
    }

    free(cptr + min_row);
}

/* allocates a 3D generic matrix with max_row-min_row + 1 rows, max_col - *
 * min_col + 1 columns, and a depth of ndmax-ndmin + 1, with each         *
 * element of size elsize. i.e. returns a pointer to a storage block      *
 * [min_row..max_row][min_col..max_col][ndmin..ndmax].  Simply cast the   *
 *  returned array pointer to the proper type.                       */
void*** alloc_matrix3(int min_row,
                      int max_row,
                      int min_col,
                      int max_col,
                      int ndmin,
                      int ndmax,
                      size_t elsize)
{
    int i, j;
    char*** cptr = (char***)my_malloc((max_row - min_row + 1) * sizeof(char**));
    cptr -= min_row;

    for (i = min_row; i <= max_row; i++) {
        cptr[i] = (char**)my_malloc((max_col - min_col + 1) * sizeof(char*));
        cptr[i] -= min_col;

        for (j = min_col; j <= max_col; j++) {
            cptr[i][j] = (char*)my_malloc((ndmax - ndmin + 1) * elsize);
            cptr[i][j] -= ndmin * elsize / sizeof(char);
        }
    }

    return ((void***)cptr);
}

void print_int_matrix3(int** *vptr,
                       int min_row,
                       int max_row,
                       int min_col,
                       int max_col,
                       int ndmin,
                       int ndmax,
                       char* file)
{
    FILE* outfile = my_fopen(file, "w");

    int i, j, k;
    for (k = min_row; k <= max_row; ++k) {
        fprintf(outfile, "Plane %d\n", k);

        for (j = min_col; j <= max_col; ++j) {
            for (i = ndmin; i <= ndmax; ++i) {
                fprintf(outfile, "%d ", vptr[k][j][i]);
            }

            fprintf(outfile, "\n");
        }

        fprintf(outfile, "\n");
    }

    fclose(outfile);
}

void free_matrix3(void* vptr,
                  int min_row,
                  int max_row,
                  int min_col,
                  int max_col,
                  int ndmin,
                  size_t elsize)
{
    int i, j;
    char*** cptr = (char***)vptr;

    for (i = min_row; i <= max_row; i++) {
        for (j = min_col; j <= max_col; j++) {
            free(cptr[i][j] + ndmin * elsize / sizeof(char));
        }

        free(cptr[i] + min_col);
    }

    free(cptr + min_row);
}

/* allocates a 3D generic matrix with max_row-min_row + 1 rows, max_col -  *
 * min_col + 1 columns, and a depth of ndmax-ndmin + 1, with each      *
 * element of size elsize. i.e. returns a pointer to a storage block *
 * [min_row..max_row][min_col..max_col][ndmin..ndmax].  Simply cast the      *
 *  returned array pointer to the proper type.                       */
void**** alloc_matrix4(int min_row,
                       int max_row,
                       int min_col,
                       int max_col,
                       int ndmin,
                       int ndmax,
                       int nemin,
                       int nemax,
                       size_t elsize)
{
    int i, j, k;
    char**** cptr = (char****)my_malloc((max_row - min_row + 1) * sizeof(char***));
    cptr -= min_row;

    for (i = min_row; i <= max_row; i++) {
        cptr[i] = (char***)my_malloc((max_col - min_col + 1) * sizeof(char**));
        cptr[i] -= min_col;

        for (j = min_col; j <= max_col; j++) {
            cptr[i][j] = (char**)my_malloc((ndmax - ndmin + 1) *
                                              sizeof(char*));
            cptr[i][j] -= ndmin;

            for (k = ndmin; k <= ndmax; k++) {
                cptr[i][j][k] = (char*)my_malloc((nemax - nemin + 1) *
                                                     elsize);
                cptr[i][j][k] -= nemin * elsize / sizeof(char); /* sizeof(char) = 1) */
            }
        }
    }

    return ((void****)cptr);
}

void free_matrix4(void* vptr,
                  int min_row,
                  int max_row,
                  int min_col,
                  int max_col,
                  int ndmin,
                  int ndmax,
                  int nemin,
                  size_t elsize)
{
    char**** cptr = (char****)vptr;

    int i, j, k;
    for (i = min_row; i <= max_row; i++) {
        for (j = min_col; j <= max_col; j++) {
            for (k = ndmin; k <= ndmax; k++) {
                free(cptr[i][j][k] + nemin * elsize / sizeof(char));
            }

            free(cptr[i][j] + ndmin * elsize / sizeof(char));
        }

        free(cptr[i] + min_col);
    }

    free(cptr + min_row);
}

/* Portable random number generator defined below.  Taken from ANSI C by  *
 * K & R.  Not a great generator, but fast, and good enough for my needs. */
#define IA 1103515245u
#define IC 12345u
#define IM 2147483648u
#define CHECK_RAND

static unsigned int current_random = 0;
static unsigned int current_rand[NUM_OF_THREADS];

void my_srandom(int seed)
{
    current_random = (unsigned int)seed;

    int a;
    for (a = 0; a < NUM_OF_THREADS; a++) {
        current_rand[a] = (unsigned int)seed;
    }
}


/* Creates a random integer between 0 and imax, inclusive. i.e.[0..imax] */
int my_irand(int imax)
{
    current_random = current_random * IA + IC;  /* Use overflow to wrap */
    int ival = current_random & (IM - 1);  /* Modulus */
    ival = (int)((double)ival * (double)(imax + 0.999) / (double)IM);

#ifdef CHECK_RAND
    if ((ival < 0) || (ival > imax)) {
        printf("Bad value in my_irand, imax = %d  ival = %d\n",
                imax, ival);
        exit(1);
    }
#endif

    return ival; /* ival >= 0 && ival <= imax */
}

/* Creates a random integer between 0 and imax, inclusive.  i.e. [0..imax] */
int my_irand_parallel(int imax, int myid)
{
    current_rand[myid] = current_rand[myid] * IA + IC;  /* Use overflow to wrap */
    int ival = current_rand[myid] & (IM - 1);   /* Modulus */
    ival = (int)((double)ival * (double)(imax + 0.999) / (double)IM);

#ifdef CHECK_RAND
    if ((ival < 0) || (ival > imax)) {
        printf("Bad value in my_irand, imax = %d  ival = %d\n", imax,
               ival);
        exit(1);
    }
#endif

    return ival; /* ival >= 0 && ival <= imax */
}

/* Creates a random double between 0 and 1.  i.e. [0..1).        */
double my_frand(void)
{
    current_random = current_random * IA + IC;  /* Use overflow to wrap */
    int ival = current_random & (IM - 1);  /* Modulus */
    double fval = (double)ival / (double)IM;

#ifdef CHECK_RAND
    if ((fval < 0) || (fval > 1.)) {
        printf("Bad value in my_frand, fval = %g\n", fval);
        exit(1);
    }
#endif

    return fval;
}

/* Creates a random double between 0 and 1.  i.e. [0..1).        */
double my_frand_parallel(int random)
{
    current_rand[random] = current_rand[random] * IA + IC; /* Use overflow to wrap */
    int ival = current_rand[random] & (IM - 1); /* Modulus */
    double fval = (double)ival / (double)IM;

#ifdef CHECK_RAND
    if ((fval < 0) || (fval > 1.)) {
        printf("Bad value in my_frand, fval = %g\n", fval);
        exit(1);
    }
#endif

    return fval;
}

