

#ifndef FILE_BUFFER_H
#define FILE_BUFFER_H

#include <stdio.h>

#include "utils.h"
#include "string-match.h"

/*
 * file_buffer.h
 * A buffer for reading and writting to files following a very similiar design to the java.nio.Buffer
 */

typedef struct buf_t BUF_T;

/*
 * buf_create
 * takes the initial size of the buffer
 */
BUF_T* buf_create(int capacity);

/*
 * buf_destroy
 * deallocates all memory for the buffer even if some contents have not been read.
 */
void buf_destroy(BUF_T *buf);

/*
 * buf_clear
 * sets the limit to the capacity and the position to zero
 */
void buf_clear(BUF_T *buf);

/*
 * buf_flip
 * sets the limit to the current position and the position to zero
 */
void buf_flip(BUF_T *buf);

/*
 * buf_rewind
 * sets the position to zero
 */
void buf_rewind(BUF_T *buf);

/*
 * buf_compact
 * moves data between the position and the limit to the begining
 * sets the position to after the data and the limit to the capacity
 */
void buf_compact(BUF_T *buf);

/*
 * buf_skip
 * moves the position forward by the specified number of characters
 * dies if skip is negative or the new position would be past the limit
 */
void buf_skip(BUF_T *buf, int skip);

/*
 * buf_mark
 * marks the current position so it can be returned to with reset
 */
void buf_mark(BUF_T *buf);

/*
 * buf_reset
 * resets the position to the mark
 */
void buf_reset(BUF_T *buf);

/*
 * buf_capacity
 * returns the capacity of the buffer
 */
int buf_capacity(BUF_T *buf);

/*
 * buf_set_capacity
 * Changes the capacity to the passed value.
 * If the new capacity is smaller than the
 * limit or position they are updated to 
 * equal it. If the new capacity is smaller
 * than the mark it is removed.
 */
void buf_set_capacity(BUF_T *buf, int capacity, BOOLEAN_T update_limit);

/*
 * buf_position
 * returns the current reading/writing position
 */
int buf_position(BUF_T *buf);

/*
 * buf_set_position
 * sets the current reading/writing position
 * the new positon must non-negative 
 * and smaller than the limit
 */
void buf_set_position(BUF_T *buf, int position);

/*
 * buf_limit
 * returns the current limit
 */
int buf_limit(BUF_T *buf);

/*
 * buf_set_limit
 * sets the limit, the new limit must be non-negative 
 * and smaller than the buffer's capacity
 */
void buf_set_limit(BUF_T *buf, int limit);

/*
 * buf_remaining
 * returns the number of chars remaining between the position and limit
 */
int buf_remaining(BUF_T *buf);

/*
 * buf_find_next
 * returns the position of the next delimiter or -1 if it is not found, 
 * does not change the buffer in any way.
 *
 * Determines where the next delimiter is by calling is_delim with the 
 * passed config on each sequential character and negating the function's
 * output if negate_is_delim is set.
 *
 * see the helper functions buf_is_delim_letter and buf_is_delim_func.
 */
int buf_find_next(BUF_T *buf, int (*is_delim)(void*, int), void* config, BOOLEAN_T negate_is_delim);

/*
 * buf_is_delim_letter
 * helper function for buf_find_next to allow specification of letters as delimiters.
 * users should specify the function as the is_delim parameter and a null-terminated
 * string of delimiters as the config parameter to the buf_find_next function.
 */
int buf_is_delim_letter(void *delims, int letter);

/*
 * buf_is_delim_space
 * helper function for buf_find_next to allow specification of whitespace as delimiters.
 * users should specify the function as the is_delim parameter, the config parameter is
 * unused.
 */
int buf_is_delim_space(void *ignored, int letter);

/*
 * buf_is_delim_digit
 * helper function for buf_find_next to allow specification of numbers as delimiters.
 * users should specify the function as the is_delim parameter, the config parameter is
 * unused.
 */
int buf_is_delim_digit(void *ignored, int letter);

/*
 * buf_is_delim_alpha
 * helper function for buf_find_next to allow specification of a-zA-Z as delimiters.
 * users should specify the function as the is_delim parameter, the config parameter is
 * unused.
 */
int buf_is_delim_alpha(void *ignored, int letter);

/*
 * buf_kmp_table
 * builds a partial match table for use by the kmp search when searching
 * for the match string.
 */
void buf_kmp_table(char *match_string, int *kmp_table, int len); 

/*
 *  buf_kmp_search
 *  uses a kmp search to move the position forward to the first possible 
 *  character of the match string, returns the number of characters 
 *  remaining to match if it gets to the end of the buffer before completely
 *  matching the string.
 */
int buf_kmp_search(BUF_T *buf, char *match_string, int *kmp_table, int len);

/*
 * buf_bm_search
 * Boyer-Moore search
 * moves the buffer forward to the first match or partial
 * match. Returns the index of the matched string or 
 * negative min chars to match.
 * If a partial match starts before a full match, the 
 * algorithm will return the partial.
 */
int buf_bm_search(BUF_T *buf, int count, BMSTR_T **strings);

/*
 * buf_getc
 * returns the next avaliable character as unsigned char 
 * or -1 if no more characters
 */
int buf_getc(BUF_T *buf);

/*
 * buf_getstr
 * gets a maximum of n bytes from the buffer to the string
 * if the buffer finishes before copying the full n bytes it
 * will null terminate the string
 * returns the number of bytes copied, not including the added
 * null byte
 */
int buf_getstr(BUF_T *buf, char *dest, int n);

/*
 * buf_get_token
 * returns a string which ends at (ie. does not include) the first character 
 * that is_delim designates is a delimeter or the end of the buffer 
 * if the done flag is set. If it reaches the end of the buffer and it 
 * doesn't find a delimiter and input isn't done it returns NULL. 
 * The parameter strlen will be set to the length of the found string or -1 
 * if the string is not found (in which case the buffer position is not changed).
 * The returned string is always null terminated.
 *
 * Allocation
 * When target is specified it will use it as a destination for the token,
 * it will be filled with up to size -1 characters from the buffer.
 * When target is NULL it will allocate a string to fit up to size -1
 * characters. If size is 0, any length is allowed to be allocated. 
 *
 * Errors
 * If the token doesn't fit in the maxlen -1 space then NULL is returned 
 * and strlen is set to the unreturned token length (so anything allocating
 * memory would need to add one). 
 *
 */
char* buf_get_token(BUF_T *buf, int (*is_delim)(void*, int), void* config, BOOLEAN_T negate_is_delim, 
    BOOLEAN_T done, char *target, int size, int* strlen);

/*
 * buf_consume
 * Skips over any delimiters that are at the front of the buffer.
 * returns the number of delimiters skipped
 */
int buf_consume(BUF_T *buf, int (*is_delim)(void*, int), void* config, BOOLEAN_T negate_is_delim);

/*
 * buf_printf
 * prints the format to the buffer
 * returns the 0 on success or the number of bytes needed
 * on failure
 */
int buf_printf(BUF_T *buf, const char *fmt, ...);

/*
 * buf_unexpected
 * checks that the specified string follows exactly. If what follows is
 * unexpected then return 1, otherwise advance the buffer to skip over 
 * the expected string and return 0.
 *
 * Errors
 * The expected string must be smaller than the buffer or this will 
 * cause a fatal error.
 */
int buf_unexpected(BUF_T *buf, char *expected, BOOLEAN_T ignorecase);

/*
 * buf_fread
 * reads data from the file ptr into the buffer
 */
int buf_fread(BUF_T *buf, FILE *fp);

/*
 * buf_fwrite
 * writes data from the buffer to the file ptr
 */
int buf_fwrite(BUF_T *buf, FILE *fp);

/*
 * buf_fgets
 * reads data into the buffer but stops if an end of line
 * is encountered.
 * returns the number of bytes read
 */
int buf_fgets(BUF_T *buf, FILE *fp);

/*
 * buf_fread_until
 * reads data from the file ptr into the buffer until
 * one of the search strings is encountered or an error
 * occurs. Will try to read from the buffer before doing
 * any IO so buffer must be ready to be read.
 * returns the number of the string discovered or zero 
 * for eof and -1 for error.
 * If outfp is not null then data will be copied to it.
 */
int buf_fread_until(BUF_T *buf, FILE *fp, FILE *outfp, int count, BMSTR_T **strings);

/*
 * buf_fread_token
 * returns a string which ends at (ie. does not include) the first character 
 * that is_delim designates is a delimeter or the end of the file. 
 * The returned string is always null terminated. 
 * The parameter strlen will be set to the length of the found string.
 *
 * Allocation
 * When target is specified it will use it as a destination for the token,
 * it will be filled with up to size -1 characters from the buffer.
 * When target is NULL it will allocate a string to fit up to size -1
 * characters. If size is 0, any length is allowed to be allocated.
 *
 * Errors
 * If the token doesn't fit in the maxlen -1 space then returns NULL with
 * strlen set to the buffer's -(offset+1) from the starting position.
 * If a file read error occurs then NULL is returned with strlen set to zero.
 */
char* buf_fread_token(BUF_T *buf, FILE *fp, int (*is_delim)(void*, int), 
    void* config, BOOLEAN_T negate_is_delim, char *target, int size, int *strlen);

/*
 * buf_fread_consume
 * consumes all the delimiters returning the number of delimiters that it consumes.
 * returns -1 on read error.
 */
int buf_fread_consume(BUF_T *buf, FILE *fp, int (*is_delim)(void*, int), void* config, BOOLEAN_T negate_is_delim);

/*
 * buf_fread_next_line
 * reads until the start of a new line is in the buffer
 * returns the number of characters skipped.
 */
int buf_fread_next_line(BUF_T *buf, FILE *fp);

/*
 * buf_fread_unexpected
 * checks that the specified string follows exactly. If what follows is
 * unexpected then return 1, otherwise advance the buffer to skip over 
 * the expected string and return 0.
 *
 * Errors
 * The expected string must be smaller than the buffer or this will 
 * cause a fatal error. If a file read error occurs than it returns -1;
 */
int buf_fread_unexpected(BUF_T *buf, FILE *fp, char *expected, BOOLEAN_T ignorecase);

#endif
