
#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "buffer.h"


/*
 * This structure defines a buffer
 */
struct buf_t {
  char* buffer;   // the allocated memory for the buffer
  int capacity;   // the total number of chars in the buffer
  int limit;      // the first unwritable/unreadable position
  int position;   // the position to write-to/read-from next 
  int mark;       // the position to jump back to on reset
};


/*
 * buf_create
 * takes the initial size of the buffer
 */
BUF_T* buf_create(int capacity) {
  BUF_T *buf = mm_malloc(sizeof(BUF_T));
  buf->buffer = mm_malloc(sizeof(char)*capacity);
  buf->capacity = capacity;
  buf->limit = capacity;
  buf->position = 0;
  buf->mark = -1;
  return buf;
}

/*
 * buf_destroy
 * deallocates all memory for the buffer even if some contents have not been read.
 */
void buf_destroy(BUF_T *buf) {
  free(buf->buffer);
  memset(buf, 0, sizeof(BUF_T));
  free(buf);
}

/*
 * buf_clear
 * sets the limit to the capacity and the position to zero
 * the mark is discarded
 */
void buf_clear(BUF_T *buf) {
  buf->limit = buf->capacity;
  buf->position = 0;
  buf->mark = -1;
}

/*
 * buf_flip
 * sets the limit to the current position and the position to zero
 * the mark is discarded
 */
void buf_flip(BUF_T *buf) {
  buf->limit = buf->position;
  buf->position = 0;
  buf->mark = -1;
}

/*
 * buf_rewind
 * sets the position to zero
 */
void buf_rewind(BUF_T *buf) {
  buf->position = 0;
  buf->mark = -1;
}

/*
 * buf_compact
 * moves data between the position and the limit to the begining
 * sets the position to after the data and the limit to the capacity
 */
void buf_compact(BUF_T *buf) {
  int len;
  char *src, *dest;
  len = (buf->limit - buf->position);
  src = (buf->buffer)+(buf->position);
  dest = buf->buffer;
  memmove(dest, src, len*sizeof(char));
  buf->position = len;
  buf->limit = buf->capacity;
  buf->mark = -1;
}

/*
 * buf_skip
 * moves the position forward by the specified number of characters
 * dies if skip is negative or the new position would be past the limit
 */
void buf_skip(BUF_T *buf, int skip) {
  if (skip < 0) {
    die("buf_skip: skip (%d) must be positive\n", skip);
  }
  if (buf->position + skip > buf->limit) {
    die("buf_skip: can't skip (%d + %d) past the limit (%d)\n", 
        buf->position, skip, buf->limit);
  }
  buf->position += skip;
}

/*
 * buf_mark
 * marks the current position so it can be returned to with reset
 */
void buf_mark(BUF_T *buf) {
  buf->mark = buf->position;
}

/*
 * buf_reset
 * resets the position to the mark
 */
void buf_reset(BUF_T *buf) {
  if (buf->mark == -1) 
    die("Mark is not set!\n");
  buf->position = buf->mark;
}

/*
 * buf_capacity
 * returns the capacity of the buffer
 */
int buf_capacity(BUF_T *buf) {
  return buf->capacity;
}

/*
 * buf_set_capacity
 * Changes the capacity to the passed value.
 * If the new capacity is smaller than the
 * limit or position they are updated to 
 * equal it. If the new capacity is smaller
 * than the mark it is removed.
 */
void buf_set_capacity(BUF_T *buf, int capacity, BOOLEAN_T update_limit) {
  buf->buffer = mm_realloc(buf->buffer, sizeof(char)*capacity);
  buf->capacity = capacity;
  if (buf->limit > capacity) { 
    buf->limit = capacity;
    if (buf->position > capacity) {
      buf->position = capacity;
      if (buf->mark > capacity)
        buf->mark = -1;
    }
  } else if (update_limit) {
    buf->limit = capacity;
  }
}

/*
 * buf_position
 * returns the current reading/writing position
 */
int buf_position(BUF_T *buf) {
  return buf->position;
}

/*
 * buf_set_position
 * sets the current reading/writing position
 * the new positon must non-negative 
 * and not larger than the limit
 */
void buf_set_position(BUF_T *buf, int position) {
  if (position < 0 || position > buf->limit) {
    die("Position out of bounds\n");
  }
  if (buf->mark > position) buf->mark = -1;
  buf->position = position;
}

/*
 * buf_limit
 * returns the current limit
 */
int buf_limit(BUF_T *buf) {
  return buf->limit;
}

/*
 * buf_set_limit
 * sets the limit, the new limit must be non-negative 
 * and smaller than the buffer's capacity
 */
void buf_set_limit(BUF_T *buf, int limit) {
  if (limit < 0 || limit > buf->capacity) {
    die("limit out of bounds\n");
  }
  if (buf->position > limit) {
    buf->position = limit;
    if (buf->mark > limit) buf->mark = -1;
  }
  buf->limit = limit;
}

/*
 * buf_remaining
 * returns the number of chars remaining between the position and limit
 */
int buf_remaining(BUF_T *buf) {
  return buf->limit - buf->position;
}

/*
 * buf_find_next
 * returns the position of the next delimiter or -1 if it is not found, 
 * does not change the buffer in any way.
 *
 * Determines where the next delimiter is by calling is_delim with the 
 * passed config on each sequential character and negating the function's
 * output if negate_is_delim is set.
 *
 * see the helper function buf_is_delim_letter
 */
int buf_find_next(BUF_T *buf, int (*is_delim)(void*, int), void* config, BOOLEAN_T negate_is_delim) {
  int i, j;
  //look for a delimiter or the end
  if (negate_is_delim) {
    for (i = buf->position; i < buf->limit; ++i) {
      if (!is_delim(config, (buf->buffer)[i])) {
        return i;
      }
    }
  } else {
    for (i = buf->position; i < buf->limit; ++i) {
      if (is_delim(config, (buf->buffer)[i])) {
        return i;
      }
    }
  }
  return -1;
}

/*
 * buf_is_delim_letter
 * helper function for buf_find_next to allow specification of letters as delimiters.
 * users should specify the function as the is_delim parameter and a null-terminated
 * string of delimiters as the config parameter to the buf_find_next function.
 */
int buf_is_delim_letter(void *delims, int letter) {
  char *letters = (char*)delims;
  while (*letters != '\0') {
    if (*letters == letter) return TRUE;
    ++letters;
  }
  return FALSE;
}

/*
 * buf_is_delim_space
 * helper function for buf_find_next to allow specification of whitespace as delimiters.
 * users should specify the function as the is_delim parameter, the config parameter is
 * unused.
 */
int buf_is_delim_space(void *ignored, int letter) {
  return isspace(letter);
}

/*
 * buf_is_delim_digit
 * helper function for buf_find_next to allow specification of numbers as delimiters.
 * users should specify the function as the is_delim parameter, the config parameter is
 * unused.
 */
int buf_is_delim_digit(void *ignored, int letter) {
  return isdigit(letter);
}

/*
 * buf_is_delim_alpha
 * helper function for buf_find_next to allow specification of a-zA-Z as delimiters.
 * users should specify the function as the is_delim parameter, the config parameter is
 * unused.
 */
int buf_is_delim_alpha(void *ignored, int letter) {
  return isalpha(letter);
}

/*
 * buf_kmp_table
 * builds a partial match table for use by the kmp search when searching
 * for the match string.
 */
void buf_kmp_table(char *match_string, int *kmp_table, int len) {
  int pos, candidate;

  kmp_table[0] = -1;
  kmp_table[1] = 0;
  pos = 2;
  candidate = 0;

  while (pos < len) {
    if (match_string[pos - 1] == match_string[candidate]) {
      kmp_table[pos] = candidate + 1;
      pos += 1;
      candidate += 1;
    }
    else if (candidate > 0) {
      candidate = kmp_table[candidate];
    } else {
      kmp_table[pos] = 0;
      pos += 1;
    }
  }
}


/*
 *  buf_kmp_search
 *  uses a kmp search to move the position forward to the first possible 
 *  character of the match string, returns the number of characters 
 *  remaining to match if it gets to the end of the buffer before completely
 *  matching the string.
 */
int buf_kmp_search(BUF_T *buf, char *match_string, int *kmp_table, int len) {
  int m, i;
  m = buf->position;
  i = 0;
  while (m + i < buf->limit) {
    if (match_string[i] == buf->buffer[m+i]) {
      i += 1;
      if (i == len) {
        buf->position = m;
        return 0;
      }
    } else {
      m = m + i - kmp_table[i];
      if (i > 0) {
        i = kmp_table[i];
      }
    }
  }
  //incomplete match, return the number of characters not matched
  buf->position = m;
  return len - i;
}

/*
 * buf_bm_search
 * Boyer-Moore search
 * moves the buffer forward to the first match or partial
 * match. Returns the index of the matched string or 
 * negative min chars to match.
 * If a partial match starts before a full match, the 
 * algorithm will return the partial.
 */
int buf_bm_search(BUF_T *buf, int count, BMSTR_T **strings) {
  int len, i, pos, pos_stm, left_i, left, left_stm;
  char *txt;
  txt = (buf->buffer)+(buf->position);
  len = (buf->limit - buf->position);
  left_i = -1;
  left = 0;//stop compilier complaining about uninitilized vars
  left_stm = 0;
  for (i = 0; i < count; ++i) {
    pos = bmstr_substring(strings[i], txt, len);
    if (pos < 0) {
      pos = -(pos + 1);
      pos_stm = bmstr_length(strings[i]) - (len - pos);
    } else {
      pos_stm = 0;
    }
    if (left_i == -1 || pos < left) {
      left_i = i;
      left = pos;
      left_stm = pos_stm;
    }
  }
  if (left_i == -1) {
    return -1;//this should only happen if count is zero, which shouldn't happen
  } else {
    buf->position += left;
    if (left_stm) return -left_stm;
    return left_i;
  }
}

/*
 * buf_getc
 * returns the next avaliable character as unsigned char 
 * or -1 if no more characters
 */
int buf_getc(BUF_T *buf) {
  int rvalue;

  if (buf->position < buf->limit) {
    rvalue = (((unsigned char*)buf->buffer)[buf->position]);
    buf->position += 1;
  } else {
    rvalue = -1;
  }
  return rvalue;
}

/*
 * buf_getstr
 * gets a maximum of n bytes from the buffer to the string
 * if the buffer finishes before copying the full n bytes it
 * will null terminate the string
 * returns the number of bytes copied, not including the added
 * null byte
 */
int buf_getstr(BUF_T *buf, char *dest, int n) {
  int i, j;
  for (i = buf->position, j = 0; i < buf->limit && j < n; ++i, ++j) {
    dest[j] = (buf->buffer)[i];
  }
  buf->position += j;
  if (j < n) {
    dest[j] = '\0';
  }
  return j;
}

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
    BOOLEAN_T done, char *target, int size, int* slen) {
  int index;
  int needed_size;
  index = buf_find_next(buf, is_delim, config, negate_is_delim);

  if (index == -1 && done) {
    index = buf->limit;
  }
  
  if (index > buf->position) {
    needed_size = index - buf->position;
    if (slen) *slen = needed_size;
    if (target == NULL) {
      if (size == 0 || size > needed_size) {
        target = (char*)mm_malloc(needed_size+1);
      } else {//doesn't fit in allowed size
        return NULL;
      }
    } else {
      if (size <= needed_size) {
        return NULL;
      }
    }

    strncpy(target, (buf->buffer)+(buf->position), needed_size);
    target[needed_size] = '\0';
    buf->position = index;
  } else { //token not found
    if (slen) *slen = -1;
    target = NULL;
  }
  return target;
}

/*
 * buf_consume
 * Skips over any delimiters that are at the front of the buffer.
 * returns the number of delimiters skipped
 */
int buf_consume(BUF_T *buf, int (*is_delim)(void*, int), void* config, BOOLEAN_T negate_is_delim) {
  int index, count;
  index = buf_find_next(buf, is_delim, config, !negate_is_delim);
  if (index == -1) index = buf->limit;
  count = index - buf->position;
  buf->position = index;
  return count;
}

/*
 * buf_get_token
 * returns an allocated string which ends at any of the specified delimiters
 * or the end of the buffer if the done flag is set. If it reaches the end of
 * the buffer and it doesn't find a delimiter and input isn't done it returns
 * NULL. Automatically skips any of the specified delimeters at the start. 
 */
/*
char* buf_get_token(BUF_T *buf, char *delimiters, int delimiter_count, BOOLEAN_T done) {
  int i, j, start, len;
  BOOLEAN_T isdelim;
  char *token;

  //skip delimiters at the start
  for (i = buf->position; i < buf->limit; ++i) {
    isdelim = FALSE;
    for (j = 0; j < delimiter_count; ++j) {
      if ((buf->buffer)[i] == delimiters[j]) {
        isdelim = TRUE;
        break;
      }
    }
    if (!isdelim) break;
  }
  buf->position = i;
  if (i == buf->limit) { //only found delimiters :(
    if (done) { //allocate and return an empty string
      token = (char*)mm_malloc(sizeof(char));
      token[0] = '\0';
      return token;
    }
    return NULL; //need more input
  }
  //record start position
  start = i;
  //look for another delimiter or the end
  for (i = start+1; i < buf->limit; ++i) {
    isdelim = FALSE;
    for (j = 0; j < delimiter_count; ++j) {
      if ((buf->buffer)[i] == delimiters[j]) {
        isdelim = TRUE;
        break;
      }
    }
    if (isdelim) break;
  }
  if (i < buf->limit || done) { // found delimiter or end of input
    buf->position = i;
    len = i - start;
    token = (char*)mm_malloc(sizeof(char)*(len+1));
    strncpy(token, (buf->buffer)+start, len);
    token[len] = '\0';
    return token;
  }
  return NULL; //need more input
}*/

/*
 * buf_printf
 * prints the format to the buffer
 * returns the 0 on success or the number of bytes needed
 * on failure
 */
int buf_printf(BUF_T *buf, const char *fmt, ...) {
  int remaining, rvalue;
  char *dest;
  va_list args;
  remaining = buf->limit - buf->position;
  dest = (buf->buffer)+(buf->position);
  va_start(args, fmt);
  rvalue = vsnprintf(dest, remaining*sizeof(char), fmt, args);
  va_end(args);
  if (rvalue <= remaining) {
    //note that \0 byte is not included in rvalue
    //hence it will be overwritten in the next call
    buf->position += rvalue; 
    return 0;
  } else {
    //return the number of bytes needed
    return rvalue - remaining;
  }
}

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
int buf_unexpected(BUF_T *buf, char *expected, BOOLEAN_T ignorecase) {
  int i;
  i = buf->position;
  while (*expected != '\0') {
    if (i == buf->limit) {
      if (buf->position == 0 && buf->limit == buf->capacity)
        die("buf_unexpected() - expected string is longer than the buffer\n");
      //can't match all characters so it is unexpected
      return 1;
    }
    if (ignorecase) {
      if (tolower(buf->buffer[i]) != tolower(*expected)) return 1;
    } else {
      if (buf->buffer[i] != *expected) return 1;
    }
    ++i;
    ++expected;
  }
  buf->position = i;
  return 0;
}

/*
 * buf_fread
 * reads data from the file ptr into the buffer
 * returns the number of bytes read
 */
int buf_fread(BUF_T *buf, FILE *fp) {
  int remaining, rvalue;
  char *dest;

  remaining = (buf->limit - buf->position);
  if (remaining == 0) return 0;
  dest = (buf->buffer)+(buf->position);
  rvalue = fread(dest, sizeof(char), remaining, fp);
  if (rvalue > 0)
    buf->position += rvalue;
  return rvalue;
}

/*
 * buf_fwrite
 * writes data from the buffer to the file ptr
 * returns the number of bytes written
 */
int buf_fwrite(BUF_T *buf, FILE *fp) {
  int remaining, rvalue;
  char *src;

  remaining = (buf->limit - buf->position);
  if (remaining == 0) return 0;
  src = (buf->buffer)+(buf->position);
  rvalue = fwrite(src, sizeof(char), remaining, fp);
  if (rvalue > 0)
    buf->position += rvalue;
  return rvalue;
}

/*
 * buf_fgets
 * reads data into the buffer but stops if an end of line
 * is encountered.
 * returns the number of bytes read
 */
int buf_fgets(BUF_T *buf, FILE *fp) {
  int remaining, i, rvalue;
  char *dest, *s;

  remaining = (buf->limit - buf->position);
  if (remaining == 0) return 0;
  dest = (buf->buffer)+(buf->position);
  s = fgets(dest, remaining, fp);
  if (s == NULL) {
    return 0;
  }
  for (i = buf->position; i < buf->limit; ++i) {
    if (buf->buffer[i] == '\0') break;
  }
  rvalue = i - buf->position;
  buf->position = i;
  return rvalue;
}

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
int buf_fread_until(BUF_T *buf, FILE *fp, FILE *outfp, int count, BMSTR_T **strings) {
  int index, limit;
  assert(buf != NULL);
  assert(fp != NULL);
  assert(count > 0);
  assert(strings != NULL);
  buf_mark(buf);
  index = buf_bm_search(buf, count, strings);
  while (TRUE) {
    if (outfp) {
      limit = buf_limit(buf);
      buf_set_limit(buf, buf_position(buf));
      buf_reset(buf);
      while (buf_remaining(buf)) {
        buf_fwrite(buf, outfp);
        if (ferror(outfp)) return -1;
      }
      buf_set_limit(buf, limit);
    }
    if (index >= 0) break;
    if (feof(fp)) return 0;
    buf_compact(buf);
    buf_fread(buf, fp);
    buf_flip(buf);
    if (ferror(fp)) return -1;
    buf_mark(buf);
    index = buf_bm_search(buf, count, strings);
  }
  return index + 1;
}

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
    void* config, BOOLEAN_T negate_is_delim, char *target, int size, int *strlen) {
  int index, block, need, offset, size_ok;
  char *token;

  assert(buf != NULL);
  assert(fp != NULL);
  assert(is_delim != NULL);
  assert(size >= 0);

  need = 0;
  offset = 0;
  
  if (target) {
    token = target;
  } else {
    token = NULL;
  } 

  while(TRUE) {
    index = buf_find_next(buf, is_delim, config, negate_is_delim);
    if (index == -1) {
      block = buf_remaining(buf);
    } else {
      block = index - buf->position;
    }
    need += block;

    size_ok = FALSE;
    if (target == NULL && (size == 0 || size > need)) { //we're allowed to allocate more
      size_ok = TRUE;
      if (token) { //the allocation can be increased
        token = mm_realloc(token, sizeof(char)*(need+1));
      } else { //the initial allocation needs to be made
        token = mm_malloc(sizeof(char)*(need+1));
      }
    } else if (target && size > need) {
      size_ok = TRUE;
    }
    if (!size_ok) { 
      //not enough space to read the token
      buf_set_position(buf, buf_limit(buf));//strlen is relative to buffer position
      if (strlen) *strlen = -(need+1); //the user can seek backwards to recover
      if (!target && token) free(token); //cleanup
      return NULL;
    }

    //copy over the data to the token
    offset += buf_getstr(buf, token+offset, block);
    assert(offset == need);
    //check if we found the end of the token
    if (index != -1 || feof(fp)) break;
    //we didn't find the end so need to keep reading
    buf_compact(buf);//empty the buffer
    //load more data into the buffer
    buf_fread(buf, fp);
    buf_flip(buf); //prepare for reading from the buffer
    if (ferror(fp)) { 
      //there's often no real way to recover from this, but the caller might want their own error message
      if (!target && token) free(token);
      if (strlen) *strlen = 0;
      return NULL;
    }
  }//end while

  token[need] = '\0';
  if (strlen) *strlen = need;
  return token;
}

/*
 * buf_fread_consume
 * consumes all the delimiters returning the number of delimiters that it consumes.
 * returns -1 on read error.
 */
int buf_fread_consume(BUF_T *buf, FILE *fp, int (*is_delim)(void*, int), void* config, BOOLEAN_T negate_is_delim) {
  int index, count;
  count = 0;
  while (TRUE) {
    index = buf_find_next(buf, is_delim, config, !negate_is_delim);
    if (index == -1) {
      count += buf_remaining(buf);
    } else {
      count += index - buf->position;
      break;
    }
    if (feof(fp)) {
      index = buf->limit;
      break;
    }
    buf_clear(buf);
    buf_fread(buf, fp);
    buf_flip(buf);
    if (ferror(fp)) {
      return -1;
    }
  }
  buf_set_position(buf, index);
  return count;
}

/*
 * is_newline
 * helper function not visible outside buffer.c
 * specifies which characters are part of a newline.
 * linux:   \n
 * windows: \r\n
 * mac:     \n\r
 * note it's not really a newline if we only see \r
 * but it can cause problems trying to be exact.
 */
static int is_newline(void *config, int letter) {
  return (letter == '\n' || letter == '\r');
}

/*
 * buf_fread_next_line
 * reads until the start of a new line is in the buffer
 * returns the number of characters skipped.
 */
int buf_fread_next_line(BUF_T *buf, FILE *fp) {
  int len, total;
  //consume the rest of the line
  len = buf_fread_consume(buf, fp, is_newline, NULL, TRUE);
  if (len == -1) return -1;
  total = len;
  len = buf_fread_consume(buf, fp, is_newline, NULL, FALSE);
  if (len == -1) return -1;
  total += len;
  return total;
}

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
int buf_fread_unexpected(BUF_T *buf, FILE *fp, char *expected, BOOLEAN_T ignorecase) {
  int i;
  i = buf->position;
  while (*expected != '\0') {
    while (i == buf->limit) {
      i -= buf->position;
      buf_compact(buf);
      if (buf_remaining(buf) == 0)
        die("buf_fread_unexpected() - expected string is longer than the buffer\n");
      if (feof(fp)) { //no more to read and still not matched
        buf_flip(buf);
        return 1;
      }
      buf_fread(buf, fp);
      buf_flip(buf);
      if (ferror(fp)) return -1;
    }
    if (ignorecase) {
      if (tolower(buf->buffer[i]) != tolower(*expected)) return 1;
    } else {
      if (buf->buffer[i] != *expected) return 1;
    }
    ++i;
    ++expected;
  }
  buf->position = i;
  return 0;
}
