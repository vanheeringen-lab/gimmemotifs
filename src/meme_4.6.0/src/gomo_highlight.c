/********************************************************************
 * FILE: gomo_highlight.c
 * AUTHOR: James Johnson 
 * CREATE DATE: 29/09/2009
 * PROJECT: MEME suite
 * COPYRIGHT: 2009, UQ
 *
 * GOMO highlight uses the GO Heirarchy to mark results
 * as being implied when a more specific term in the heirarcy is
 * also a significant result. It also adds additional information
 * relating to the term, such as its name grouping and relative 
 * position in the go hierarchy.
 *
 ********************************************************************/

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "binary-search.h"
#include "buffer.h"
#include "gomo_highlight.h"
#include "utils.h"

/*
 * Stores the structure of the
 * Gene Ontology
 */
typedef struct GONode GONODE_T;
struct GONode {
  char *group;
  char *name;
  char *id;
  int nabove;
  int nbelow;
  GONODE_T **parents;
  int pcount;
  GONODE_T **children;
  int ccount;
  BOOLEAN_T marked;
};

/*
 * Stores the link between a GO identifier
 * and its place in the structure of the
 * Gene Ontology
 */
typedef struct GOTerm GOTERM_T;
struct GOTerm {
  char *id;
  GONODE_T *node;
};

/*
 * Stores the structure and IDs of the GO.
 */
struct GODAG {
  GONODE_T root;
  int node_count;
  GONODE_T *nodes;
  int term_count;
  GOTERM_T *terms;
};

/*
 * Stores a pattern for string matching
 */
typedef struct match MATCH_T;
struct match {
  int length;
  char *string;
  int table[];
};

/*
 * String match patterns used by rewrite_xml
 */
/*
MATCH_T m_group = {
  .length = 8,
  .string = "group=\"\"",
  .table = {-1,0,0,0,0,0,0,0}
};
MATCH_T m_nabove = {
  .length = 9,
  .string = "nabove=\"\"",
  .table = {-1,0,0,0,0,0,0,0,0}
};
MATCH_T m_nbelow = {
  .length = 9,
  .string = "nbelow=\"\"",
  .table = {-1,0,0,0,0,0,0,0,0}
};
MATCH_T m_implied = {
  .length = 11,
  .string = "implied=\"u\"",
  .table = {-1,0,0,0,0,1,0,0,0,0,0}
};
MATCH_T m_name = {
  .length = 7,
  .string = "name=\"\"",
  .table = {-1,0,0,0,0,0,0}
};
MATCH_T m_motif_start = {
  .length = 6,
  .string = "<motif",
  .table = {-1,0,0,0,0,0}
};
MATCH_T m_id_attr = {
  .length = 4,
  .string = "id=\"",
  .table = {-1,0,0,0}
};
MATCH_T m_endtag = {
  .length = 1,
  .string = ">",
  .table = {-1, 0}
};
MATCH_T m_endtag2 = {
  .length = 2,
  .string = "/>",
  .table = {-1, 0}
};
MATCH_T m_motif_end = {
  .length = 7,
  .string = "</motif",
  .table = {-1, 0, 0, 0, 0, 0, 0}
};
MATCH_T m_goterm = {
  .length = 7,
  .string = "<goterm",
  .table = {-1, 0, 0, 0, 0, 0, 0}
};
*/
/*
 * get another token or die trying
 */
/*
char* next_token(BUF_T *buffer, FILE *fp, char *delimiters, int dcount) {
  char *token;
  //try to get a token without doing any reading
  if ((token = buf_get_token(buffer, delimiters, dcount, feof(fp)))) {
    return token;
  }
  //ok so we need to read, prepare the buffer for reading
  buf_compact(buffer);
  while (!feof(fp)) {
    buf_fread(buffer, fp);
    if (ferror(fp)) {
      die("Error occured while reading the file, error was given as: %s\n", strerror(ferror(fp)));
    }
    //have another go at reading the token
    buf_flip(buffer);
    if ((token = buf_get_token(buffer, delimiters, dcount, feof(fp)))) {
      return token;
    }
    //maybe we can read more into the buffer?
    buf_compact(buffer);
    if (!buf_remaining(buffer)) {
      die("Buffer full, but can't read line\n");
    }
  }
  buf_flip(buffer);
  //got to the end of the file
  return buf_get_token(buffer, delimiters, dcount, TRUE);
}
*/
char* next_token(BUF_T *buffer, FILE *fp, int (*is_delim)(void*, int)) {
  char *token;
  int strlen;
  token = buf_fread_token(buffer, fp, is_delim, NULL, FALSE, NULL, 0, &strlen);
  if (token == NULL) {
    if (ferror(fp)) {
      die("Error occured while reading the file, error was given as: %s\n", strerror(ferror(fp)));
    }
    die("Error while reading next token\n");
  }
  return token;
}

void consume(BUF_T *buffer, FILE *fp, int (*is_delim)(void*, int)) {
  int num;
  num = buf_fread_consume(buffer, fp, is_delim, NULL, FALSE);
  if (num < 0) {
    if (ferror(fp)) {
      die("Error occured while reading the file, error was given as: %s\n", strerror(ferror(fp)));
    }
    die("Error while consuming delimiter\n");
  }
}

/*
 * read a specified number of bytes
 */
char* next_chunk(BUF_T *buffer, FILE *fp, int len) {
  int i, c;
  char *chunk;

  if (len <= 0) return NULL;
  chunk = mm_malloc(sizeof(char)*(len+1));
  for (i = 0; i < len; ++i) {
    while ((c = buf_getc(buffer)) == -1) {
      if (feof(fp)) {
        die("File ends with %d bytes still expected\n", (len - i));
      }
      buf_compact(buffer);
      buf_fread(buffer, fp);
      if (ferror(fp)) {
          die("Error occured while reading the file, error was given as: %s\n", strerror(ferror(fp)));
      }
      buf_flip(buffer);
    }
    chunk[i] = (char)c; //narrowing conversion
  }
  chunk[i] = '\0';
  return chunk;
}

/*
 * parse an integer
 */
int parse_pos_int(char *input) {
  long value;
  char *end;
  value = strtol(input, &end, 10);
  assert(*end == '\0');
  if (*end != '\0') {
    die("Not an int, got: \"%s\"\n", input);
  }
  if (value > INT_MAX) {
    die("Too big for int\n");
  }
  if (value < 0) {
    die("Expected positive int\n");
  }
  return (int)value;
}


inline static int is_comment_start(void *config, int letter) {
  return (letter == '#');
}
inline static int is_newline(void *config, int letter) {
  return (letter == '\n' || letter == '\r');
}
inline static int is_tab(void *config, int letter) {
  return (letter == '\t');
}
inline static int is_tab_or_nl(void *config, int letter) {
  return (letter == '\t' || letter == '\n' || letter == '\r');
}
inline static int is_quote(void *config, int letter) {
  return (letter == '"');
}

/*
 * skips over comments indicated by a # at the start of the line
 * expects the buffer to be at the start of a line
 */
void skip_comments(BUF_T *buffer, FILE *fp) {
  //while there's a # character to consume
  while (buf_fread_consume(buffer, fp, is_comment_start, NULL, FALSE) > 0) {
    //consume anything that isn't a newline character
    if (buf_fread_consume(buffer, fp, is_newline, NULL, TRUE) == -1) break;
    //consume the newline
    if (buf_fread_consume(buffer, fp, is_newline, NULL, FALSE) == -1) break;
  }
  if (ferror(fp)) {
    die("Error reading file");
  }
}



/*
 * loads a godag from file
 */
GODAG_T* load_go_dag(char *file) {
  GODAG_T *godag;
  char *token;
  int i, j, token_count, node_index;
  BUF_T *buffer;
  FILE *fp;
  BOOLEAN_T primary;

  if (!file_exists(file)) {
    die("Couldn't load GO DAG as passed file \"%s\" does not exist\n", file);
  }

  godag = (GODAG_T*)mm_malloc(sizeof(GODAG_T));
  memset(&(godag->root), 0, sizeof(GONODE_T));
  godag->root.group = "obsolete";
  godag->root.name = "";
  fp = fopen(file, "r");
  if (fp == NULL) {
    die("Error opening file \"%s\", error was given as: %s\n", file, strerror(errno));
  }
  buffer = buf_create(100);
  buf_flip(buffer);// prepare the buffer for reading because that's what next_token expects
  //skip comments, which are only allowed at the top of the file
  skip_comments(buffer, fp);
  //read the node count
  consume(buffer, fp, is_newline);
  token = next_token(buffer, fp, is_newline);
  consume(buffer, fp, is_newline);
  godag->node_count = parse_pos_int(token);
  free(token);
  godag->nodes = (GONODE_T*)mm_malloc(sizeof(GONODE_T)*godag->node_count);
  for (i = 0; i < godag->node_count; ++i) {
    godag->nodes[i].marked = FALSE;
    godag->nodes[i].id = NULL;
    //read the grouping of the node
    godag->nodes[i].group = next_token(buffer, fp, is_newline);
    consume(buffer, fp, is_newline);
    //parse the name of the node
    token = next_token(buffer, fp, is_tab);
    token_count = parse_pos_int(token);//read the size of the name
    free(token);
    //skip the tab
    buf_skip(buffer, 1);
    //read the name
    godag->nodes[i].name = next_chunk(buffer, fp, token_count);
    //parse the position of the node in the DAG
    consume(buffer, fp, is_newline);
    token = next_token(buffer, fp, is_tab);
    consume(buffer, fp, is_tab);
    godag->nodes[i].nabove = parse_pos_int(token);//read the number of nodes above 
    free(token);
    token = next_token(buffer, fp, is_newline);
    consume(buffer, fp, is_newline);
    godag->nodes[i].nbelow = parse_pos_int(token);//read the number of nodes below
    free(token);
    //read the token count for the parents
    token = next_token(buffer, fp, is_tab_or_nl);
    token_count = parse_pos_int(token);
    free(token);
    if (token_count != 0) {
      consume(buffer, fp, is_tab);
      godag->nodes[i].pcount = token_count;
      godag->nodes[i].parents = mm_malloc(sizeof(GONODE_T*)*token_count);
    } else { // root node, join it with our dummy root node
      consume(buffer, fp, is_newline);
      if (godag->root.ccount == 0) {
        godag->root.ccount = 1;
        godag->root.children = (GONODE_T**)mm_malloc(sizeof(GONODE_T*));
      } else {
        godag->root.ccount += 1;
        godag->root.children = (GONODE_T**)mm_realloc(godag->root.children, sizeof(GONODE_T*)*godag->root.ccount);
      }
      godag->root.children[godag->root.ccount -1] = godag->nodes+i;
      godag->nodes[i].pcount = 1;
      godag->nodes[i].parents = mm_malloc(sizeof(GONODE_T*));
      godag->nodes[i].parents[0] = &(godag->root);
    }
    //load all except the last token (as then we don't have to scan for \r or \n)
    for (j = 0; j < token_count -1; ++j) {
      token = next_token(buffer, fp, is_tab);
      consume(buffer, fp, is_tab);
      node_index = parse_pos_int(token);
      free(token);
      godag->nodes[i].parents[j] = godag->nodes+node_index;
    }
    if (token_count > 0) {//load the last token
      token = next_token(buffer, fp, is_newline);
      consume(buffer, fp, is_newline);
      node_index = parse_pos_int(token);
      free(token);
      godag->nodes[i].parents[j] = godag->nodes+node_index;
    }
    //read the token count for the children
    token = next_token(buffer, fp, is_tab_or_nl);
    consume(buffer, fp, is_tab_or_nl);
    token_count = parse_pos_int(token);
    free(token);
    godag->nodes[i].ccount = token_count;
    godag->nodes[i].children = mm_malloc(sizeof(GONODE_T*)*token_count);
    //load all except the last token (as then we don't have to scan for \r or \n)
    for (j = 0; j < token_count -1; ++j) {
      token = next_token(buffer, fp, is_tab);
      consume(buffer, fp, is_tab);
      node_index = parse_pos_int(token);
      free(token);
      godag->nodes[i].children[j] = godag->nodes+node_index;
    }
    if (token_count > 0) {//load the last token
      token = next_token(buffer, fp, is_newline);
      consume(buffer, fp, is_newline);
      node_index = parse_pos_int(token);
      free(token);
      godag->nodes[i].children[j] = godag->nodes+node_index;
    }
  }//end for each node
  //read in the go term id count
  token = next_token(buffer, fp, is_newline);
  consume(buffer, fp, is_newline);
  godag->term_count = parse_pos_int(token);
  free(token);
  godag->terms = (GOTERM_T*)mm_malloc(sizeof(GOTERM_T)*godag->term_count);
  for (i = 0; i < godag->term_count; ++i) {
    token = next_token(buffer, fp, is_tab);
    consume(buffer, fp, is_tab);
    //detect if it is a primary or alternate id
    if (strncmp(token, ">", 1) == 0) {
      primary = TRUE;
    } else if (strncmp(token, "+", 1) == 0) {
      primary = FALSE;
    } else {
      die("Unexpected token value \"%s\", expected either \">\" or \"+\"\n", token);
    }
    free(token);
    godag->terms[i].id = next_token(buffer, fp, is_tab);
    consume(buffer, fp, is_tab);
    token = next_token(buffer, fp, is_newline);
    consume(buffer, fp, is_newline);
    node_index = parse_pos_int(token);
    free(token);
    if (node_index == 0) { //term is obsolete
      godag->terms[i].node = &(godag->root);
    } else {
      godag->terms[i].node = godag->nodes+(node_index-1);
      godag->terms[i].node->id = godag->terms[i].id;
    }
  }
  buf_destroy(buffer);
  fclose(fp);
  return godag;
}

/*
 * frees memory associated with the godag
 */
void destroy_go_dag(GODAG_T *godag) {
  int i;
  for (i = 0; i < godag->term_count; ++i) {
    free(godag->terms[i].id);
  }
  free(godag->terms);
  for (i = 0; i < godag->node_count; ++i) {
    free(godag->nodes[i].parents);
    free(godag->nodes[i].children);
    free(godag->nodes[i].name);
    free(godag->nodes[i].group);
  }
  free(godag->root.children);
  free(godag->nodes);
}

/*
 * compares a goid and a GOTERM_T 
 */
int compare(const void *key, const void *value) {
  char *key_str = (char*)key;
  GOTERM_T *value_term = (GOTERM_T*)value;
  return strcmp(key_str, value_term->id);
}

/*
 * finds the node for the goid
 */
GONODE_T* find(GODAG_T *godag, char *goid) {
  int pos;
  pos = binary_search(goid, godag->terms, godag->term_count, sizeof(GOTERM_T), compare);
  if (pos < 0) return NULL;
  return godag->terms[pos].node;
}

/*
 * recursively marks the node and its parents
 */
void recursive_mark(GONODE_T *node) {
  int i;
  if (!node->marked) {
    node->marked = TRUE;
    for (i = 0; i < node->pcount; ++i) {
      recursive_mark(node->parents[i]);
    }
  }
}

/*
 * marks the node in the tree (if it exists)
 */
void mark(GONODE_T *node) {
  if (node != NULL) recursive_mark(node);
}

/*
 * recursively clears child nodes
 */
void recursive_clear(GONODE_T *node) {
  int i;
  if (node->marked) {
    node->marked = FALSE;
    for (i = 0; i < node->ccount; ++i) {
      recursive_clear(node->children[i]);
    }
  }
}

/*
 * clears any marked nodes in the tree
 */
void clear(GODAG_T *godag) {
  recursive_clear(&(godag->root));
}

/*
 * get_type
 * returns  y for implied nodes
 *          n for non-implied
 *          u for unknown
 */
char get_type(GONODE_T *node) {
  int i;
  if (node == NULL) return 'u';
  if (node->marked) {
    for (i = 0; i < node->ccount; ++i) {
      if (node->children[i]->marked) 
        return 'y'; //marked child so not a leaf 
    }
    return 'n'; //marked and no marked children so a leaf
  }
  return 'y'; //not marked so not a leaf
}

/*
 * get_group
 */
char* get_group(GONODE_T *node) {
  if (node == NULL) {
    return "";
  }
  return node->group;
}


/*
 * get_nabove
 * gets the precalculated number of nodes above
 */
int get_nabove(GONODE_T *node) {
  if (node == NULL) {
    return 0;
  }
  return node->nabove;
}

/*
 * get_nbelow
 * gets the precalculated number of nodes below
 */
int get_nbelow(GONODE_T *node) {
  if (node == NULL) {
    return 0;
  }
  return node->nbelow;
}

/*
 * get_name
 * gets the name of the node
 */
char* get_name(GONODE_T *node) {
  if (node == NULL) {
    return "";
  }
  return node->name;
}

/*
 * read_until_match
 * Reads until the first match and returns a number indicating which,
 * if a match can't be found returns 0
 * If outfp is non-NULL then the read data is writen to it
 */
/*
int read_until_match(BUF_T *buffer, FILE *fp, FILE *outfp, MATCH_T **matches, int num) {
  int i, pos, remain, best, besti;
  BOOLEAN_T match;
  for (i = 0; i < num; ++i) {
    if (buf_capacity(buffer) < matches[i]->length) {
      die("Buffer too small to complete match\n");
    }
  }
  //try to get a match without doing any reading
  buf_mark(buffer);
  best = -1;
  besti = -1;
  match = FALSE;
  for (i = 0; i < num; ++i) { //find the match or possible match that starts at the earliest position
    buf_reset(buffer);
    remain = buf_kmp_search(buffer, matches[i]->string, matches[i]->table, matches[i]->length);
    pos = buf_position(buffer);
    if (pos < best || best == -1) {
      best = pos;
      besti = i;
      match = (remain == 0);
    }
  }
  //if outfp is set, write the non-matching data out
  if (outfp) {
    int limit = buf_limit(buffer);
    buf_reset(buffer);
    buf_set_limit(buffer, best);
    while (buf_remaining(buffer)) {
      buf_fwrite(buffer, outfp);
      if (ferror(outfp)) {
        die("Error occured while writing, error was given as %s\n", strerror(ferror(outfp)));
      }
    }
    buf_set_limit(buffer, limit);
  }
  buf_set_position(buffer, best);
  if (match) {//the earliest match was a complete match so we can return
    return besti + 1;
  }
  //ok so we need to read, prepare the buffer for reading
  buf_compact(buffer);
  while (!feof(fp)) {
    buf_fread(buffer, fp);
    if (ferror(fp)) {
      die("Error occured while reading, error was given as %s\n", strerror(ferror(fp)));
    }
    //have another go at finding the match 
    buf_flip(buffer);
    buf_mark(buffer);
    best = -1;
    besti = -1;
    match = FALSE;
    for (i = 0; i < num; ++i) { //find the match or possible match that starts at the earliest position
      buf_reset(buffer);
      remain = buf_kmp_search(buffer, matches[i]->string, matches[i]->table, matches[i]->length);
      pos = buf_position(buffer);
      if (pos < best || best == -1) {
        best = pos;
        besti = i;
        match = (remain == 0);
      }
    }
    //if outfp is set, write the non-matching data out
    if (outfp) {
      int limit = buf_limit(buffer);
      buf_reset(buffer);
      buf_set_limit(buffer, best);
      while (buf_remaining(buffer)) {
        buf_fwrite(buffer, outfp);
        if (ferror(outfp)) {
          die("Error occured while writing, error was given as %s\n", strerror(ferror(outfp)));
        }
      }
      buf_set_limit(buffer, limit);
    }
    buf_set_position(buffer, best);
    if (match) {//the earliest match was a complete match so we can return
      return besti + 1;
    }
    //maybe we can read more into the buffer?
    buf_compact(buffer);
  }
  buf_flip(buffer);
  //got to the end of the file
  buf_mark(buffer);
  best = -1;
  besti = -1;
  for (i = 0; i < num; ++i) { //find the match that starts at the earliest position as we have nothing more to read we ignore partial matches
    buf_reset(buffer);
    if (buf_kmp_search(buffer, matches[i]->string, matches[i]->table, matches[i]->length) == 0) {
      pos = buf_position(buffer);
      if (pos < best || best == -1) {
        best = pos;
        besti = i;
      }
    }
  }
  match = (best != -1);
  if (!match) best = buf_limit(buffer);
  //if outfp is set, write the non-matching data out
  if (outfp) {
    int limit = buf_limit(buffer);
    buf_reset(buffer);
    buf_set_limit(buffer, best);
    while (buf_remaining(buffer)) {
      buf_fwrite(buffer, outfp);
      if (ferror(outfp)) {
        die("Error occured while writing, error was given as %s\n", strerror(ferror(outfp)));
      }
    }
    buf_set_limit(buffer, limit);
  }
  buf_set_position(buffer, buf_limit(buffer));
  return (match) ? besti + 1 : 0;
}
*/

/*
 * calls buf_fread_until with one string to match
 */
int read_until_match1(BUF_T *buffer, FILE *fp, FILE *outfp, BMSTR_T *match1) {
  BMSTR_T *matches[1];
  matches[0] = match1;
  return buf_fread_until(buffer, fp, outfp, 1, matches);
}

/*
 * calls buf_fread_until with two strings to match
 */
int read_until_match2(BUF_T *buffer, FILE *fp, FILE *outfp, BMSTR_T *match1, BMSTR_T *match2) {
  BMSTR_T *matches[2];
  matches[0] = match1;
  matches[1] = match2;
  return buf_fread_until(buffer, fp, outfp, 2, matches);
}


/*
 * Rewrites the gomo xml out to a new temp file and then renames to overwrite the original
 */
void rewrite_gomo_xml(const char *tmp_dir, char *xml_file, GODAG_T *godag) {
  GONODE_T *node;
  char *tmp_file, *motifid, *termid, *desc;
  char implied;
  int tmpfd, which;
  FILE *tmpfp, *xmlfp, *outfp;
  struct stat stat_buf;
  long motif_start;
  BMSTR_T *goterm_strs[6];

  //create all search strings
  BMSTR_T *bm_motif_start   = bmstr_create("<motif");
  BMSTR_T *bm_id            = bmstr_create("id=\"");
  BMSTR_T *bm_goterm_start  = bmstr_create("<goterm");
  BMSTR_T *bm_group         = bmstr_create("group=\"\"");
  BMSTR_T *bm_nabove        = bmstr_create("nabove=\"\"");
  BMSTR_T *bm_nbelow        = bmstr_create("nbelow=\"\"");
  BMSTR_T *bm_implied       = bmstr_create("implied=\"u\"");
  BMSTR_T *bm_name          = bmstr_create("name=\"\"");
  BMSTR_T *bm_endtag        = bmstr_create(">");
  BMSTR_T *bm_endtag2       = bmstr_create("/>");
  BMSTR_T *bm_motif_end     = bmstr_create("</motif");

  goterm_strs[0] = bm_group;
  goterm_strs[1] = bm_nabove;
  goterm_strs[2] = bm_nbelow;
  goterm_strs[3] = bm_implied;
  goterm_strs[4] = bm_name;
  goterm_strs[5] = bm_endtag;
  //generate the temp file path
  if (tmp_dir == NULL) {
    tmp_dir = "/tmp/";
  }
  tmp_file = make_path_to_file(tmp_dir, "gomo_XXXXXX");
  //find out the file permissions on the xml file, we want to maintain them
  if (stat(xml_file, &stat_buf) == -1) {
    die("Unable to stat xml file \"%s\", error was given as: %s\n", xml_file, strerror(errno)); 
  }
  //open the xml file
  if ((xmlfp = fopen(xml_file, "r")) == NULL) {
    die("Error opening xml file \"%s\", error was given as: %s\n", xml_file, strerror(errno));
  }
  //open the temp file where we write out the changed xml file
  if ((tmpfd = mkstemp(tmp_file)) == -1) {
    die("Error creating temporary file, error was given as: %s\n", strerror(errno));
  }
  //set the file permissions the same as the xml file
  if (fchmod(tmpfd, stat_buf.st_mode) == -1) {
    die("Error setting file permissions on temporary file \"%s\" descriptor (%d), error was given as: %s\n", tmp_file, tmpfd, strerror(errno));
  }
  //open the file descriptor as a stream
  if ((tmpfp = fdopen(tmpfd, "w")) == NULL) {
    die("Error opening temporary file \"%s\" descriptor (%d) as a stream, error was given as: %s\n", tmp_file, tmpfd, strerror(errno));
  }

  BUF_T *buffer = buf_create(100);
  buf_flip(buffer); //puts the buffer into reading mode as read_util_match expects it
  //look for the start of a motif
  while (read_until_match1(buffer, xmlfp, tmpfp, bm_motif_start)) {
    //skip over '<motif'
    fprintf(tmpfp, "<motif");
    buf_skip(buffer, bmstr_length(bm_motif_start));
    //read until 'id="'
    if (!read_until_match1(buffer, xmlfp, tmpfp, bm_id)) die("Expected id attribute\n");
    //skip over 'id="'
    buf_skip(buffer, bmstr_length(bm_id));
    //get token before '"' which is the id of the motif
    motifid = next_token(buffer, xmlfp, is_quote);
    //skip over '"'
    buf_skip(buffer, 1);
    //output the id attribute to the temp file
    fprintf(tmpfp, "id=\"%s\"", motifid);
    //read until '>'
    which = read_until_match2(buffer, xmlfp, tmpfp, bm_endtag, bm_endtag2); 
    if (!which) die("Expected motif tag close bracket\n");
    if (which == 2) continue; //no content in motif tag
    //found '>' so skip it and then process the contents
    buf_skip(buffer, bmstr_length(bm_endtag));
    fprintf(tmpfp, ">");
    fseek(xmlfp, -buf_remaining(buffer), SEEK_CUR);
    buf_clear(buffer);
    buf_flip(buffer);
    //record the position returned from ftell
    motif_start = ftell(xmlfp);
    //now scan though the go terms marking them all the the DAG
    //read until '<goterm' or '</motif'
    while ((which = read_until_match2(buffer, xmlfp, NULL, bm_goterm_start, bm_motif_end)) <= 1) {
      if (which == 1) { //<goterm
        buf_skip(buffer, bmstr_length(bm_goterm_start));
        //read until 'id="'
        if (read_until_match2(buffer, xmlfp, NULL, bm_id, bm_endtag) != 1) die("Expected id attribute\n");
        //skip over 'id="'
        buf_skip(buffer, bmstr_length(bm_id));
        //get token before '"' which is the id of the go term 
        termid = next_token(buffer, xmlfp, is_quote);
        mark(find(godag, termid));
        free(termid);
      } else { //eof?
        die("Expected closing motif tag\n");
      }
    }
    // process again, and this time write out to file
    buf_skip(buffer, bmstr_length(bm_motif_end));
    fseek(xmlfp, motif_start, SEEK_SET);
    buf_clear(buffer);
    buf_flip(buffer);
    //read until '<goterm' or '</motif'
    while ((which = read_until_match2(buffer, xmlfp, tmpfp, bm_goterm_start, bm_motif_end)) <= 1) {
      if (which == 1) { //<goterm
        fprintf(tmpfp, "<goterm");
        buf_skip(buffer, bmstr_length(bm_goterm_start));
        //read until 'id="'
        if (read_until_match2(buffer, xmlfp, tmpfp, bm_id, bm_endtag) != 1) die("Expected id attribute\n");
        //skip over 'id="'
        buf_skip(buffer, bmstr_length(bm_id));
        //get token before '"' which is the id of the go term
        //this makes the assumption that the term id doesn't have any " characters in it (gomo makes the same assumption)
        termid = next_token(buffer, xmlfp, is_quote);
        //skip the '"'
        buf_skip(buffer, 1);
        //output the id attribute
        fprintf(tmpfp, "id=\"%s\"", termid);
        //lookup the node
        node = find(godag, termid);
        if (node == NULL) {
          fprintf(stderr, "Warning: node \"%s\" is unknown\n", termid);
        } else if (node == &(godag->root)) {
          fprintf(stderr, "Warning: node \"%s\" is obsolete\n", termid);
        }
        //making the assumption that the id will be before the implied or description attributes
        //read until 'group=""', 'nabove=""', 'nbelow=""', 'implied="u"', 'name=""', '>'
        while ((which = buf_fread_until(buffer, xmlfp, tmpfp, 6, goterm_strs)) <= 5) {
          switch (which) {
            case 1: //found 'group=""'
              buf_skip(buffer, bmstr_length(bm_group));
              fprintf(tmpfp, "group=\"%s\"", get_group(node));
              break;
            case 2: //found 'nabove=""'
              buf_skip(buffer, bmstr_length(bm_nabove));
              fprintf(tmpfp, "nabove=\"%d\"", get_nabove(node));
              break;
            case 3: //found 'nbelow=""'
              buf_skip(buffer, bmstr_length(bm_nbelow));
              fprintf(tmpfp, "nbelow=\"%d\"", get_nbelow(node)); 
              break;
            case 4://found 'implied="u"'
              buf_skip(buffer, bmstr_length(bm_implied));
              fprintf(tmpfp, "implied=\"%c\"", get_type(node));
              break;
            case 5://found 'name=""'
              buf_skip(buffer, bmstr_length(bm_name));
              fprintf(tmpfp, "name=\"%s\"", get_name(node));
              break;
            default://eof?
              die("Expected end of goterm tag\n");
          }
        }
        free(termid);
      } else { //eof?
        die("Expected closing motif tag\n");
      }
    }
    clear(godag);
    free(motifid);
  }
  //write out the remainder of the file
  buf_compact(buffer);
  while (!feof(xmlfp)) {
    buf_fread(buffer, xmlfp);
    if (ferror(xmlfp)) {
      die("Error occured while reading the file \"%s\", error was given as: %s\n", xml_file, strerror(ferror(xmlfp)));
    }
    buf_flip(buffer);
    buf_fwrite(buffer, tmpfp);
    if (ferror(tmpfp)) {
      die("Error occured while writing the file \"%s\", error was given as: %s\n", tmp_file, strerror(ferror(tmpfp)));
    }
    buf_compact(buffer);
  }
  buf_flip(buffer);
  while (buf_remaining(buffer)) {
    buf_fwrite(buffer, tmpfp);
    if (ferror(tmpfp)) {
      die("Error occured while writing the file \"%s\", error was given as: %s\n", tmp_file, strerror(ferror(tmpfp)));
    }
  }
  //close the files
  fclose(xmlfp);
  fclose(tmpfp);
  //replace xml file with the rewriten temp file
  if (rename(tmp_file, xml_file) == -1) {
    die("Error renaming temp file \"%s\" to \"%s\", error was given as: %s\n", tmp_file, xml_file, strerror(errno));
  }

  //clean up
  free(tmp_file);
}

#ifdef MAIN

/*
 * Reads a GO Tree from file
 * Sequentally reads and tweaks the gomo xml outputs to highlight the
 * most specific terms.
 */
int main(int argc, char** argv) {
  GODAG_T *godag;
  int i;
  if (argc < 2) {
    printf("USAGE:\n\tgomo_highlight <GO DAG> <GOMO XML>+\n");
    return 0;
  }
  godag = load_go_dag(argv[1]);
  for (i = 2; i < argc; ++i) {
    rewrite_gomo_xml(NULL, argv[i], godag);
  }

  //cleanup
  destroy_go_dag(godag);
  return 0;
}

#endif

