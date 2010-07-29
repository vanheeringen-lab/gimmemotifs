#include <stdio.h>
#include <ctype.h>
#include <stdlib.h> /* added for exit */
#include <string.h> /* missing in original */

/* #include <malloc.h> also in stdlib.h */

#include "getfa.h"

#define NINC   56

#define TMPNAME  "tmp.seg"

void abo(mess)
char *mess;
{
	printf("Error: %s.\n", mess);
	exit(1);
}

void *mymalloc(size)
int size;
{
	void *buf;

	if ((buf = malloc(size)) == NULL) {
		abo("Not enough memory");
	}
	return buf;
}

FILE *myfopen(name, mode)
char *name, *mode;
{
	FILE *fp;

	if (name == NULL) {
		return (mode != NULL && *mode == 'r') ? stdin : stdout;
	}
	if ((fp = fopen(name, mode)) == NULL) {
		if (mode != NULL && *mode == 'r') {
			printf("No such file: \"%s\"\n", name);
		} else {
			printf("Failed to open file: \"%s\"\n", name);
		}
		exit(1);
	}
	return fp;
}

static void bad_format()
{
	abo("Not a FASTA file");
}

static void check_eof(c)
int c;
{
	if (c == EOF) {
		bad_format();
	}
}

FASTA *getfa(name)
char *name;
{
	int c, len, buflen;
	FILE *fp;
	char *head, *sq, *buf;
	FASTA *fa;

	fp = myfopen(name, "r");
	c = getc(fp);
	check_eof(c);
	if (c != '>') {
		bad_format();
	}
	len = 0;
	buflen = 1;
	head = mymalloc(1);
	while ((c = getc(fp)) != EOF && c != '\n') {
		if (len >= buflen-1) {
			buf = mymalloc(len + NINC);
			if (len > 0) {
				memcpy(buf, head, len);
				free(head);
			}
			head = buf;
			buflen = len + NINC;
		}
		head[len] = c;
		len++;
	}
	check_eof(c);
	head[len++] = '\0';
	len = 0;
	buflen = 1;
	sq = mymalloc(1);
	while ((c = getc(fp)) != EOF && c != '>') {
		if (isspace(c)) {
			continue;
		}
		if (len >= buflen-1) {
			buf = mymalloc(len + NINC);
			if (len > 0) {
				memcpy(buf, sq, len);
				free(sq);
			}
			sq = buf;
			buflen = len + NINC;
		}
		sq[len] = c;
		len++;
	}
	sq[len] = '\0';
	fa = mymalloc(sizeof(fa[0]));
	fa->len = len;
	fa->header = head;
	fa->seq = sq;
	if (fp != stdin) {
		fclose(fp);
	}
	return fa;
}

void putfa(fa, name)
FASTA *fa;
char *name;
{
	FILE *fp;
	char *s;
	int i, i60;
	size_t ignore; /* keep compiler happy */

	fp = myfopen(name, "w");
	s = (fa->header == NULL) ? "ANONYMOUS" : fa->header;
	if (*s != '>') {
		putc('>', fp);
	}
	ignore = fwrite(s, 1, strlen(s), fp); /* keep compiler happy */
	i60 = 50;
	for (i=0, s = fa->seq; i < fa->len; i++, s++, i60++) {
		if (i60 == 50) {
			putc('\n', fp);
			i60 = 0;
		}
		putc(*s, fp);
	}
	putc('\n', fp);
	if (fp != stdout) {
		fclose(fp);
	}
}

void getfafun(name, fun)
char *name;
void (*fun)();
{
	int c, len, headlen, sqlen, iseq;
	FILE *fp;
	char *head, *sq, *buf;
	FASTA *fa;

	fp = myfopen(name, "r");
	c = getc(fp);
	check_eof(c);
	if (c != '>') {
		bad_format();
	}
	headlen = 1;
	sqlen = 1;
	head = mymalloc(1);
	sq = mymalloc(1);
	fa = mymalloc(sizeof(fa[0]));
	iseq = 0;
	do {
		len = 0;
		while ((c = getc(fp)) != EOF && c != '\n') {
			if (len >= headlen-1) {
				buf = mymalloc(len + NINC);
				if (len > 0) {
					memcpy(buf, head, len);
					free(head);
				}
				head = buf;
				headlen = len + NINC;
			}
			head[len] = c;
			len++;
		}
		check_eof(c);
		head[len++] = '\0';
		len = 0;
		while ((c = getc(fp)) != EOF && c != '>') {
			if (isspace(c)) {
				continue;
			}
			if (len >= sqlen-1) {
				buf = mymalloc(len + NINC);
				if (len > 0) {
					memcpy(buf, sq, len);
					free(sq);
				}
				sq = buf;
				sqlen = len + NINC;
			}
			sq[len] = c;
			len++;
		}
		sq[len] = '\0';
		fa->len = len;
		fa->header = head;
		fa->seq = sq;
		fun(len, sq);
		putfa(fa, NULL);
		iseq++;
#if 0
		if (iseq % 10 == 0) {
			fprintf(stderr, "\r%d", iseq);
		}
#endif
	} while (c == '>');
#if 0
	fprintf(stderr, "\r%d sequences\n", iseq);
#endif
	free(head);
	free(sq);
	free(fa);
	if (fp != stdin) {
		fclose(fp);
	}
}

/*void main()
{
	FASTA *fa;

	fa = getfa(NULL);
	printf("Header: %s\n", fa->header);
	printf("Sequence (%d): %.*s\n", fa->len, fa->len, fa->seq);
	putfa(fa, NULL);
}  */

