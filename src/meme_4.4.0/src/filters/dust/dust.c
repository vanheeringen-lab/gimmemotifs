#include <math.h>
#include <ctype.h>

#define TRUE   1
#define FALSE  0

static int word = 3; 
static int window = 64; 
static int window2 = 32; 
static int level = 20;

static int mv, iv, jv;

void set_dust_level(value)
int value;
{
	level = value;
}

void set_dust_window(value)
int value;
{
	window = value;
	window2 = window / 2;
}

void set_dust_word(value)
int value;
{
	word = value;
}

static void wo1(len, s, ivv)
int len;
char *s;
int ivv;
{
	int i, ii, j, v, t, n, n1, sum;
	static int counts[32*32*32];
	static int iis[32*32*32];
	int js, nis;

	n = 32 * 32 * 32;
	n1 = n - 1;
	nis = 0;
	i = 0;
	ii = 0;
	sum = 0;
	v = 0;
	for (j=0; j < len; j++, s++) {
		ii <<= 5;
		if (isalpha(*s)) {
			if (islower(*s)) {
				ii |= *s - 'a';
			} else {
				ii |= *s - 'A';
			}
		} else {
			i = 0;
			continue;
		}
		ii &= n1;
		i++;
		if (i >= word) {
			for (js=0; js < nis && iis[js] != ii; js++) ;
			if (js == nis) {
				iis[nis] = ii;
				counts[ii] = 0;
				nis++;
			}
			if ((t = counts[ii]) > 0) {
				sum += t;
				v = 10 * sum / j;
				if (mv < v) {
					mv = v;
					iv = ivv;
					jv = j;
				}
			}
			counts[ii]++;
		}
	}
}

static int wo(len, s, beg, end)
int len;
char *s;
int *beg, *end;
{
	int i, j, l1, v;

	l1 = len - word + 1;
	if (l1 < 0) {
		*beg = 0;
		*end = len - 1;
		return 0;
	}
	mv = 0;
	iv = 0;
	jv = 0;
	for (i=0; i < l1; i++) {
		wo1(len-i, s+i, i);
	}
	*beg = iv;
	*end = iv + jv;
	return mv;
}

void dust(len, s)
int len;
char *s;
{
	int i, j, l, from, to, a, b, v;

	from = 0;
	to = -1;
	for (i=0; i < len; i += window2) {
		from -= window2;
		to -= window2;
		l = (len > i+window) ? window : len-i;
		v = wo(l, s+i, &a, &b);
		for (j = from; j <= to; j++) {
			if (isalpha(s[i+j]))
				s[i+j] = 'N';
		}
		if (v > level) {
			for (j = a; j <= b && j < window2; j++) {
				if (isalpha(s[i+j]))
					s[i+j] = 'N';
			}
			from = j;
			to = b;
		} else {
			from = 0;
			to = -1;
		}
	}
}

#define MAXREG    1001

typedef struct {
	int from;
	int to;
	int score;
} REGION;

static REGION reg[MAXREG];
static int nreg;

REGION *dust_segs(len, s)
int len;
char *s;
{
	int i, j, l, from, to, a, b, v;

	nreg = 0;
	from = 0;
	to = -1;
	for (i=0; i < len; i += window2) {
		l = (len > i+window) ? window : len-i;
		v = wo(l, s+i, &a, &b);
		if (v > level) {
			if (nreg > 0 && reg[nreg-1].to+1 >= a+i) {
				reg[nreg-1].to = b + i;
			} else if (nreg < MAXREG-1) {
				reg[nreg].from = a + i;
				reg[nreg].to = b + i;
				nreg++;
			}
			if (b < window2) {
				i += b - window2;
			}
		}
	}
	reg[nreg].from = 0;
	reg[nreg].to = -1;
	return reg;
}

