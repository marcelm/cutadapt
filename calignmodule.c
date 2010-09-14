/*
Copyright (c) 2010 Marcel Martin <marcel.martin@tu-dortmund.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <Python.h>

#ifdef __GNUC__
#define UNUSED __attribute__((__unused__))
#else
#define UNUSED
#endif

/** Reverse a string in-place (helper function) */
static void reverse_string(char* s, int len) {
	int i;
	for (i = 0; i < len/2; ++i) {
		char tmp = s[i];
		s[i] = s[len-1-i];
		s[len-1-i] = tmp;
	}
}

PyDoc_STRVAR(semiglobalalign__doc__,
"semiglobalalign(string1, string2) -> (r1, r2, start1, stop1, start2, stop2, errors)\n\n\
\n\
Compute an end-gap free alignment (also called free-shift alignment or\n\
semiglobal alignment) of strings s1 and s2.\n\
\n\
The alignment uses the following score function:\n\
match    +1\n\
mismatch -1\n\
indel    -1\n\
\n\
Return a tuple (r1, r2, start1, stop1, start2, stop2, errors)\n\
where r1 and r2 are sequences of the same length containing the alignment\n\
(an INDEL is marked by '-').\n\
\n\
start1 is the position within r1 at which the part of s1, that is aligned, starts.\n\
start2 is the position within r1 at which the part of s1, that is aligned, ends.\n\
The same holds for start2, stop2.\n\
\n\
It is always the case that at least one of start1 and start2 is zero.\n\
\n\
It is always the case that either stop1==len(r1) or stop2==len(r2) or both\n\
(note that len(r1)==len(r2)). This is a property of semiglobal alignments.\n\
\n\
errors is the number of errors in the alignment.\n\
\n\
For example, semiglobalalign(\"SISSI\", \"MISSISSIPPI\") returns:\n\
\n\
r1 = [ '-', '-', '-', 'S', 'I', 'S', 'S', 'I', '-', '-', '-']\n\
r2 = [ 'M', 'I', 'S', 'S', 'I', 'S', 'S', 'I', 'P', 'P', 'I']\n\
start1, stop1 = 0, 5\n\
start2, stop2 = 3, 8\n\
errors = 0\n\
\n\
This corresponds to the following alignment:\n\
   SISSI\n\
   |||||\n\
MISSISSIPPI\n");

static PyObject *
py_semiglobalalign(PyObject *self UNUSED, PyObject *args)
{
	const char *s1;
	const char *s2;
	int m, n;

	if (!PyArg_ParseTuple(args, "s#s#", &s1, &m, &s2, &n))
		return NULL;

	/*
	DP Matrix:
	           s2 (j)
	        ----------> n
	       |
	s1 (i) |
	       |
	       V
	       m
	*/

	/*
	NOTE
	With C99, one could use this to access matrix elements:
	double (*matrix)[n];
	matrix = malloc ...

	matrix[i][j] = ...
	*/

	// direction constants for backtrace table
	enum { LEFT = 1, UP = 2, DIAG = 3 };

	// structure for a DP matrix entry
	typedef struct { int score; int backtrack; } Entry;

	// the DP matrix is stored column-major
	Entry* columns;
	columns = (Entry*)malloc((m+1)*(n+1)*sizeof(Entry));
	if (columns == NULL)
		return NULL;

	int i, j;

	// initialize first column
	for (i = 0; i <= m; ++i) {
		columns[i].score = 0;
		columns[i].backtrack = UP; // TODO never read
	}
	// initialize first row
	for (j = 0; j <= n; ++j) {
		columns[j*(m+1)].score = 0;
		columns[j*(m+1)].backtrack = LEFT; // TODO never read
	}

	// fill DP matrix
	// outer loop goes over columns
	Entry* cur_column;
	Entry* prev_column = columns;
	for (j = 1; j <= n; ++j) {
		cur_column = columns + j*(m+1);
		for (i = 1; i <= m; ++i) {
			int bt = DIAG;
			int score = prev_column[i-1].score + ((s1[i-1] == s2[j-1]) ? 1 : -1);
			int tmp = cur_column[i-1].score - 1;
			if (tmp > score) {
				bt = UP;
				score = tmp;
			}
			tmp = prev_column[i].score - 1;
			if (tmp > score) {
				bt = LEFT;
				score = tmp;
			}
			cur_column[i].score = score;
			cur_column[i].backtrack = bt;
		}
		prev_column = cur_column;
	}

	// find position with highest score in last column or last row
	int best = -1, best_j = -1, best_i = -1;
	for (j = 0; j <= n; ++j) { // last row
		if (columns[j*(m+1)+m].score >= best) {
			best_i = m;
			best_j = j;
			best = columns[j*(m+1)+m].score;
		}
	}
	assert(best_j != -1);
	Entry* last_column = columns + n*(m+1);
	for (i = 0; i <= m; ++i) {
		if (last_column[i].score >= best) {
			best_i = i;
			best_j = n;
			best = last_column[i].score;
		}
	}
	assert(best_i != -1);

	// trace back
	char* alignment1 = malloc((m+n+4)*sizeof(char));
	char* alignment2 = malloc((m+n+4)*sizeof(char));

	char* p1 = alignment1;
	char* p2 = alignment2;

	i = m;
	j = n;

	// first, walk from the lower right corner to the
	// position where we found the maximum score
	if (i == best_i) { // we are in the last row
		while (j > best_j) {
			*p1++ = '-';
			*p2++ = s2[--j];
		}
	}
	else { // we are in the last column
		while (i > best_i) {
			*p1++ = s1[--i];
			*p2++ = '-';
		}
	}

	int errors = 0;

	// the actual backtracing.
	// we build reverse sequences while backtracing and
	// reverse them afterwards.
	while (i > 0 && j > 0) {
		int direction = columns[j*(m+1)+i].backtrack;
		if (direction == DIAG) {
			if (s1[--i] != s2[--j])
				errors++;
			*p1++ = s1[i];
			*p2++ = s2[j];
		} else if (direction == LEFT) {
			errors++;
			*p1++ = '-';
			*p2++ = s2[--j];
		} else if (direction == UP) {
			*p1++ = s1[--i];
			*p2++ = '-';
			errors++;
		}
	}

	free(columns);

	int start1 = i;
	int start2 = j;


	while (j > 0) {
		*p1++ = '-';
		*p2++ = s2[--j];
	}
	while (i > 0) {
		*p1++ = s1[--i];
		*p2++ = '-';
	}
	assert(i == 0 && j == 0);

	assert(columns[best_j*(m+1)+best_i].score == length - 2*errors);

	// reverse result
	reverse_string(alignment1, p1-alignment1);
	reverse_string(alignment2, p2-alignment2);
	*p1 = '\0';
	*p2 = '\0';

	//return (r1, r2, start1, stop1, start2, stop2, errors)
	PyObject* o = Py_BuildValue("s#s#iiiii", alignment1, p1-alignment1, alignment2, p2-alignment2,
		start1, best_i, start2, best_j, errors);
	free(alignment1);
	free(alignment2);

	return o;
}

/* module initialization */

static PyMethodDef methods[] = {
	{ "semiglobalalign", py_semiglobalalign, METH_VARARGS, semiglobalalign__doc__ },
	{NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC
initcalign(void)
{
    (void) Py_InitModule("calign", methods);
}
