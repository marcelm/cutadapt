/*
Copyright (c) 2010,2011 Marcel Martin <marcel.martin@tu-dortmund.de>

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

enum FLAG {
	START_WITHIN_SEQ1 = 1,
	START_WITHIN_SEQ2 = 2,
	STOP_WITHIN_SEQ1 = 4,
	STOP_WITHIN_SEQ2 = 8,
	SEMIGLOBAL = 15,
	ALLOW_WILDCARD_SEQ1 = 1,
	ALLOW_WILDCARD_SEQ2 = 2,
};

// insertion means: inserted into seq1 (does not appear in seq2)
#define SCORE_MATCH 1
#define SCORE_MISMATCH -2
#define SCORE_DELETION -2
#define SCORE_INSERTION -2

/** Reverse a string in-place (helper function) */
static void reverse_string(char* s, int len) {
	int i;
	for (i = 0; i < len/2; ++i) {
		char tmp = s[i];
		s[i] = s[len-1-i];
		s[len-1-i] = tmp;
	}
}

// TODO fix docstring!

PyDoc_STRVAR(globalalign__doc__,
"globalalign(string1, string2, flags=0) -> (r1, r2, start1, stop1, start2, stop2, errors)\n\n\
\n\
Compute an end-gap free alignment (also called free-shift alignment or\n\
semiglobal alignment) of strings s1 and s2.\n\
\n\
The alignment uses the following score function:\n\
match    +1\n\
mismatch -1\n\
indel    -2\n\
\n\
Return a tuple (r1, r2, start1, stop1, start2, stop2, errors)\n\
where r1 and r2 are strings of the same length containing the alignment\n\
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
py_globalalign(PyObject *self UNUSED, PyObject *args)
{
	const char *s1;
	const char *s2;
	int m, n;
	int flags = 0;

	if (!PyArg_ParseTuple(args, "s#s#|i", &s1, &m, &s2, &n, &flags))
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
	typedef struct { int score; int backtrace; } Entry;

	// the DP matrix is stored column-major
	Entry* columns;
	columns = (Entry*)malloc((m+1)*(n+1)*sizeof(Entry));
	if (columns == NULL)
		return NULL;

	int i, j;

	// initialize first column
	for (i = 0; i <= m; ++i) {
		columns[i].score = (flags & START_WITHIN_SEQ1) ? 0 : i * SCORE_DELETION;
		columns[i].backtrace = UP;
	}

	// initialize first row
	for (j = 0; j <= n; ++j) {
		columns[j*(m+1)].score = (flags & START_WITHIN_SEQ2) ? 0 : j * SCORE_INSERTION;
		columns[j*(m+1)].backtrace = LEFT;
	}

	// fill the entire DP matrix
	// outer loop goes over columns
	Entry* cur_column;
	Entry* prev_column = columns;
	for (j = 1; j <= n; ++j) {
		cur_column = columns + j*(m+1);
		for (i = 1; i <= m; ++i) {
			int bt = DIAG;
			int score = prev_column[i-1].score + ((s1[i-1] == s2[j-1]) ? SCORE_MATCH : SCORE_MISMATCH);
			int tmp = cur_column[i-1].score + SCORE_INSERTION;
			if (tmp > score) {
				bt = UP;
				score = tmp;
			}
			tmp = prev_column[i].score + SCORE_DELETION;
			if (tmp > score) {
				bt = LEFT;
				score = tmp;
			}
			cur_column[i].score = score;
			cur_column[i].backtrace = bt;
		}
		prev_column = cur_column;
	}

	// initialize best score and its position to the bottomright cell
	int best_i = m; // also: s1stop
	int best_j = n; // also: s2stop
	int best = columns[(n+1)*(m+1)-1].score;

	if (flags & STOP_WITHIN_SEQ2) {
		// search also in last row
		for (j = 0; j <= n; ++j) {
			if (columns[j*(m+1)+m].score >= best) {
				best = columns[j*(m+1)+m].score;
				best_i = m;
				best_j = j;
			}
		}
	}
	if (flags & STOP_WITHIN_SEQ1) {
		// search also in last column
		Entry* last_column = columns + n*(m+1);
		for (i = 0; i <= m; ++i) {
			if (last_column[i].score >= best) {
				best_i = i;
				best_j = n;
				best = last_column[i].score;
			}
		}
	}

/*
	printf("s1: %s\ns2: %s\n", s1, s2);
	printf("m: %d n: %d\n", m, n);
	for (i = 0; i <= m; ++i) {
		printf("i=%5d   %c", i, (i>0)?s1[i-1]:'X');
		for (j = 0; j <= n; ++j) {
			printf("%3d ", columns[j*(m+1) + i].score); //, columns[j*(m+1) + i].backtrace);
		}
		printf("\n");
	}
	printf("best: (%d, %d)\n", best_i, best_j);
*/

	// trace back
	char* alignment1 = malloc((m+n+4)*sizeof(char));
	char* alignment2 = malloc((m+n+4)*sizeof(char));

	char* p1 = alignment1;
	char* p2 = alignment2;

	i = m;
	j = n;

	// first, walk from the lower right corner to the
	// position where we found the maximum score

	int errors = 0;

/*	printf("best_i: %d   best_j: %d\n", best_i, best_j);
	printf("i: %d   j: %d\n", i, j);*/

	int gaps_are_errors; // if gaps are currently errors, this is 1, 0 otherwise
	gaps_are_errors = (flags & STOP_WITHIN_SEQ2) ? 0 : 1;
	if (i == best_i) { // we are in the last row
		while (j > best_j) {
			*p1++ = '-';
			*p2++ = s2[--j];
			errors += gaps_are_errors;
		}
	}
	else { // we are in the last column
		gaps_are_errors = (flags & STOP_WITHIN_SEQ1) ? 0 : 1;
		while (i > best_i) {
			*p1++ = s1[--i];
			*p2++ = '-';
			errors += gaps_are_errors;
		}
	}
// 	printf("i: %d   j: %d\n", i, j);

	assert(i == best_i && j == best_j);

	// the actual backtracing
	// The alignments are constructed in reverse
	// and this is undone afterwards.
	while (i > 0 && j > 0) {
		int direction = columns[j*(m+1)+i].backtrace;
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

	int start1 = (flags & START_WITHIN_SEQ1) ? i : 0;
	int start2 = (flags & START_WITHIN_SEQ2) ? j : 0;

	errors += (i - start1) + (j - start2);

	while (j > 0) {
		*p1++ = '-';
		*p2++ = s2[--j];
	}
	while (i > 0) {
		*p1++ = s1[--i];
		*p2++ = '-';
	}
	assert(i == 0 && j == 0);

	//assert(columns[best_j*(m+1)+best_i].score == /*length*/ - 2*errors);

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

#define DELETION_COST 1
#define INSERTION_COST 1
#define MISMATCH_COST 1
#define MATCH_COST 0

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))


PyDoc_STRVAR(globalalign_locate__doc__,
"globalalign_locate(string1, string2, max_error_rate, flags=SEMIGLOBAL) -> (start1, stop1, start2, stop2, matches, errors)\n\n\
\n\
Locate one string within another by computing an optimal semiglobal alignment between string1 and string2.\n\
\n\
The alignment uses unit costs, which means that mismatches, insertions and deletions are\n\
counted as one error.\n\
\n\
The alignment is semiglobal, which means that an arbitrary number of characters in the beginning and end of\n\
string1 and in the beginning and end of string2 may be skipped at no cost. These\n\
skipped parts are described with two interval (start1, stop1), (start2, stop2).\n\
For example, the semiglobal alignment of SISSI and MISSISSIPPI looks like this:\n\
\n\
---SISSI---\n\
MISSISSIPPI\n\
\n\
start1, stop1 = 0, 5\n\
start2, stop2 = 3, 8\n\
(with zero errors)\n\
\n\
The aligned parts are string1[start1:stop1] and string2[start2:stop2].\n\
\n\
The alignment itself is not returned, only the tuple\n\
(start1, stop1, start2, stop2, matches, errors), where the first four fields have the\n\
meaning as described, matches is the number of matches and errors is the number of \n\
errors in the alignment.\n\
\n\
The error_rate is: errors / length where length is (stop1 - start1).\n\
\n\
(TODO length is computed on string1 only! It could also be min(or max)(stop1-start1, stop2-start2).)\n\
\n\
An optimal alignment fulfills all of these criteria:\n\
- error_rate <= max_error_rate\n\
- Among those alignments with error_rate <= max_error_rate, the alignment contains\n\
  a maximal number of matches (there is no alignment with more matches).\n\
- If there are multiple alignments with the same no. of matches, then one that\n\
  has minimal no. of errors is chosen.\n\
\n\
It is always the case that at least one of start1 and start2 is zero.\n\
\n\
The flags parameter allows to compute semiglobal alignments in which initial\n\
or trailing gaps in only one of the strings are penalized.\n");

static PyObject *
py_globalalign_locate(PyObject *self UNUSED, PyObject *args)
{
	const char *s1;
	const char *s2;
	int m, n;
	int flags = SEMIGLOBAL;
	double error_rate;
	int degenerate = 0;

	if (!PyArg_ParseTuple(args, "s#s#d|ii", &s1, &m, &s2, &n, &error_rate, &flags, &degenerate))
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

	// structure for a DP matrix entry
	typedef struct {
		int cost;
		int matches; // no. of matches in this alignment
		int origin; // where the alignment originated: negative for positions within seq1, positive for pos. within seq2
	} Entry;

	// only a single column of the DP matrix is stored
	Entry* column;
	column = (Entry*)malloc((m+1)*sizeof(Entry));
	if (column == NULL)
		return NULL;

	int i, j, best_i, best_j, best_cost, best_matches, best_origin;

	//initialize first column
	for (i = 0; i <= m; ++i) {
		column[i].matches = 0;
		column[i].cost = (flags & START_WITHIN_SEQ1) ? 0 : i * DELETION_COST;
		column[i].origin = (flags & START_WITHIN_SEQ1) ? -i : 0;
	}

	best_i = m;
	best_j = 0;
	best_cost = column[m].cost;
	best_matches = 0;
	best_origin = column[m].origin;

	// iterate over columns
	for (j = 1; j <= n; ++j) {
		// remember first entry
		Entry tmp_entry = column[0];

		// fill in first entry in this column
		if (flags & START_WITHIN_SEQ2) {
			column[0].cost = 0;
			column[0].origin = j;
			column[0].matches = 0;
		} else {
			column[0].cost = j * INSERTION_COST;
			column[0].origin = 0;
			column[0].matches = 0;
		}
		for (i = 1; i <= m; ++i) {
			int match = (s1[i-1] == s2[j-1])
						|| ((degenerate & ALLOW_WILDCARD_SEQ1) && (s1[i-1] == 'N'))
						|| ((degenerate & ALLOW_WILDCARD_SEQ2) && (s2[j-1] == 'N'));
			int cost_diag = tmp_entry.cost + (match ? MATCH_COST : MISMATCH_COST);
			int cost_deletion = column[i].cost + DELETION_COST;
			int cost_insertion = column[i-1].cost + INSERTION_COST;

			int origin, cost, matches;
			if (cost_diag <= cost_deletion && cost_diag <= cost_insertion) {
				// MATCH or MISMATCH
				cost = cost_diag;
				origin = tmp_entry.origin;
				matches = tmp_entry.matches + match;
			} else if (cost_insertion <= cost_deletion) {
				// INSERTION
				cost = cost_insertion;
				origin = column[i-1].origin;
				matches = column[i-1].matches;
			} else {
				// DELETION
				cost = cost_deletion;
				origin = column[i].origin;
				matches = column[i].matches;
			}

			// remember current cell for next iteration
			tmp_entry = column[i];

			column[i].cost = cost;
			column[i].origin = origin;
			column[i].matches = matches;
		}
		// column finished

		// if requested, find best match in last row
		if (flags & STOP_WITHIN_SEQ2) {
			// length of the aligned part of string1
			int length = m + min(column[m].origin, 0);
			int cost = column[m].cost;
			int matches = column[m].matches;
			if (cost <= length * error_rate && (matches > best_matches || (matches == best_matches && cost <= best_cost))) {
				// update
				best_matches = matches;
				best_cost = cost;
				best_origin = column[m].origin;
				best_i = m;
				best_j = j;
			}
		}
	}

	if (flags & STOP_WITHIN_SEQ1) {
		// search in last column
		for (i = 0; i <= m; ++i) {
			int length = i + min(column[i].origin, 0);
			int cost = column[i].cost;
			int matches = column[i].matches;
			if (cost <= length * error_rate && (matches > best_matches || (matches == best_matches && cost <= best_cost))) {
				// update best
				best_matches = matches;
				best_cost = cost;
				best_origin = column[i].origin;
				best_i = i;
				best_j = n;
			}
		}
	}

	free(column);
	int start1, start2;
	if (best_origin >= 0) {
		start1 = 0;
		start2 = best_origin;
	} else {
		start1 = -best_origin;
		start2 = 0;
	}

	// return (start1, stop1, start2, stop2, matches, errors)
	PyObject* o = Py_BuildValue("iiiiii", start1, best_i, start2, best_j, best_matches, best_cost);
	return o;
}

/* module initialization */

static PyMethodDef methods[] = {
	{"globalalign", (PyCFunction)py_globalalign, METH_VARARGS, globalalign__doc__},
	{"globalalign_locate", (PyCFunction)py_globalalign_locate, METH_VARARGS, globalalign_locate__doc__},
	{NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"calign",
	NULL, // module docstring
	-1,   // additional memory (-1 if not needed)
	methods,
	NULL,
	NULL, // GC traversal function
	NULL, // GC clear function
	NULL
};

#define INITERROR return NULL

PyObject *
PyInit_calign(void)

#else
#define INITERROR return

void
initcalign(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
	PyObject *module = PyModule_Create(&moduledef);
#else
	PyObject *module = Py_InitModule("calign", methods);
#endif

	if (module == NULL)
		INITERROR;

#if PY_MAJOR_VERSION >= 3
	return module;
#endif
}
