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

#define GAPCHAR '\0'

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

	// maximum no. of errors
	int k = error_rate * m;
	int last = k + 1;
	if (flags & START_WITHIN_SEQ1) {
		last = m;
	}
	// iterate over columns
	for (j = 1; j <= n; ++j) {
		// remember first entry
		Entry tmp_entry = column[0];

		// fill in first entry in this column TODO move out of loop
		if (flags & START_WITHIN_SEQ2) {
			column[0].cost = 0;
			column[0].origin = j;
			column[0].matches = 0;
		} else {
			column[0].cost = j * INSERTION_COST;
			column[0].origin = 0;
			column[0].matches = 0;
		}
		for (i = 1; i <= last; ++i) {
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
		while (column[last].cost > k) {
			last--;
		}
		if (last < m) {
			last++;
		} else {
			// found
			// if requested, find best match in last row
			if (flags & STOP_WITHIN_SEQ2) {
				// length of the aligned part of string1
				int length = m + min(column[m].origin, 0);
				int cost = column[m].cost;
				int matches = column[m].matches;
				if (cost <= length * error_rate && (matches > best_matches || (matches == best_matches && cost < best_cost))) {
					// update
					best_matches = matches;
					best_cost = cost;
					best_origin = column[m].origin;
					best_i = m;
					best_j = j;
				}
			}

		}
		// column finished
	}

	if (flags & STOP_WITHIN_SEQ1) {
		// search in last column // TODO last?
		for (i = 0; i <= m; ++i) {
			int length = i + min(column[i].origin, 0);
			int cost = column[i].cost;
			int matches = column[i].matches;
			if (cost <= length * error_rate && (matches > best_matches || (matches == best_matches && cost < best_cost))) {
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
