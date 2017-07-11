#include <Python.h>
#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include "bwamem.h"

static PyMethodDef module_functions[] = {
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) void init##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(bwapy)
{
    PyObject *m;

    MOD_DEF(m, "bwalib", "High-level binding to bwa mem",
            module_functions)

    if (m == NULL)
        return MOD_ERROR_VAL;
    return MOD_SUCCESS_VAL(m);

}

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#  ifdef MODULE_API_EXPORTS
#    define MODULE_API __declspec(dllexport)
#    define restrict __restrict
#  else
#    define MODULE_API __declspec(dllimport)
#  endif
#else
#  define MODULE_API
#endif

MODULE_API int module_init();

#ifdef __cplusplus
}
#endif


typedef struct { size_t n; mem_aln_t *aln; } mem_aln_v;


mem_aln_v *new_mem_aln_v (size_t n) {
	// Allocate a mem_aln_t vector
	mem_aln_v *alns = malloc (sizeof (mem_aln_v));
	if (alns == NULL)
		return NULL;
	alns->aln = malloc (n * sizeof (mem_aln_t));
	if (alns->aln == NULL) {
		free (alns);
		return NULL;
	}
	alns->n = n;
	return alns;
}

void free_mem_aln_v (mem_aln_v *alns) {
	// free mem_aln_v and all its submembers
	if (alns != NULL) {
		for (size_t i = 0; i < alns->n; ++i) {
			free(alns->aln[i].cigar);
		}
        	free(alns);
	}
}

bwaidx_t *bwa_idx_load_all(const char *hint){
	return bwa_idx_load(hint, BWA_IDX_ALL);
}

mem_aln_v *align(mem_opt_t * opt, bwaidx_t * idx, char * seq) {
	//mem_opt_t *opt;
        const size_t seq_len = strlen(seq);

	//opt = mem_opt_init(); // default values
	mem_alnreg_v ar;
	ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, seq_len, seq);

	// We won't take secondary alignments
	size_t primary = 0;
        for (size_t i = 0; i < ar.n; ++i) {
		if (ar.a[i].secondary >= 0) continue;
		primary++;
	}

        if(primary == 0){
		free(ar.a);
		return NULL;
	} else {
		mem_aln_v *alns = new_mem_aln_v(primary);
		alns->n = primary;

		primary = 0;
        	for (size_t i = 0; i < ar.n; ++i) {
			if (ar.a[i].secondary >= 0) continue;
			alns->aln[primary++] = mem_reg2aln(opt, idx->bns, idx->pac, seq_len, seq, &ar.a[i]);
		}

		free(ar.a);
		return alns;
	}
}
