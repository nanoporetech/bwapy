import argparse
from collections import namedtuple
import importlib
import imp
import os
import sys

from cffi import FFI
ffi = FFI()

"""High-level interface to bwa (mem) aligner."""


def get_shared_lib(name):
    """Cross-platform resolution of shared-object libraries, working
    around vagueries of setuptools.
    :param name: name of shared library to find.
    
    :returns: FFI shared library object.
    """
    try:
        # after 'python setup.py install' we should be able to do this
        lib_file = importlib.import_module(name).__file__
    except Exception as e:
        try:
            # after 'python setup.py develop' this should work
            lib_file = imp.find_module(name)[1]
        except Exception as e:
            raise ImportError('Cannot locate C library "{}".'.format(name))
        else:
            lib_file = os.path.abspath(lib_file)
    finally:
        library = ffi.dlopen(lib_file)
    return library


libbwa = get_shared_lib('bwalib')

ffi.cdef("""
  ////////////////////////////////
  // Alignment hit list structures
  //
  typedef struct {
    int64_t rb, re; // [rb,re): reference sequence in the alignment
    int qb, qe;     // [qb,qe): query sequence in the alignment
    int rid;        // reference seq ID
    int score;      // best local SW score
    int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
    int sub;        // 2nd best SW score
    int alt_sc;
    int csub;       // SW score of a tandem hit
    int sub_n;      // approximate number of suboptimal hits
    int w;          // actual band width used in extension
    int seedcov;    // length of regions coverged by seeds
    int secondary;  // index of the parent hit shadowing the current hit; <0 if primary
    int secondary_all;
    int seedlen0;   // length of the starting seed
    int n_comp:30, is_alt:2; // number of sub-alignments chained together
    float frac_rep;
    uint64_t hash;
  } mem_alnreg_t;

  typedef struct { size_t n, m; mem_alnreg_t *a; } mem_alnreg_v;

  typedef struct {   // This struct is only used for the convenience of API.
    int64_t pos;     // forward strand 5'-end mapping position
    int rid;         // reference sequence index in bntseq_t; <0 for unmapped
    int flag;        // extra flag
    uint32_t is_rev:1, is_alt:1, mapq:8, NM:22; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
    int n_cigar;     // number of CIGAR operations
    uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
    char *XA;        // alternative mappings
  
    int score, sub, alt_sc;
  } mem_aln_t;

  typedef struct { size_t n; mem_aln_t *aln; } mem_aln_v;

  void free_mem_aln_v (mem_aln_v *alns);

  ///////////////////////
  // bwa index structures
  //
  typedef uint64_t bwtint_t;

  typedef struct {
    bwtint_t primary; // S^{-1}(0), or the primary index of BWT
    bwtint_t L2[5]; // C(), cumulative count
    bwtint_t seq_len; // sequence length
    bwtint_t bwt_size; // size of bwt, about seq_len/4
    uint32_t *bwt; // BWT
    // occurance array, separated to two parts
    uint32_t cnt_table[256];
    // suffix array
    int sa_intv;
    bwtint_t n_sa;
    bwtint_t *sa;
  } bwt_t;

  typedef struct {
    int64_t offset;
    int32_t len;
    int32_t n_ambs;
    uint32_t gi;
    int32_t is_alt;
    char *name, *anno;
  } bntann1_t;

  typedef struct {
    int64_t offset;
    int32_t len;
    char amb;
  } bntamb1_t;

  typedef struct {
    int64_t l_pac;
    int32_t n_seqs;
    uint32_t seed;
    bntann1_t *anns; // n_seqs elements
    int32_t n_holes;
    bntamb1_t *ambs; // n_holes elements
    FILE *fp_pac;
  } bntseq_t;

  typedef struct {
    bwt_t    *bwt; // FM-index
    bntseq_t *bns; // information on the reference sequences
    uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base

    int    is_shm;
    int64_t l_mem;
    uint8_t  *mem;
  } bwaidx_t;

  bwaidx_t *bwa_idx_load_all(const char *hint); 
  void bwa_idx_destroy(bwaidx_t *idx);

  ///////////////////
  // Run an alignment
  //
  mem_aln_v *align(bwaidx_t * index, char * seq);
""")


Alignment = namedtuple('Alignment', [
    'rname', 'orient', 'pos', 'mapq', 'cigar', 'NM'
])

class BwaAligner(object):
    def __init__(self, index:str):
        self.index_base = index.encode()
        self.cigchar = "MIDSH"
        self.index = libbwa.bwa_idx_load_all(self.index_base)
        if self.index == ffi.NULL:
            raise ValueError('Failed to load bwa index.')

    def __del__(self):
        libbwa.bwa_idx_destroy(self.index)

    def _build_alignment(self, aln):
        cigar = aln.cigar
        cigar = ''.join(
            # oplen + op
            str(cigar[k]>>4) + self.cigchar[cigar[k] & 0xf]
            for k in range(aln.n_cigar)
        )
        return Alignment(
            ffi.string(self.index.bns.anns[aln.rid].name).decode(),
            '+-'[aln.is_rev], aln.pos, aln.mapq, cigar, aln.NM
        )


    def align_seq(self, seq:str):
        """Align a sequence to the index.

        :param seq: base sequence to align

        :returns: tuple of :class:`Alignment`
        """
        alns = libbwa.align(self.index, seq.encode())
        if alns == ffi.NULL:
            alignments = tuple()
        else:
            alignments = tuple(self._build_alignment(alns.aln[i]) for i in range(alns.n))
            libbwa.free_mem_aln_v(alns)
        return alignments


def get_parser():
    parser = argparse.ArgumentParser('Align a sequence with bwa mem.')
    parser.add_argument('index', help='bwa index base path.')
    parser.add_argument('sequence', nargs='+', help='base sequence')
    return parser


def main():
    args =  get_parser().parse_args()
    aligner = BwaAligner(args.index)
    for i, seq in enumerate(args.sequence, 1):
        alignments = aligner.align_seq(seq)
        print('Found {} alignments for input {}.'.format(len(alignments), i))
        for aln in alignments:
            print('  ', aln)
 
