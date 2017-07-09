Welcome to bwapy's documentation!
==================================

Python bindings to `bwa mem` aligner; sufficient to load and index and perform
alignments of sequences to the index to obtain basic statistics.

These python bindings are licensed under Mozilla Public License 2.0, bwa is licenced
under GNU General Public License v3.0.

Installation
------------

The git source repository contains bwa as a submodule. The repository should therefore
be cloned using the recursive option.

The package `setup.py` script requires `libbwa.a` to have been built in the submodule
directory before running. To build and install the package one should therefore run:

.. code-block:: bash

    git clone --recursive https://git/research/bwapy.git
    cd bwapy/bwa && make libbwa.a && cd ..
    python setup.py install


Performing Alignments
---------------------

The `BwaAligner` class provides a pythonic interface to `bwa mem` aligner. It
takes as input a bwa index fileset on construction and can then be used to find
alignments of sequences given as strings:

.. code-block:: python

    from bwapy import BwaAligner
    aligner = BwaAligner(args.index)
    alignments = aligner.align_seq(seq)
    print('Found {} alignments for input {}.'.format(len(alignments), i))
    for aln in alignments:
        print('  ', aln)

The alignments are returned as a named tuple, e.g.:

.. code-block:: python

    Alignment(rname='yeast', orient='+', pos=0, mapq=60, cigar='915M3D29M3D27M3D13M', NM=12)


Contents
--------

.. toctree::
   :maxdepth: 2


Full API reference
------------------

.. toctree::
   :maxdepth: 3
      
   bwapy

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

