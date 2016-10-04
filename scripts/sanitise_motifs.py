'''
sanitise_motifs.py
===================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

The motifs geneset from the misgdb is full of redundant genesets i.e. there are multiple
genesets for transcription factors that differ by only a few genes/ bases in the motif.
This script aims to sanitise that geneset by taking the union of genesets for each transcription
factor. It uses the V$ nomenclature to define a gene set.


.. Overall purpose and function of the script>

Usage
-----

.. Example use case

Example::

   python sanitise_motifs.py

Type::

   python sanitise_motifs.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import collections
import itertools

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    tf_genes = collections.defaultdict(set)
    just_once = collections.defaultdict(int)
    
    ngenesets = set()
    nunknown = set()
    nmir = set()
    nredundant = set()
    
    for line in sys.stdin.readlines():
        data = line[:-1].split("\t")
        name, gene, ontology = data[2], data[1], data[0], 

        # add to genesets for counting
        ngenesets.add(name)
        
        # get rid of unknown
        if "UNKNOWN" in name:
            nunknown.add(name)
            continue
        # get rid of miRNAs at this point
        if len(name.split(",")) >= 1 and "MIR" in name:
            nmir.add(name)
            continue

        else:
            # collect all remaining
            # add transcription factor name here in tuple
            tf = name.split("$")[1].split("_")[0]
            tf_genes[(name, tf)].add(gene)

    for motif in tf_genes.keys():
        just_once[motif[1]] += 1

    just_once = [i for i in just_once if just_once[i] == 1]
    
    # iterate over sets and remove any redundant
    # ones
    keep = set()

    for name1, name2 in itertools.combinations([i[0] for i in tf_genes.keys()], 2):
        tf1 = name1.split("$")[1].split("_")[0]
        tf2 = name2.split("$")[1].split("_")[0]

        # collect those that only have one matrix present
        if tf1 in just_once:
            keep.add(name1)
        if tf2 in just_once:
            keep.add(name2)

        # check for redundancy
        if tf1 == tf2:
            if tf_genes[(name1, tf1)] == tf_genes[(name2, tf2)]:
                nredundant.add((name1, name2))

                # only keep one of the redundant sets
                keep.add(name1)
            else:
                keep.add(name1)
                keep.add(name2)

    for motif_id, genes in tf_genes.iteritems():
        name, tf = motif_id
        if name in keep:
            for gene in genes:
                sys.stdout.write("%(ontology)s\t%(gene)s\t%(name)s\t%(name)s\t%(ontology)s\n" % locals())

    E.info("started with %i genesets" % len(ngenesets))
    E.info("removed %i UNKNOWN" % len(nunknown))
    E.info("removed %i miRNA" % len(nmir))
    E.info("removed %i completely redundant" % len(nredundant))
    E.info("wrote %i motif genesets" % len(keep))

    
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
