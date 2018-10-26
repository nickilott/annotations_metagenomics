'''
gmt2tsv.py
===========

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Convert .gmt file into .tsv file suitable for analysis with runGO.py.
.gmt files are usually downloaded from the MSigDatabase at the broad.

Usage
-----

Example::

   python gmt2tsv.py 

Type::

   python gmt2tsv.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse

import cgatcore.experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-o", "--ontology", dest="ontology", type="string",
                      help="ontology label"  )

    parser.add_option("-f", "--filter", dest="filter", action="store_true",
                      help="filter out genesets")

    parser.add_option("-l", "--filter-list", dest="filter_list", type="string",
                      help="list of pathways to keep"  )


    parser.set_defaults(ontology=None
                        , filter = False
                        , filter_list = None)

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.start( parser, argv = argv )

    if options.filter:
        assert options.filter_list, "must specify a list of pathways to keep"
        filter_set = set()
        for line in open(options.filter_list).readlines():
            filter_set.add(line[:-1])

    inf = options.stdin
    for line in inf.readlines():
        data = line[:-1].split("\t")
        name, description, evidence = data[0], data[0], data[1]
        if options.filter:
            if name not in filter_set: continue
        genes = data[2:]
        for gene in genes:
            options.stdout.write("\t".join([options.ontology, gene, name, description, evidence])+ "\n")
    
        

    ## write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
