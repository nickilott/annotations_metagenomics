'''
eggnog2categories.py
=====================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Output a file that maps eggNOG identifiers to functional categories.

Usage
-----

Example::

   python eggnog2categories.py 

Type::

   python eggnog2categories.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.IOTools as IOTools


############################################
############################################
############################################


def readCategories(infile):
    '''
    Return a dictionary that maps a Category ID:Description
    '''
    result = {}
    inf = IOTools.openFile(infile)
    for line in inf.readlines():
        # category ids are stored in [] in the
        # input file
        if "[" in line:
            identifier = line[:-1].split("[")[1].split("]")[0]
            description = line[:-1].split("] ")[1]
            result[identifier] = description
        else:
            continue
    for x, y in result.iteritems():
        print x, y
    return result


############################################
############################################
############################################


def readNogs(infile):
    '''
    return a dictionary of NOG identifiers 
    mapping to category identifers
    '''
    result = {}
    inf = IOTools.openFile(infile)
    for line in inf.readlines():
        data = line[:-1].split("\t")
        nog, identifier = data
        result[nog] = identifier
    return result


############################################
############################################
############################################


def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-f", "--functions", dest="functions", type="string",
                      help="supply eggNOG functional categories i.e. ID:description file"  )
    parser.add_option("-n", "--nogs", dest="nogs", type="string",
                      help="supply nog2category mapping file"  )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    E.info("Reading functional categories")
    category2description = readCategories(options.functions)

    E.info("Reading NOG 2 category mapping")
    nog2category = readNogs(options.nogs)

    E.info("Assigning NOGs to functional category descriptions")
    for nog, category in nog2category.iteritems():
        # if the NOG does not have a category it
        # goes into function unknown (category S)
        if category == '':
            category = "S"
        # catch those with multiple annotations
        if len(category) > 1:
            for c in list(category):
                description = category2description[c]
                sys.stdout.write("NOG\t%(nog)s\t%(c)s\t%(description)s\tNA\n" % locals())
        else:
            description = category2description[category]
            sys.stdout.write("NOG\t%(nog)s\t%(category)s\t%(description)s\tNA\n" % locals())

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
