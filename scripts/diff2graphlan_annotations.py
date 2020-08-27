'''
diff2graphlan_annotations.py
=============================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Create graphical parameters for input to graphlan.

.. Overall purpose and function of the script>

Usage
-----

.. Example use case

Example::

   python diff2graphlan_annotations.py

Type::

   python diff2graphlan_annotations.py --help

for command line help.

Command line options
--------------------

'''

import sys
import math
import CGAT.IOTools as IOTools
import collections
import CGAT.Experiment as E
from matplotlib import colors
import six
import random
import string

def readTree(infile, level="family"):
    '''
    read the tab separated tree file and map nodes
    to use specified highest node
    '''
    tree = collections.defaultdict(set)
    inf = IOTools.openFile(infile)
    # header line
    inf.readline()
    for line in inf.readlines():
        data = line.strip("\n").split("\t")
        # columns with taxa
        d, k, p, c, o, f, g, s = data
        if level == "domain":
            l = 0
        elif level == "kingdom":
            l = 0
        elif level == "phylum":
            l = 1
        elif level == "class":
            l = 2
        elif level == "order":
            l = 3
        elif level == "family":
            l = 4
        elif level == "genus":
            l = 5
        elif level == "species":
            l = 6
        else:
            raise ValueError("level must be one of domain, Kingdom, phylum, class, order, family")
        for t in [data[x] for x in range(l+1, len(data))]:
            tree[data[l]].add(t)
    inf.close()
    return tree
        

#############################################################
#############################################################
#############################################################

def writeInput(infile, tree):
    '''
    write the mapping file in the form that is taken
    by graphlan
    '''
    inf = IOTools.openFile(infile)
    # header line
    inf.readline()
    outf = IOTools.openFile("input.txt", "w")
    done = set()
    for line in inf.readlines():
        data = line.strip("\n").split("\t")
        # columns with taxa
        d, k, p, c, o, f, g, s = data
        new_string = ".".join([k, p, c, o, f, g, s])
        for taxon in tree.keys():
            if taxon in new_string:
                done.add(new_string)
    for record in done:
        outf.write(record + "\n")
    outf.close()
        
#############################################################
#############################################################
#############################################################

def getColours(ncols):
    '''
    get colours - lowe clades inherit from "highest level"
    clade
    '''
    colors_ = list(six.iteritems(colors.cnames))

    # Add the single letter colors.
    for name, rgb in six.iteritems(colors.ColorConverter.colors):
        hex_ = colors.rgb2hex(rgb)
        colors_.append((name, hex_))

    # Transform to hex color values.
    hex_ = [color[1] for color in colors_]
    cols = hex_[0:ncols]

    return cols
    
#############################################################
#############################################################
#############################################################
    

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-u", "--use", dest="use", type="choice",
                      choices=("pval", "padj"),
                      help="Type of p-value to use for clade size")

    parser.add_option("-m", "--taxa-map", dest="taxa_map", type="string",
                      help="Taxa mapping file - basically a text tree")

    parser.add_option("-l", "--highest-level", dest="highest_level", type="string",
                      help="highest taxonomic level to visualise")

    parser.add_option("-f", "--filter", dest="filter", action="store_true",
                      help="do you want to filter? will filter based on highest-level")

    parser.add_option("-k", "--keep", dest="keep", type="string",
                      help="keep all clades below these")
    
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # read tree
    tree = readTree(options.taxa_map, options.highest_level)

    # filter if neccessary
    new_tree = {}
    keep = set()
    if options.filter:
        assert options.keep, "must specify which clades to keep"
        to_keep = options.keep.split(",")
        for t in to_keep:
            new_tree[t] = tree[t]
    else:
        assert len(tree.keys()) < 159, "not enough colours to support n = %i clades, please filter" % len(tree)
        new_tree = tree

    # get those to keep
    keep = set()
    for taxon, rest in new_tree.iteritems():
        keep.add(taxon)
        for r in rest:
            keep.add(r)

    # get colours
    ncols = len(new_tree.keys())
    colours = getColours(ncols)
    taxon2colour = {}
    for i in range(ncols):
        taxon2colour[new_tree.keys()[i]] = colours[i]

    # read diff and output annotations
    result = collections.defaultdict(list)
    ps = []
    taxa = []
    colours = []
    shapes = []
    sig = []
    for line in options.stdin.readlines():

        # skip header
        if "taxa" and "logFC" in line:
            continue
        data = line.strip("\n").split("\t")
        taxon = data[-1].replace('"', '')
        if taxon not in keep:
            continue

        # assign colour
        if taxon in taxon2colour.keys():
            colour = taxon2colour[taxon]
        
        else:
            colour = "NA"
        colours.append(colour)
            
        # append taxon to list
        taxa.append(taxon)
        
        # use -log10 p-value for clade size
        if options.use == "pval":
            column = 3
        elif options.use == "padj":
            column = 4
        else:
            raise ValueError("must use pval or padj, not %s for clade size" % options.use)
        
        p = -math.log10(float(data[column]))
        if p >= 1.3:
            shapes.append("*")
            sig.append(taxon)
        else:
            shapes.append("o")

        p = p*100
        result[taxon].append(p)
        ps.append(p)
        
    # get maximum pvalue and make the rest a % of that
#    maxp = max(ps)
#    ps = [(x/maxp)*100 for x in ps]

    # output annotations
    for t, p in zip(taxa, ps):
        if t in sig:
            options.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_color", "r"))
            options.stdout.write("%s\t%s\t%f\n" % (t, "clade_marker_size", 200))
            options.stdout.write("%s\t%s\t%s\n" % (t, "annotation_background_color", "r"))
        else:
            options.stdout.write("%s\t%s\t%f\n" % (t, "clade_marker_size", p))

    for t, s in zip(taxa, shapes):
        options.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_shape", s))
    for t, c in taxon2colour.iteritems():
        options.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_color", c))
    for t, c in taxon2colour.iteritems():
        options.stdout.write("%s\t%s\t%s\n" % (t, "annotation_background_color", c))

    # only out put annotation for highest-level and
    # sig taxa
    for t, p in zip(taxa, ps):
        if t in sig or t in new_tree.keys():
            if "_" in t:
                a =  random.sample(list(string.lowercase),1)[0] + ":" + t
            else:
                a = t
            options.stdout.write("%s\t%s\t%s\n" % (t, "annotation", a))
            
    # output the file for input i.e. the tree for taxa we want
    writeInput(options.taxa_map, new_tree)
    
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
