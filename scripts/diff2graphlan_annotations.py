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
import cgatcore.iotools as IOTools
import collections
import cgatcore.experiment as E
from matplotlib import colors
import six
import random
import string

def readTree(infile, highest_level="family"):
    '''
    read the tree file and map nodes
    to use specified highest node.
    '''
    taxa = collections.defaultdict(set)
    inf = IOTools.open_file(infile)
    inf.readline()
    for line in inf.readlines():
        data = line.strip("\n")

        # use the length of the taxon name
        # to determine the taxonomic level
        # and add additional lower ranks
        # depending if they contain this string
        name_length = len(data.split("."))
        if name_length == 1:
            level = "kingdom"
        elif name_length == 2:
            level = "phylum"
        elif name_length == 3:
            level = "class"
        elif name_length == 4:
            level = "order"
        elif name_length == 5:
            level = "family"
        elif name_length == 6:
            level = "genus"
        elif name_length == 7:
            level = "species"
        else:
            try:
                line = line.replace("._", "_")
            except ValueError("malformatted line %s" % line):
                continue

        level_indices = {"kingdom": 0, "phylum": 1, "class": 2, "order": 3, "family": 4, "genus": 5, "species": 6}
        highest_level_index = level_indices[highest_level]
        if level_indices[level] == highest_level_index:
            taxa[data].add(data)
        elif level_indices[level] > highest_level_index:
            taxon = ".".join(data.split(".")[0:highest_level_index+1])
            taxa[taxon].add(data)
        else:
            continue
    return(taxa)

    
    # tree = collections.defaultdict(set)
    # inf = IOTools.open_file(infile)
    # # header line
    # inf.readline()
    # for line in inf.readlines():
    #     data = line.strip("\n").split(".")
    #     assert len(data) == 7, "malformatted taxon mapping file, must have exactly 7 columns"
    #     # columns with taxa
    #     k, p, c, o, f, g, s = data
    #     if level == "kingdom":
    #         l = 0
    #     elif level == "phylum":
    #         l = 1
    #     elif level == "class":
    #         l = 2
    #     elif level == "order":
    #         l = 3
    #     elif level == "family":
    #         l = 4
    #     elif level == "genus":
    #         l = 5
    #     elif level == "species":
    #         l = 6
    #     else:
    #         raise ValueError("level must be one of Kingdom, phylum, class, order, family")
    #     for t in [data[x] for x in range(l+1, len(data))]:
    #         tree[data[l]].add(t)
    # inf.close()
    # return (tree)
        

#############################################################
#############################################################
#############################################################

def writeInput(infile, tree):
    '''
    write the mapping file in the form that is taken
    by graphlan
    '''
    inf = IOTools.open_file(infile)
    # header line
    inf.readline()
    outf = IOTools.open_file("input.txt", "w")
    done = set()
    for line in inf.readlines():
        data = line.strip("\n").split("\t")
        # columns with taxa
        k, p, c, o, f, g, s = data
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
    get colours - lower clades inherit from "highest level"
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

    parser.add_option("--additional-labels", dest="additional_labels", type="string",
                      help="by default just the highest level labels are shown. Here you can add additional labels")
    
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

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
        assert len(list(tree.keys())) < 159, "not enough colours to support n = %i clades, please filter" % len(tree)
        new_tree = tree

    #get those to keep
    keep = set()
    for taxon, rest in new_tree.items():
        keep.add(taxon)
        for r in rest:
            keep.add(r)

    # get colours
    ncols = len(list(new_tree.keys()))
    colours = getColours(ncols)
    taxon2colour = {}
    for i in range(ncols):
        taxon2colour[list(new_tree.keys())[i]] = colours[i]

    # add in colours for all clade nodes that is
    # based on the highest node
    for h, r in new_tree.items():
        for taxon in r:
            taxon2colour[taxon] = taxon2colour[h]

    # read diff and output annotations
    result = collections.defaultdict(list)
    ps = []
    fcs =[]
    taxa = []
    colours = []
    shapes = []
    sig = []
    for line in options.stdin.readlines():

        # skip header
        if "taxa" and "log2FoldChange" in line:
            continue

        data = line.strip("\n").split("\t")
        taxon = data[0]
        if taxon not in keep:
            continue

        # assign colour
        if taxon in list(taxon2colour.keys()):
            colour = taxon2colour[taxon]
        else:
            colour = "NA"
        colours.append(colour)
            
        # append taxon to list
        taxa.append(taxon)
        
        # use -log10 p-value for clade size
        if options.use == "pval":
            column = 5
        elif options.use == "padj":
            column = 6
        else:
            raise ValueError("must use pval or padj, not %s for clade size" % options.use)

        # catch NA pvalues
        if data[column] == "NA":
            data[column] = 1
        
        p = -math.log10(float(data[column]))
        if p >= 1.3:
            shapes.append("*")
            sig.append(taxon)
        else:
            shapes.append("o")

        p = p*100
        result[taxon].append(p)
        ps.append(p)

        # fold changes
        if data[2] == "NA":
            fc = 0
        else:
            fc = data[2]
        fcs.append(fc)
        
    # output annotations
    options.stdout.write("%s\t%s\n" % ("clade_separation", "0.9"))
    for t, s in zip(taxa, shapes):
         options.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_shape", s))
    for t, c in taxon2colour.items():
         options.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_color", c))
    for t, c in taxon2colour.items():
        options.stdout.write("%s\t%s\t%s\n" % (t, "annotation_background_color", c))
    for t, c in taxon2colour.items():
        options.stdout.write("%s\t%s\t%s\n" % (t, "annotation_font_size", 12))

    for t, p, f in zip(taxa, ps, fcs):
        if t in sig and float(f) > 0:
            options.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_color", "r"))
            options.stdout.write("%s\t%s\t%f\n" % (t, "clade_marker_size", 200))
            options.stdout.write("%s\t%s\t%s\n" % (t, "annotation_background_color", "r"))
            options.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_height", 1, f))
            options.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_color", 1, "r"))
            options.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_alpha", 1, 0.5))

        elif t in sig and float(f) < 0:
            options.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_color", "b"))
            options.stdout.write("%s\t%s\t%f\n" % (t, "clade_marker_size", 200))
            options.stdout.write("%s\t%s\t%s\n" % (t, "annotation_background_color", "b"))
            options.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_height", 1, float(f)*-1))
            options.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_color", 1, "b"))
            options.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_alpha", 1, 0.5))
            
        elif t not in sig:
            options.stdout.write("%s\t%s\t%f\n" % (t, "clade_marker_size", p))

    # only output annotation for highest-level and
    # additional labels
    if options.additional_labels:
        additional_labels = options.additional_labels.split(",")
    else:
        additional_labels = []
    for t, p in zip(taxa, ps):
        #if t in list(new_tree.keys()):
        if t in additional_labels or t in list(tree.keys()):
           if "_" in t:
               a =  random.sample(list(string.ascii_lowercase),1)[0] + ":" + t
           else:
               a = t.split(".")[-1]
               options.stdout.write("%s\t%s\t%s\n" % (t, "annotation", a))

    # write the tree out
    outf = open("input_tree.txt", "w")
    for x, y in new_tree.items():
        for taxon in y:
            outf.write(taxon + "\n")
    outf.close()
            
           
    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
