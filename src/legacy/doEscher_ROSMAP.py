################################################################################### METADATA ####
#
# doEscher_ROSMAP.m
# Version: 2.0.1 (2020-10-20)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Visualisation of CellFie onto Escher maps.
#
# Citations:
#
# [1] Escher:
# King ZA  et al,  "Escher:  A web  application  for building,  sharing, and  embedding data-rich
# visualizations  of biological  pathways."  PLOS Computational  Biology, 11(8):  e1004321.  DOI:
# 10.1371/journal.pcbi.1004321.
#
# [2] Function readKeyValues:
# Ayesh                                      Alshukri,                                      2017,
# https://ayeshalshukri.co.uk/category/dev/python-script-to-extract-key-value-pairs-from-csv-file/
# Accessed 2020-11-02.
#
###################################################################################### LEGAL ####
#
# Copyright © 2020 King's College London
#
# This work is licensed under the Creative Commons Attribution 4.0 International Licence. To view
# a copy of this license,  visit http://creativecommons.org/licences/by/4.0/  or send a letter to
# Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
#
# Permission is hereby granted,  free of charge, to any person  obtaining a copy of this software
# and  associated  documentation  files  (the  "Software"),  to  deal  in  the  Software  without
# restriction,  including without  limitation the  rights to use,  copy, modify,  merge, publish,
# distribute, and/or sell copies  of the Software, and to permit persons  to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above  copyright notice  and this  permission notice  shall be included  with all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED  "AS IS", WITHOUT WARRANTY OF ANY KIND,  EXPRESS OR IMPLIED, INCLUDING
# BUT NOT  LIMITED TO THE WARRANTIES  OF MERCHANTABILITY,  FITNESS FOR A  PARTICULAR PURPOSE, AND
# NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS  OR COPYRIGHT HOLDERS BE  LIABLE FOR ANY CLAIM,
# DAMAGES, OR  OTHER LIABILITY,  WHETHER IN AN  ACTION OF CONTRACT,  TORT, OR  OTHERWISE, ARISING
# FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OF OR DEALING IN THE SOFTWARE.
#
##################################################################################### INPUTS ####
#
# [1] supplementary_ROSMAP/Results/CellFie/CellFie.csv
#     CellFie results.
#
#################################################################################### OUTPUTS ####
#
# [1] supplementary_ROSMAP/Results/Escher/
#     Escher results.
#
############################################################################# INITIALISATION ####
import escher
from escher import Builder
import cobra
import os

import csv

outpath = "supplementary_ROSMAP/Results/Escher/"
inpath = "supplementary_ROSMAP/Results/CellFie/CellFie.csv"
escherMapsDir = "Escher/Data/Escher maps/"

################################################################################## FUNCTIONS ####

def readKeyValues(inputCsvFile,keyColumn,ValueColumn):
    #input file reader
    infile = open(inputCsvFile, "r")
    read = csv.reader(infile)
    headers = next(read) # skip header

    returnDictionary={}
    returnList=[]

    #for each row
    for row in read:
        key   = row[keyColumn]
        value = row[ValueColumn]

        #Add to dictionary (note, will overwrite and only store single occurrences)
        returnDictionary[key] = value

        #Add to list (note, will store multiple occurrences)
        returnList.append([key,value])

    return returnDictionary
#\readKeyValues()

####################################################################################### MAIN ####
if not os.path.exists(outpath):
    os.mkdir(outpath)

infile = open(inpath, "r")
read = csv.reader(infile)
header = next(read)
infile.close()

for i in range(1, len(header)):
    filename = header[i]
    keyColumn    = 0
    ValueColumn  = i

    if not os.path.exists(outpath + filename):
        os.mkdir(outpath + filename)
    data = readKeyValues(inpath,keyColumn,ValueColumn)

    escher.rc['never_ask_before_quit'] = True

    for map in os.listdir(escherMapsDir):
        builder = Builder(map_json = escherMapsDir + map,
                          reaction_data = data,
                          reaction_scale_preset = "GaBuRd",
                          scroll_behavior = "zoom")
        builder.save_html(outpath + filename + "/" + filename + "_" + map + ".html")
