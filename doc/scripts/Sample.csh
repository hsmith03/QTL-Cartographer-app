#!INSERTCSHCOMMAND
#           Sample.csh 
#   Copyright (C) 2000 Christopher J. Basten 
# Usage governed by the terms of the GNU  General Public License,  version 2 or higher
#  See the file COPYING in this directory
#
#   This file was meant as a simple example.  You  will need to edit it 
#   to work on your particular system with your data files.
#
#  Start by setting the variables needed.  
# 
set stem=sample                         #  filename stem
set mapfile=sample.mps                  #  result of running Mapmaker on sample.raw
set datafile=sample.raw                 #  mapmaker raw file
set qbin=/usr/local/bin                 #  where are the QTL Cart binaries
#
#   Now do the work...make sure you are in the same directory as the sample.mps and
#   sample.raw files
#
$qbin/Rmap -A -V -X $stem -i $mapfile   # convert the map
$qbin/Rcross -A -V -i $datafile         # convert the raw file
$qbin/Qstats -A -V                      # do basic stats
$qbin/LRmapqtl -A -V                    # single marker analysis
$qbin/SRmapqtl -A -V                    # stepwise regression
$qbin/Zmapqtl -A -V                     # interval mapping
$qbin/Eqtl -A -V  -o $stem.eqtl.3       # summarize best positions
$qbin/Zmapqtl -A -V -M 6                # do composite interval mappping
$qbin/Eqtl -A -V   -o $stem.eqtl.6      # summarize best model
