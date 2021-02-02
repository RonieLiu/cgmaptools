#!/usr/bin/env python

"""
    cgmaptools - CGmapToRegion.py

    Copyright (C) Weilong Guo
    Contact: Weilong Guo <guoweilong@126.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
"""

import sys
#import os
#import os.path
import re

import gzip

#
def Get_key (str) :
    # This function should be consistent with the one in Sort_chr_pos.py
    match = re.match(r"^chr(\d+)", str, re.I)
    if match :
        chr = int(match.group(1))
    else :
        chr = str
    #
    return [chr, str]
    #re.match(r"^chr(\d+)", "chr019_random").group(1)
#
def CGmapToRegion (CGmap_fn, region_fn):
    try:
        if CGmap_fn :
            if CGmap_fn.endswith(".gz") :
                CGMAP = gzip.open(CGmap_fn, 'rb')
            else :
                CGMAP = open(CGmap_fn, 'r')
            #
        else :
            CGMAP = sys.stdin
        #
    except IOError:
        print ("\n[Error]:\n\t File cannot be open: " % CGmap_fn)
        exit(-1)
    #
    try:
        if region_fn.endswith(".gz") :
            REGION = gzip.open(region_fn, 'rb')
        else :
            REGION = open(region_fn, 'r')
        #
    except IOError:
        print ("\n[Error]:\n\t File cannot be open: %s" % region_fn)
        exit(-1)
    #
    line_cg_map = CGMAP.readline()
    line_region = REGION.readline()
    #
    mC_sum = 0  # cg_map文件中和region文件中染色体位置有overlap的且 NC_cg_map > 0的cg_map的行中 ( NmC_cg_map / NC_cg_map )  的和

    #cg_map_line_count的值是连续行符合这个条件的数量 chr_id_cg_map == chr_id_region && pos_start_region <= chr_cg_map <= pos_end_region && NC_cg_map > 0
    #  其实就是 cg_map染色体位置有overlap的且 NC_cg_map > 0的cg_map的行数
    cg_map_line_count = 0

    NmC_sum_cg_map = 0 # cg_map文件中和region文件中染色体位置有overlap的且 NC_cg_map > 0的cg_map的行中 NmC_cg_map  的和
    NC_sum_cg_map = 0 # cg_map文件中和region文件中染色体位置有overlap的且 NC_cg_map > 0的cg_map的行中 NC_cg_map  的和
    #
    try :
        chr_region, pos_start_region, pos_end_region = line_region.strip().split()[0:3]
    except ValueError :
        print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % region_fn)
        exit(-1)
    #
    while line_cg_map and line_region :
        #print "\t"+line_cg_map, "\t"+line_region
        #print "========"
        try :
            chr_cg_map, nuc_cg_map, pos_cg_map, pattern_cg_map, dinuc_cg_map, methyl_cg_map, NmC_cg_map, NC_cg_map = line_cg_map.strip().split()
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % CGmap_fn)
            exit(-1)
        #
        try :
            chr_region, pos_start_region, pos_end_region = line_region.strip().split()[0:3]
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % region_fn)
            exit(-1)
        #
        #print chr_region, pos_start_region, pos_end_region
        chr_id_cg_map = Get_key(chr_cg_map)
        chr_id_region = Get_key(chr_region)
        if chr_id_cg_map < chr_id_region :
            # cg_map中的chr比region中的chr更小，这一行的cg map不处理，继续读入下一行的cg map,
            line_cg_map = CGMAP.readline()
        elif chr_id_cg_map > chr_id_region :
            # count不再累加
            if cg_map_line_count > 0 :
                # 输出结果中的行
                print ("%s\t%s\t%s\t%.4f\t%d\t%.4f\t%d" % (chr_region, pos_start_region, pos_end_region, mC_sum/cg_map_line_count, cg_map_line_count, float(NmC_sum_cg_map)/NC_sum_cg_map, NC_sum_cg_map) )
            else :
                # 输出结果的行，count不大于0的，第4列和第6列的值是NA
                print ("%s\t%s\t%s\tNA\t%d\tNA\t%d" % (chr_region, pos_start_region, pos_end_region, cg_map_line_count, NC_sum_cg_map) )
            #
            # init
            mC_sum = 0 # # ( NmC_cg_map / NC_cg_map )  的和
            cg_map_line_count = 0
            NmC_sum_cg_map = 0
            NC_sum_cg_map = 0
            # do the next
            # 继续读入下一行的region, cg_map还在这一行
            line_region = REGION.readline()
        else : # chr_cg_map == chr_region
            if int(pos_cg_map) < int(pos_start_region) :
                line_cg_map = CGMAP.readline()
            elif int(pos_cg_map) > int(pos_end_region) :
                if cg_map_line_count > 0 :
                    # 输出结果中的行
                    print ("%s\t%s\t%s\t%.4f\t%d\t%.4f\t%d" % (chr_region, pos_start_region, pos_end_region, mC_sum/cg_map_line_count, cg_map_line_count, float(NmC_sum_cg_map)/NC_sum_cg_map, NC_sum_cg_map) )
                else :
                    # 输出结果的行，count不大于0的，第4列和第6列的值是NA
                    print ("%s\t%s\t%s\tNA\t%d\tNA\t%d" % (chr_region, pos_start_region, pos_end_region, cg_map_line_count, NC_sum_cg_map) )
                #
                # init
                mC_sum = 0
                cg_map_line_count = 0
                NmC_sum_cg_map = 0
                NC_sum_cg_map = 0
                line_region = REGION.readline()
            else : # pos_start_region <= chr_cg_map <= pos_end_region
                NC_cg_map = int(NC_cg_map)
                if NC_cg_map > 0 :
                    mC_sum = mC_sum + ( float(NmC_cg_map) / NC_cg_map )
                    #  chr_id_cg_map == chr_id_region  && pos_start_region <= chr_cg_map <= pos_end_region && NC_cg_map > 0
                    cg_map_line_count = cg_map_line_count + 1
                    NmC_sum_cg_map = NmC_sum_cg_map + int(NmC_cg_map)
                    NC_sum_cg_map = NC_sum_cg_map + int(NC_cg_map)
                #
                line_cg_map = CGMAP.readline()
            #
        #
    # End for reading files
    #
    while line_region :
        if cg_map_line_count > 0 :
            # 输出结果中的行
            print ("%s\t%s\t%s\t%.2f\t%d\t%.2f\t%d" % (chr_region, pos_start_region, pos_end_region, mC_sum/cg_map_line_count, cg_map_line_count, float(NmC_sum_cg_map)/NC_sum_cg_map, NC_sum_cg_map) )
        else :
            # 输出结果的行，count不大于0的，第4列和第6列的值是NA
            print ("%s\t%s\t%s\tNA\t%d\tNA\t%d" % (chr_region, pos_start_region, pos_end_region, cg_map_line_count, NC_sum_cg_map) )
        #
        # do the next
        line_region = REGION.readline()
    #
    REGION.close()
    if CGMAP is not sys.stdin:
        CGMAP.close()
    #
#
from optparse import OptionParser

# todo:
# * Provide the build-in sorting function

# ===========================================
def main():
    usage = "Usage: cgmaptools mtr [-i <CGmap>] -r <region> [-o <output>]\n" \
            "      (aka CGmapToRegion)\n" \
            "Description: Calculated the methylation levels in regions in two ways.\n" \
            "Notice: The region file should be sorted and non-overlapped regions.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2019-03-24\n" \
            "Format of Region file:\n" \
            "  #chr    start_pos  end_pos\n" \
            "   chr1   8275       8429\n" \
            "Output file format:\n" \
            "   #chr  start_pos  end_pos  mean(mC)  #_C  #read(C)/#read(T+C)  #read(T+C)\n" \
            "   chr1   8275       8429     0.34     72         0.40             164\n" \
            "Note: The two input CGmap files should be sorted by Sort_chr_pos.py first.\n" \
            "      This script would not distinguish CG/CHG/CHH contexts."
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmapFile",
                      help="File name end with .CGmap or .CGmap.gz. "
                           "If not specified, STDIN will be used.",
                      metavar="FILE")
    parser.add_option("-r", dest="regionFile",
                      help="Filename for region file, support *.gz", metavar="FILE")
    parser.add_option("-o", dest="outfile", default=None,
                      help="To standard output if not specified.\n")
    #
    (options, args) = parser.parse_args()
    #
    if (options.regionFile is None) :
        parser.print_help()
        exit(-1)
    #
    if (options.outfile is not None) :
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')
        #
    #
    CGmapToRegion(options.CGmapFile, options.regionFile)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
