#!/home/users/astar/gis/hanjian/scratch-LORNAGS/hanjian/.conda/envs/PYTHON39/bin/python3.10
__doc__="""
    the script is used for extract the light(antisnese) strand
    """
__version__="v1.0"
__author__="hanjian"
__last_modify__="1-Feb-2025"

import argparse
import pysam
import numpy as np
import os

# ------------------------------------------
#PAIRED and PROPER_PAIR and ({READ1 and MREVERSE and not REVERSE} or {READ2 and REVERSE and not MREVERSE}) and not (UNMAP or MUNMAPorSECONDARY or QCFAIL or DUP or SUPPLEMENTARY).

def classify_strands(input_sam, output_light, output_heavy):
    """sense/antisense"""
    with pysam.AlignmentFile(input_sam, "r") as infile, \
         pysam.AlignmentFile(output_light, "w", template=infile) as light, \
         pysam.AlignmentFile(output_heavy, "w", template=infile) as heavy:

        for read in infile:
            if not (read.is_unmapped or read.mate_is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate or read.is_supplementary):

                if read.is_paired and read.is_proper_pair and ((read.is_read1 and read.mate_is_reverse and not read.is_reverse) or (read.is_read2 and read.is_reverse and not read.mate_is_reverse)):
                    light.write(read)
                else:
                    heavy.write(read)
        light.close()
        heavy.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i","--input",help="inputfile of mutrate.txt.gz file")
    parser.add_argument("-l","--output_light",default="-",help="outputfile of light file")
    parser.add_argument("-o","--output_heavy",default="-",help="outputfile of heavy file")
    args = parser.parse_args()
    classify_strands(args.input, args.output_light, args.output_heavy) 