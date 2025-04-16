""" 
extractBarcodeReads.py

Dane Vassiliadis. Peter Mac. Jan 2020

Extracts reads of interest from fastq files using regex

From a fastq input, filter out reads that meet the following criteria:
- Exact match to the specified 5' and 3' constant regions
- [OPTIONAL] A correct sample index
- Average phred quality >= 30 across barcode length. Can be changed with --minqual
- No ambiguous N residues
- Barcode length at least 30. Can be changed with --minlength option

Optionally prefix barcodes with the sample index. 
Specify --index as tab delimited file with columns ID and Sequences. e.g.
ID	            Sequences
Sample_1	    ATCACG

"""

#------------------
# Load packages
#------------------

import re
import sys
import time
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from datetime import datetime
import gzip
import shutil
from optparse import OptionParser


startTime = datetime.now()

#------------------
# Define options
#------------------
usage = "USAGE: filter_trim_barcode_reads.py -i [input-fastq] -o [output-fastq] --upconstant [5' constant-region] --downconstant [3' constant region] --minlength 30 --minqual 30"

parser = OptionParser()
parser.add_option('-i', '--input', type = 'string', dest = 'input', help="Path to input fastq file.")
parser.add_option('-o', '--output', type = 'string', dest = 'outfile', help="Path to output fastq file.")
parser.add_option('--upconstant', type = 'string', dest = 'upconstant', help="Upstream (i.e. 5') constant region.")
parser.add_option('--downconstant', type = 'string', dest = 'downconstant', help="Downstream (i.e. 3') constant region.")
parser.add_option("-q", "--minqual", type = "int", dest = "qualityCutoff", help="Specify the minimum average phredscore of valid barcodes.")
parser.add_option("-r", "--regex", type = "string", dest = "regexPattern", help="Specify the regex pattern to search for in the reads.")
parser.add_option("-l", "--minlength", type = "int", dest = "lengthCutoff", help="Specify the minimum length of valid barcodes.")

(options,args) = parser.parse_args()

#------------------
# Parse inputs
#------------------

print("")
print("extractBarcodeReads.py")
print("")

# parse input fastq
if options.input is not None:
    try:
        os.path.isfile(options.input)
        path = options.input
        print("Input file: ", os.path.basename(path))
    except:
        print("no input fastq defined. exiting")
        print(usage)
        sys.exit(1)

# define output 
if options.outfile is not None:
    print("Writing output to ", options.outfile)
    outfile_name = options.outfile
else:
    print("no output file defined. exiting")
    print(usage)
    sys.exit(1)

# parse constant regions
# upconstant
if options.upconstant is not None:
    upstream_constant = options.upconstant
    print("5' constant region: ", upstream_constant.upper())
else:
    upstream_constant = "CGGATCctgaccatgtacgattgacta"
    print("no 5' constant region given, setting constant region to: ", upstream_constant.upper())

# downconstant
if options.downconstant is not None:
    downstream_constant = options.downconstant
    print("3' constant region: ", downstream_constant.upper())
else:
    downstream_constant = "CGGATCctgaccatgtacgattgacta"
    print("no 3' constant region given, setting constant region to: ", downstream_constant.upper())

# parse regex pattern
if options.regexPattern is not None:
    regex = options.regexPattern
    print("Regex set to: ", str(regex))
    barcode_re = re.compile(regex)

else:
    print('no regex pattern given. defaulting to SPLINTR GFP pattern: ([ATCG][ATCG][GC][AT][GC][ATCG][ATCG][AT][GC][AT]){3,6}')
    regex = "([ATCG][ATCG][GC][AT][GC][ATCG][ATCG][AT][GC][AT]){3,6}"
    barcode_re = re.compile(regex)

# parse qualityCutoff
if options.qualityCutoff is not None:
    qualityCutoff = int(options.qualityCutoff)
    print("Quality Cutoff set to: " + str(qualityCutoff))
else: 
    print("No Quality Cutoff given. Defaulting to 30")
    qualityCutoff = 30

# parse lengthCutoff
if options.lengthCutoff is not None:
    lengthCutoff = int(options.lengthCutoff)
    print("Length Cutoff set to: " + str(lengthCutoff))
else: 
    print("No Length Cutoff given. Defaulting to 30")
    lengthCutoff = 30

#------------------
# Filtering function
#------------------

def filter_and_trim_reads(fastq_generator, upstream_constant, downstream_constant, regex, qualityCutoff, lengthCutoff):
        # collect counts of total and filtered reads
            total_count = 0
            filtered_count = 0
            filtered_N = 0
            filtered_qual = 0
            filtered_length = 0
            filtered_upconstant_mismatch = 0
            filtered_downconstant_mismatch = 0
            filtered_pattern_mismatch = 0
        
            # main parsing block
            for record in fastq_generator:
                total_count += 1

            # take only complete reads with correct start and end constant regions
                if str(upstream_constant) in str(record.seq): # TO-DO implement fuzzy matching here. also add in downconstant matching
                    read = str(record.seq)
                    start = read.index(str(upstream_constant)) + len(upstream_constant)
                    end = len(read)
                    trimmed_read = record[start:end]
                    
                    # take as much of the barcode portion as possible i.e. up to 60bp for SPLINTR
                    if(len(trimmed_read.seq)>60):
                        trimmed_read = trimmed_read[0:60]
                        #print(trimmed_read)

                    # check length is OK and return barcode only
                    if len(trimmed_read.seq) >= lengthCutoff:
                        read_quals = trimmed_read.letter_annotations["phred_quality"]
                        avg_qual = sum(trimmed_read.letter_annotations["phred_quality"]) / len(trimmed_read)
                            
                        # barcode must be above average quality threshold
                        if avg_qual >= qualityCutoff:
                                
                            # barcodes cannot contain N residues
                            if "N" not in str(trimmed_read.seq):
                                    
                                # barcode must match pattern -  TO-DO implement fuzzy matching here.
                                if barcode_re.search(str(trimmed_read.seq)):
                                    yield trimmed_read
                                    filtered_count += 1
                                else:
                                    filtered_pattern_mismatch += 1
                            else:
                                filtered_N += 1
                        else:
                            filtered_qual += 1
                    else:
                        filtered_length += 1
                else:
                    filtered_upconstant_mismatch += 1                        

            print("extractBarcodeReads.py")
            print("")
            print("Total reads parsed: ", total_count)
            print("Number of reads passing filters: ", filtered_count)
            print("Number of reads removed - N residues: ", filtered_N)
            print("Number of reads removed - Quality: ", filtered_qual)
            print("Number of reads removed - Length: ", filtered_length)
            print("Number of reads removed - upstream constant mismatch: ", filtered_upconstant_mismatch)
            print("Number of reads removed - downstream constant mismatch: ", filtered_downconstant_mismatch)
            print("Number of reads removed - barcode pattern mismatch: ", filtered_pattern_mismatch)
            
            print("Percentage of reads kept: ", round((filtered_count / total_count) * 100,2))
    
        

#---------------------
# Run Filtering script
#---------------------
def parse_fastq(file, outfile_name, upstream_constant, downstream_constant):

    # get sample name from filename
    filename = os.path.basename(file)
    filename = filename.split("_")[0]
    filename = filename.split(".")[0]
    print("")
    print("Parsing", filename)
    
    # parse constant regions              
    upstream_constant = upstream_constant.upper()
    downstream_constant = downstream_constant.upper()
    
    # setup generator to parse fastq
    with gzip.open(file, 'rt') as fq_handle:
        generator = SeqIO.parse(fq_handle, "fastq")
        
        # parse and write filtered and trimmed reads to file
        with open(outfile_name, "w") as out_handle:
            SeqIO.write(filter_and_trim_reads(generator, upstream_constant, downstream_constant, regex, qualityCutoff, lengthCutoff), out_handle, "fastq")


    # gzip the outfile, can probably do this in one step?
    with open(outfile_name, 'rb') as f_in:
        with gzip.open(outfile_name + '.gz' , 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

parse_fastq(path, outfile_name, upstream_constant, downstream_constant)

# finish and print runtime
print("Script runtime: ", datetime.now() - startTime)
