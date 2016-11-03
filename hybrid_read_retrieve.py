import string
import datetime

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# hybrid_read_retrieve.py
# 
# Use to obtain hybrid reads from fastq files. That is, find reads that 
# are derived from vector:insert junctions. As written, this program only 
# retrieves reads obtained from the 5' hybrid junction (vector:5' end of 
# insert). 
# 
# Provide the program with the name of the input file and the sequence of 
# the vector you require the reads to contain (called "upvector"). Only 
# reads that contain the entire upvector sequence will be considered 
# hybrid reads.
#
# The output file will contain a _hybrids.fastq extension. The output 
# contains only the insert portion of the reads, but with all quality 
# information and header maintained. In cases where a particular read
# contains the reverse complement of the upvector sequence, the output 
# file reports the reverse complement of the original read. This method
# allows the hybrid reads to be mapped in a strand-specific manner to
# the reference genome. 
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# Notes about use:
#
# At present, the user must enter the name of the input file, the upvector
# sequence, and the cutoff value manually each time. In the future, this
# could be streamlined by retrieving those from the command line or an 
# input file.


# Notes about upvector sequence derived from Gateway donor vector pDONR222 
# compared to the SK1 genome (version SK1-20151029-withDonorVec):
#
# Sequence CAAAAAAGTTGG from the junction does not uniquely specify 
# hybrid reads (it appears twice in the endogenous genome). 
# To be cautious, try using GTACAAAAAAGTTGG, which adds three nucleotides 
# to ensure specificity.


# User-Entered Inputs
inputfilename = "EUK001A_S1_L001_R1_001.fastq" # Filename MUST contain only one . character!
upvector = "GTACAAAAAAGTTGG" # confirmed entry vector 15-mer
cutoff = 15 # number of nucleotides of the insert that must be present 
#                to be included in output fastq file

upvectorlength = len(upvector)

def revcomp(seq):
    # returns reverse complement of a sequence containing g/G, c/C, a/A, t/T,
    # and n/N. Caution: No error reporting in cases of non-DNA sequence input!
    seq = string.replace(seq, "g", "Z")
    seq = string.replace(seq, "c", "g")
    seq = string.replace(seq, "Z", "c")
    seq = string.replace(seq, "a", "Z")
    seq = string.replace(seq, "t", "a")
    seq = string.replace(seq, "Z", "t")
    
    seq = string.replace(seq, "G", "Z")
    seq = string.replace(seq, "C", "G")
    seq = string.replace(seq, "Z", "C")
    seq = string.replace(seq, "A", "Z")
    seq = string.replace(seq, "T", "A")
    seq = string.replace(seq, "Z", "T")
    
    return seq[::-1]

def gethybrid(seq, upvector, upvectorlength):
    # input: sequence of one read, upvector sequence, upvector length
    # returns: (1) the insert part of the hybrid read if it's a hybrid read; 
    #           otherwise return False, (2) position in the string where the 
    #           insert starts if a hybrid read; else returns False, (3) the
    #           orientation (forward or reverse) as a string
    # variables:
    # upvector = sequence 5' to the insert that must be matched
    
    # looks for upvector in the forward direction of the read (called seq)
    vectorlocation = string.find(seq, upvector)
    orientation = "none"
    if vectorlocation != -1: # if forward sequence contains vector
        insertstart = vectorlocation + upvectorlength
        hybrid = seq[insertstart:] # sequence of the insert part of the hybrid
        orientation = "forward"
        return hybrid, insertstart, orientation
    else: # else, try looking in the reverse complement, then give up
        revseq = revcomp(seq)
        revvectorlocation = string.find(revseq, upvector)
        if revvectorlocation != -1: # if reverse sequence contains vector
            insertstart = revvectorlocation + upvectorlength
            hybrid = revseq[insertstart:]
            orientation = "reverse"
            return hybrid, insertstart, orientation # note that insertstart is 
        #                          calculated from the reverse orientation
        else:
            return False, False, False

def run(inputfilename, upvector, upvectorlength, cutoff):
    
    # input: name of file
    # returns: none (makes a new fastq file with only hybrid reads in 
    #   same directory)
    
    # reads the input fastq file
    f = open(inputfilename, "r")
    contents = f.read()
    f.close()
    
    # creates output file name and empty output string to append to
    basefilename = string.split(inputfilename, ".")[0]
    outputfilename = basefilename + "_hybrids.fastq"
    output = ""
    
    # creates report file
    reportfilename = basefilename + "_hybrids_report.txt"
    
    # creates a list of fastq reads by splitting at the @ sign
    reads = string.split(contents, "@")
    reads.pop(0) # removes the first entry, which is an empty string
    
    # for each read, tests if a hybrid read. if hybrid, then adds to output
    # calculates statistics for hybrid reads to report in a report text file
    numlonghybrid = 0 # hybrids meeting the cutoff
    numshorthybrid = 0 # hybrids not meeting the cutoff
    numnothybrid = 0 # number of reads that are not hybrid
    for r in reads:
        r = string.rstrip(r) # removes trailing \n or spaces
        data = string.split(r, "\n")
        header = data[0] # header line
        read = data[1] # read line
        delimiter = data[2] # delimiter line
        quality = data[3] # Phred quality score line
        hybrid, insertstart, orientation = gethybrid(read, upvector, upvectorlength)
        if hybrid:
            if len(hybrid) >= cutoff: # if hybrid and meets cutoff
                # Modifies Quality Scores to match hybrid read
                if orientation == "forward":
                    quality = quality[insertstart:]
                elif orientation == "reverse":
                    quality = quality[::-1]
                    quality = quality[insertstart:]                
                # appends to output and updates statistics
                output += "\n".join(["\n@"+header, hybrid, delimiter, quality])
                numlonghybrid += 1
            else: # else hybrid does not meet length cutoff
                numshorthybrid += 1
        else: # else not a hybrid
            numnothybrid += 1
    
    # writes to output file
    output = string.strip(output)
    outputfile = open(outputfilename, "w")
    outputfile.write(output)
    outputfile.close()
    
    # makes timestamp
    now = datetime.datetime.now()
    timestamp = now.isoformat()
    
    # creates report file with input file name and # of hybrid/nonhybrid reads
    reportfile = open(reportfilename, "w")
    reportline1 = "Report for sequence file:\t%s" % inputfilename
    reportline2 = "Operation completed:\t%s" % timestamp
    reportline4 = "total hybrid reads:\t%s" % (numlonghybrid + numshorthybrid)
    reportline5 = "hybrid reads meeting cutoff:\t%s" % numlonghybrid
    reportline6 = "hybrid reads below cutoff:\t%s" % numshorthybrid
    reportline7 = "non-hybrid reads:\t%s" % numnothybrid
    reportline3 = "total reads:\t%s" % (numlonghybrid + numshorthybrid + \
                                        numnothybrid)
    reportline8 = "insert length cutoff:\t%s" % cutoff
    reportline9 = "upvector sequence:\t%s" % upvector
    report = "\n".join([reportline1, reportline2, reportline3, reportline4,\
                        reportline5, reportline6, reportline7, reportline8, \
                        reportline9])
    reportfile.write(report)
    reportfile.close()
    
run(inputfilename, upvector, upvectorlength, cutoff)
print "DONE processing file %s with upvector %s and cutoff %d nt" % \
      (inputfilename, upvector, cutoff)