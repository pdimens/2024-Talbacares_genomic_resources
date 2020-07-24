
#!/usr/bin/env python2

from collections import defaultdict
from optparse import OptionParser
from sets import Set
import os, sys
import re

def setup_options():
    parser = OptionParser()
    parser.add_option("-b", "--backbone", dest="backbone_filename", \
            help="Backbone fasta file containing contigs with gaps.", metavar="FILE")
    parser.add_option("-c", "--consensus", dest="consensus_filename", \
            help="Fasta file where header corresponds to backbone id and entries correspond to read ids. Typically 'DBG2OLC_Consensus_info.txt' in pipeline.")
    parser.add_option("-r", "--reads", dest="reads_filename", \
            help="Fasta file of the reads. Typically 'ctg_pb.fasta' in pipeline. ")
    parser.add_option("-o", "--output-dir", dest="output_dir", \
            help="Output directory of split fasta files.")
    parser.add_option("-p", "--prefix", dest="prefix", default="backbone", \
            help="Default prefix for split fasta.")


    (options, args) = parser.parse_args()

    return (options,args)


## SeqIO.py stuff
'''
 ParseFasta - class for reading a fasta-like formatted file.
 
 The assumption is that the file contains a set of multi-line records
 separated by single-line headers starting with a specific separator 
 (by default >)
  
 - Taken from ParseFasta.pm - Mihai Pop
 
 @author Christopher Hill
'''
import time,sys
from UserString import MutableString

class ParseFasta:
    ''' Constructor '''
    def __init__(self, file, head=">", linesep=""):
        # Head/Record separator, default is ">"
        self.headSep = head
        
        # Line seperator used when concatenating the lines in the input
        # forming the body of each record.  Useful for .qual files
        self.lineSep = linesep
        
        # String buffer
        self.buffer = None
        
        #Represents the line buffer
        self.file = open(file)
        
        # If the file doesn't start with a header, exit
        self.buffer = self.file.readline()
        if self.buffer == None or not self.buffer.startswith(self.headSep):
            return None
    
    '''Return head separator'''        
    def getHeadSep(self):
        return self.headSep
    
    ''' Return line separator '''
    def getLineSep(self):
        return self.lineSep
    
    '''            
    Reads a record and returns the head and data in a String array.
     array[0] = head
     array[1] = data
    '''
    def getRecord(self):
        # Stores the head entry
        head = ""
        
        # Stores the data
        data = ""
        
        # If the buffer doesn't start with a header
        if not self.buffer or not self.buffer.startswith(self.headSep):
            return None
            
        # Set the header
        head = self.buffer[len(self.headSep):]
        head = head.rstrip()
        self.buffer = self.file.readline();
        
        # Set the data by continously looping through the record until a new record
        # is reached or the end of file is reached
        while self.buffer and not self.buffer.startswith(self.headSep):
            data += (self.buffer.rstrip() + self.lineSep ) # Might have to add trim
            self.buffer = self.file.readline()
            
        # Prepare the record
        results = [head, data]
        
        return results;

    ''' Close Stream '''
    def closeStream(self):
        self.file.close()


# Modified from: http://scipher.wordpress.com/2010/05/06/simple-python-fastq-parser/
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
        
    def __iter__(self):
        return self
    
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
        
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line #%s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file and try again **" % (self._hdSyms[0])
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file and try again **" % (self._hdSyms[1])
        
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
    
    def getNextReadSeq(self):
        """Convenience method: calls self.getNext and returns only the readSeq."""
        try:
            record = self.next()
            return record[1]
        except StopIteration:
            return None
'''
a = time.time()
pf = ParseFasta('G:\\Work\\input.fasta')
tuple = pf.getRecord()
while tuple is not None:
    #print tuple[1]
    tuple = pf.getRecord()
    
b= time.time()
print b-a
'''
## end of SeqIO.py stuff

def split_backbone(options):
    """ 
    Split backbone fasta file into chunks. 
    Returns dictionary of backbone -> id.
    """

    backbone_to_id = {}
    id_counter = 0

    # Write all backbone files to their own fasta file.
    pf = ParseFasta(options.backbone_filename)
    tuple = pf.getRecord() # [seqname, sequence]
    while tuple is not None:
        print tuple[0]

        split_backbone = open(options.output_dir + '/' + options.prefix + '-' + str(id_counter) + '.fasta', 'w')
        split_backbone.write('>' + tuple[0] + '\n' + tuple[1])
        split_backbone.close()

        # backbone_to_id[seqname] = fileprefix-number
        backbone_to_id[tuple[0]] = options.prefix + '-' + str(id_counter)

        id_counter += 1
        tuple = pf.getRecord()

    return backbone_to_id


def build_reads_to_backbone_dict(options):
    """
    Return dictionary of reads to the corresponding backbone.
    """

    """
    Usually in DBG2OLC_Consensus_info.txt
    >Backbone_1
    Contig_59850
    Contig_61226
    Contig_87853
    Contig_42901
    Contig_308247
    Contig_180195
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/167/0_4524
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1566/0_9679
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1888/9474_18333
    """
    reads_to_backbone = defaultdict(list)

    curr_backbone = None
    with open(options.consensus_filename, 'r') as cf:
        for line in cf:
            line = line.strip()
            if line.startswith('>'):
                curr_backbone = line[1:]

            else:
                if curr_backbone:
                    print line + '\t' + curr_backbone
                    ## reads_to_backbone[read_or_contig_name].append(backbone_name)
                    reads_to_backbone[line].append(curr_backbone)

    return reads_to_backbone

def build_backbone_to_reads_dict(options):
    """
    Return dictionary of backbones to the corresponding reads.
    """

    """
    Usually in DBG2OLC_Consensus_info.txt
    >Backbone_1
    Contig_59850
    Contig_61226
    Contig_87853
    Contig_42901
    Contig_308247
    Contig_180195
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/167/0_4524
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1566/0_9679
    m130605_000141_42207_c100515142550000001823076608221372_s1_p0/1888/9474_18333
    """
    backbone_to_reads = defaultdict(list)

    curr_backbone = None
    with open(options.consensus_filename, 'r') as cf:
        for line in cf:
            line = line.strip()
            if line.startswith('>'):
                curr_backbone = line[1:]

            else:
                if curr_backbone:
##                    print line + '\t' + curr_backbone
                    ## backbone_to_reads[backbone_name].append(read_or_contig_name)
                    backbone_to_reads[curr_backbone].append(line)

    return backbone_to_reads


def build_readseqs_dict(options):
    readseqs = {}
    pf = ParseFasta(options.reads_filename) ## typically 'ctg_pb.fasta' in pipeline
    tuple = pf.getRecord()
    while tuple is not None:
        readseqs[tuple[0]] = tuple[1]
        tuple = pf.getRecord()
    return readseqs
        

def ensure_dir(f):
    print f
    d = os.path.dirname(f)
    print d
    if not os.path.exists(d):
        os.makedirs(d)


def main():
    (options, args) = setup_options()

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)

    # Split backbone fasta file into chunks.
    backbone_to_id = split_backbone(options)
    # use: backbone_to_id[seqname] = fileprefix-number

    # Find the reads that correspond to a given backbone.
    reads_to_backbone = build_reads_to_backbone_dict(options)
    ## use: reads_to_backbone[read_or_contig_name].append(backbone_name)
    
    
    ## New approach:
    ## open and close files as you need them
    ## A little slower than original

    pf = ParseFasta(options.reads_filename)
    tuple = pf.getRecord()
    id = None
    while tuple is not None:

        if len(reads_to_backbone[tuple[0]]) > 0:
            for backbone in reads_to_backbone[tuple[0]]:
                id = backbone_to_id[backbone]

                fname = options.output_dir + '/' + str(id) + '.reads.fasta'
                
                if os.path.isfile(fname):
                    with open(fname, 'a') as f:
                        f.write('>' + tuple[0] + '\n' + tuple[1] + '\n')

                else:
                    with open(fname, 'w') as f:
                        f.write('>' + tuple[0] + '\n' + tuple[1] + '\n')

        tuple = pf.getRecord()

    
if __name__ == '__main__':
    main()
