################################################################################
#
# Copyright (c) 2012, Karin Lagesen, karin.lagesen@bio.uio.no
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
################################################################################

import os, sys, re
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
from optparse import OptionParser

def read_midfile(midfile):
    # Here we read in the mid tagfile. The expected format is 
    # mid-name in col1, mid in col2, the rest is not used
    fh = open(midfile, "r")
    mids = {}
    for line in fh:
        midid = line.split()[0]
        mid = line.split()[1]
        mids[mid] =  midid
    return mids

def doRecord(record, middict, search_for, maxsize, varsize):
    #These two are for reporting purposes, and will be returned with the record
    reportFMid = "NoF"
    reportRMid = "NoR"
    
    # Iff we do not have a forward mid, we discard sequence and just return None
    editedRecord = None

    # Here we examine each record to find the forward and backward mids
    # The way we do this is to have a regexp object (search_for) that
    # is composed of all of the mids that we have taken in. We then
    # search for these in both ends of the record.
    match_objForward = search_for.search(str(record.seq[:maxsize]))
    match_objReverse = search_for.search(str(record.reverse_complement().seq[:maxsize]))
        
    if match_objForward is None:
        return reportFMid, reportRMid, editedRecord

    # Here, we have at least found a forward mid, and we get what it is, and
    # where it is.
    whichmidForward = match_objForward.group()
    forwardMidSize = len(whichmidForward)
    indexForward = match_objForward.start()
        
    reportFMid = middict[whichmidForward]
    
    if match_objReverse is None:
        # If there is no reverse, I set it to the forward one since
        # there is no logical difference between missing and the same
        whichmidReverse=whichmidForward
        reverseMidSize = 0
        indexReverse=0

    else:
        # Here, we have a reverse mid, and we use those
        whichmidReverse = match_objReverse.group()
        reverseMidSize = len(whichmidReverse)
        indexReverse = match_objReverse.start()

        reportRMid = middict[whichmidReverse]

    if indexForward <= varsize and indexReverse <= varsize:
        if whichmidForward == whichmidReverse:
            # Here we test first, if the mids are the same, second, that they
            # are where we expect them to be.
            firstcut = indexForward
            lastcut = len(record.seq) - indexReverse
            editedRecord = record[firstcut:lastcut]
            editedRecord.id = editedRecord.id + "_" + reportFMid + "_" + reportRMid
    else:
        if indexForward > varsize:
            reportFMid = "LongF"
        if indesReverse > varsize:
            reportRMid = "LongR"
            
    return reportFMid, reportRMid, editedRecord
        
def split_on_mid(middict, fastafile, qualfile, varsize, otag):
    # Here we go through all of the records.
    
    # First, I create a regexp object consisting of all of the mids
    # I can then in doRecord search for them. I also figure out
    # the max size of the region that I accept it to be in, by
    # adding the length of the longest mid and the variable size
    # that the user supplied.
    mids = middict.keys()
    search_for = re.compile("|".join(mids))
    maxsize = max([len(x) for x in mids]) + varsize

    output_handleF = open(otag + ".fsa", "w")

    # Here we go through the file(s). It is done like this to allow
    # for not having a qual file.
    totalseqs = 0
    noseqs = 0
    reportDict = {}
    if qualfile != "None":
        output_handleQ = open(otag + ".qual", "w")
        for record in PairedFastaQualIterator(open(fastafile), open(qualfile)):
            totalseqs += 1
            reportFMid, reportRMid, editedRecord = doRecord(record, middict, search_for, maxsize, varsize)
            reportDict[(reportFMid,reportRMid)] = reportDict.get((reportFMid, reportRMid), 0) + 1
            if editedRecord == None:
                continue
            SeqIO.write(editedRecord, output_handleF, "fasta")
            SeqIO.write(editedRecord, output_handleQ, "qual")
            noseqs += 1

        output_handleQ.close()

    else:
        for record in SeqIO.parse(open(fastafile, "rU"), "fasta"):
            totalseqs += 1
            reportFMid, reportRMid, editedRecord = doRecord(record, middict, search_for, maxsize, varsize)
            reportDict[(reportFMid,reportRMid)] = reportDict.get((reportFMid, reportRMid), 0) + 1
            if editedRecord == None:
                continue        
            SeqIO.write(editedRecord, output_handleF, "fasta")
            noseqs += 1
    
    output_handleF.close()
        
    print "A total of %s sequences written, out of a start set of %s" %(noseqs, totalseqs)
    reportFile = open(otag + ".log", "w")
    reportFile.write("Forward\tReverse\tNumber\n" )
    for key in reportDict:
        category = "Ma"
        if key[0] == "NoF":
            category = "NF"
        elif key[1] == "NoR":
            category = "NR"
        elif key[0] != key[1]:
            category = "MM"
        elif "Long" in key[0] or "Long" in key[1]:
            category = "LM"
    
        reportFile.write("%s\t%s\t%s\t%s\n" % (category, key[0],key[1], reportDict[key]))
    reportFile.close()
    
if __name__ == '__main__':

    parser = OptionParser()
    
    parser.add_option("-f", "--fastafile", dest = "fastafile",
                      help = "fasta filename", metavar = "FILE")
    parser.add_option("-q", "--qualfile", dest = "qualfile", default = "None",
                      help = "qual file, optional", metavar = "FILE")
    parser.add_option("-m", "--midfile", dest = "midfile",
                      help = "midfile name", metavar = "FILE")
    parser.add_option("-v", "--varsize", dest="varsize", type="int",
                      help = "variable region size", metavar = "INT")
    parser.add_option("-o", "--outnametag", dest="otag",
                      help = "output filename tag", metavar="STRING")
    
    (options, args) = parser.parse_args()
    
    mids = read_midfile(options.midfile)
    split_on_mid(mids, options.fastafile, options.qualfile, \
                 options.varsize, options.otag)