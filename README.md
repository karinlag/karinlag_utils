karinlag_utils
==============

Repository for HTS utility scripts

Contents:

###same_mid_quantitate.py

##Purpose: detect,trim and quantify MID tags in both ends of 454 sequences.

```
Options:
  -h, --help            show this help message and exit
  -f FILE, --fastafile=FILE
                        fasta filename
  -q FILE, --qualfile=FILE
                        qual file, optional
  -m FILE, --midfile=FILE
                        midfile name
  -v INT, --varsize=INT
                        variable region size
  -o STRING, --outnametag=STRING
                        output filename tag

```

This script takes as input a fasta file, and optionally a qual file, together with a file listing MID tag names and sequences. The script goes through each sequence and tries to find any of the listed tags in either end of the sequence (in the 3´ end, the reverse complement is searched for). The tag should be found within the set variable size (-v) in either end, otherwise it is not considered as found. If a tag is found (in either end) within the varsize limit, any trailing characters are removed, and the sequenced is output with the tag name(s) added to the fasta/qual description line. Note: a sequence will be included in the output only if a forward MID is found. Additionally, the sequence should in the end either have the same mid reverse complemented, or no mid should be found. 



