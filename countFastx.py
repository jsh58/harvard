#!/usr/bin/python

# JMG 8/2017

# Count ACGTNX in a fasta/fastq file.

import sys
import gzip

def procSeq(seq, d):
  '''Count nucs in a seq.'''
  for char in seq:
    d[char] = d.get(char, 0) + 1

def main():
  args = sys.argv[1:]
  if len(args) < 1:
    sys.stderr.write('Need fasta/fastq file on CL\n')
    sys.exit(-1)

  if args[0][-3:] == '.gz':
    f = gzip.open(args[0], 'rb')
  else:
    f = open(args[0], 'rU')
  d = dict()

  # determine file format
  line = f.readline()
  if line[0] not in ['>', '@']:
    sys.stderr.write('Error! File not in fasta/fastq format\n')
    sys.exit(-1)
  char = line[0]

  # fasta file
  count = 0
  if char == '>':
    seq = ''
    while line:
      if line[0] == '>':
        if seq:
          procSeq(seq, d)
          count += 1
        seq = ''
      else:
        seq += line.rstrip()
      line = f.readline()
    if seq:
      procSeq(seq, d)
      count += 1

  # fastq file
  else:
    while line:
      if line[0] != '@':
        sys.stderr.write('Error! File not in fastq format\n')
        sys.exit(-1)
      for i in range(3):
        line = f.readline()
        if i == 0:
          procSeq(line.rstrip(), d)
          count += 1
      line = f.readline()

  print 'Reads processed:', count
  for nuc in d:
    print nuc + '\t' + str(d[nuc])

if __name__ == '__main__':
  main()
