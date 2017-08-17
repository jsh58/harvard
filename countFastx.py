#!/usr/bin/python

# JMG 8/2017

# Count ACGTNX in a fasta/fastq file.

import sys

def procSeq(seq, d):
  '''Count nucs in a seq.'''
  for char in seq:
    d[char] = d.get(char, 0) + 1

def main():
  args = sys.argv[1:]
  if len(args) < 1:
    sys.stderr.write('Need fasta/fastq file on CL\n')
    sys.exit(-1)
  f = open(args[0], 'rU')
  d = dict()

  # determine file format
  line = f.readline()
  if line[0] not in ['>', '@']:
    sys.stderr.write('Error! File not in fasta/fastq format\n')
    sys.exit(-1)
  char = line[0]

  # fastq file
  if char == '>':
    seq = ''
    while line:
      if line[0] == '>':
        if seq:
          procSeq(seq, d)
        seq = ''
      else:
        seq += line.rstrip()
      line = f.readline()
    if seq:
      procSeq(seq, d)

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
      line = f.readline()

  for nuc in d:
    print nuc + '\t' + str(d[nuc])

if __name__ == '__main__':
  main()
