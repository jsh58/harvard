#!/usr/bin/python

# JMG 8/2017
# Extract sequence differences from a
#   SeqPrep alignment file.

import sys
import gzip
import re

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for reading\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wb')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

def printDiff(f, read, seq1, seq2):
  '''
  Print differences between two seqs.
  '''
  # determine length of leading/trailing gaps ('---')
  leadGap = 0
  while seq1[leadGap] == '-' or seq2[leadGap] == '-':
    leadGap += 1
  tailGap = min(len(seq1), len(seq2)) - 1
  while seq1[tailGap] == '-' or seq2[tailGap] == '-':
    tailGap -= 1
  for j in xrange(leadGap, tailGap + 1):
    if seq1[j] != seq2[j] and seq1[j] != ' ' and seq2[j] != ' ':
      f.write('\t'.join([read, str(j - leadGap), seq1[j], '!', seq2[j], '!']) + '\n')

def parseFile(fIn, fOut):
  '''
  Produce list of mismatches from SeqPrep alignment file.
  '''
  line = fIn.readline()
  read = ''
  subj = quer = ''
  count = 0
  while line:
    spl = line.rstrip()
    if spl[0:3] == 'ID:':
      read = spl[3:].lstrip().split(' ')[0]
    elif spl[0:5] == 'SUBJ:':
      subj = spl[6:]
    elif spl[0:6] == 'READ1:':
      subj = spl[7:]
    elif spl[0:5] == 'QUER:':
      quer = spl[6:]
      printDiff(fOut, read, subj, quer)
      read = subj = quer = ''
      count += 1
    elif spl[0:6] == 'READ2:':
      quer = spl[7:]
      printDiff(fOut, read, subj, quer)
      read = subj = quer = ''
      count += 1
    line = fIn.readline()
  sys.stderr.write('Reads analyzed: %d\n' % count)

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python SeqPrepDiff.py <in> <out>\n')
    sys.exit(-1)

  # process file
  fIn = openRead(args[0])
  fOut = openWrite(args[1])
  parseFile(fIn, fOut)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
