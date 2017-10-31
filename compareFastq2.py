#!/usr/bin/python

# JMG 8/2017

# Reconstructing alignments not matching stitch,
#   checking for mismatches.

import sys
import gzip

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

def revComp(dna):
  '''
  Reverse-complements the given DNA sequence.
  '''
  rc = ''
  for nuc in dna[::-1]:
    comp = ''
    if nuc == 'A': comp = 'T'
    elif nuc == 'C': comp = 'G'
    elif nuc == 'G': comp = 'C'
    elif nuc == 'T': comp = 'A'
    elif nuc == 'N': comp = 'N'
    else:
      print 'Error! Unknown nucleotide: %s' % nuc
    rc += comp
  return rc

def loadFastq(f1):
  '''Return dict of read lengths in a fastq file.'''
  d = dict()
  count = 0
  line = f1.readline()
  while line:
    head = line.split(' ')[0]
    if head[0] != '@':
      sys.stderr.write('Error! Not FASTQ format\n')
      sys.exit(-1)
    for i in xrange(3):
      line = f1.readline()
      if i == 0:
        seq = line.rstrip()
      #elif i == 2:
      #  qual = line
    d[head] = len(seq)
    count += 1
    line = f1.readline()
  f1.close()
  return d, count

def compareFastq(f2, d):
  '''Compare read lengths. Return dict of unmatched reads.'''
  d2 = dict()
  count = miss = diff = 0
  line = f2.readline()
  while line:
    head = line.split(' ')[0]
    if head[0] != '@':
      sys.stderr.write('Error! Not FASTQ format\n')
      sys.exit(-1)
    for i in xrange(3):
      line = f2.readline()
      if i == 0:
        seq = line.rstrip()
      elif i == 2:
        qual = line.rstrip()
    if head not in d:
      miss += 1
      d2[head] = (seq, qual)
      #print head, 'not in d'
    elif d[head] != len(seq):
      diff += 1
      d2[head] = (seq, qual)
      #print d[head], head, seq
      #raw_input()

    count += 1
    line = f2.readline()
  f2.close()
  return d2, count, miss, diff

def retrieveReads(r1, r2, d):
  '''Retrieve reads from r1/r2 based on headers in d.'''
  raw = dict()
  count = len(d)
  line1 = r1.readline().split(' ')[0]
  line2 = r2.readline().split(' ')[0]
  while line1 and line2 and count:
    match = False
    head = line1
    if head[0] != '@':
      sys.stderr.write('Error! Not FASTQ format\n')
      sys.exit(-1)
    if line1 != line2:
      sys.stderr.write('Error! R1/R2 files do not match:\n' + line1 + line2)
      sys.exit(-1)
    if head in d:
      match = True
    for i in xrange(3):
      line1 = r1.readline()
      line2 = r2.readline()
      if match:
        if i == 0:
          seq1 = line1.rstrip()
          seq2 = line2.rstrip()
        elif i == 2:
          qual1 = line1.rstrip()
          qual2 = line2.rstrip()
    if match:
      raw[head] = [seq1, qual1, revComp(seq2), qual2[::-1]]
      count -= 1
    line1 = r1.readline().split(' ')[0]
    line2 = r2.readline().split(' ')[0]

  return raw

def findDiffs(d2, raw):
  '''Return dict of alignment diffs and Ns.'''
  d = dict()
  for head in raw:
    read = head[1:].rstrip()
    d[read] = list()

    # position alignment according to length in d2 dict
    offset = len(d2[head][0]) - len(raw[head][2])
    if offset < 0:
      raw[head][2] = raw[head][2][-offset:]
      raw[head][3] = raw[head][3][-offset:]
    if offset > 0:
      raw[head][2] = ' ' * (offset) + raw[head][2]
      raw[head][3] = ' ' * (offset) + raw[head][3]

    # check for sequence diffs, Ns
    for i in range(min(len(raw[head][0]), len(raw[head][2]))):
      if (raw[head][0][i] != raw[head][2][i] and raw[head][2][i] != ' ') \
          or raw[head][0][i] == 'N' or raw[head][2][i] == 'N':
        d[read].append('\t'.join([str(i), raw[head][0][i], raw[head][1][i], \
          raw[head][2][i], raw[head][3][i]]))
  return d

def printOutput(f, fOut, d):
  '''Produce output. Copy stitch diffs from f, except
     for reads in d. Append diffs from d.'''
  for line in f:
    spl = line.split('\t')
    if spl[0] not in d:
      fOut.write(line)
  for r in d:
    for i in range(len(d[r])):
      fOut.write(r + '\t' + d[r][i] + '\n')
  f.close()

def main():
  args = sys.argv[1:]
  if len(args) < 6:
    sys.stderr.write('Usage: python compareFastq2.py <ref.fastq> <query.fastq> \ \n' \
      + '  <raw_R1.fastq> <raw_R2.fastq> <stitchDiff.txt> <output>\n')
    sys.exit(-1)

  # load read lengths
  f1 = openRead(args[0])
  d, count = loadFastq(f1)
  sys.stderr.write('Reads in %s: %d\n' % (args[0], count))

  # find reads that do not match
  f2 = openRead(args[1])
  d2, count, miss, diff = compareFastq(f2, d)
  sys.stderr.write('Reads in %s: %d\n' % (args[1], count))
  sys.stderr.write('  Missing  : %d\n' % (miss))
  sys.stderr.write('  DiffLen  : %d\n' % (diff))

  # retrieve original reads
  r1 = openRead(args[2])
  r2 = openRead(args[3])
  raw = retrieveReads(r1, r2, d2)
  sys.stderr.write('Reads retrieved: %d\n' % (len(raw)))

  # determine differences
  diff = findDiffs(d2, raw)

  # print output
  f3 = openRead(args[4])
  fOut = openWrite(args[5])
  printOutput(f3, fOut, diff)

  fOut.close()

if __name__ == '__main__':
  main()
