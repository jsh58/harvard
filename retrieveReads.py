#!/usr/bin/python

# JMG 8/2017

# Retrieve a subset of reads.
#   - load headers of mapped reads from a SAM
#   - copy those reads from fastq files to output fastq files

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

def loadReads(f):
  '''
  Load read headers for aligned reads.
  '''
  d = dict()
  #count = 0
  for line in f:
    if line[0] == '@': continue
    spl = line.rstrip().split('\t')
    if len(spl) < 11: continue
    d[spl[0]] = 1
    #count += 1
  sys.stderr.write('Reads loaded: %d\n' % len(d))
  return d

def main():
  args = sys.argv[1:]
  if len(args) < 5:
    sys.stderr.write('Usage: python retrieveReads.py  ' \
      + '<inSAM>  <inFQ1> <inFQ2>  <outFQ1> <outFQ2>\n')
    sys.exit(-1)

  # open SAM file
  fIn = openRead(args[0])
  d = loadReads(fIn)
  if fIn != sys.stdin:
    fIn.close()

  # open fastq files
  fIn1 = openRead(args[1])
  fIn2 = openRead(args[2])
  fOut1 = openWrite(args[3])
  fOut2 = openWrite(args[4])

  # check reads from 1st input file
  total = count = 0
  flag = True
  while flag:
    match = False
    for i in range(4):
      line = fIn1.readline()
      if not line:
        flag = False
        break
      if i == 0:
        if line.split(' ')[0][1:].rstrip() in d:
          match = True
      if match:
        fOut1.write(line)
    if flag:
      total += 1
      if match:
        count += 1
  if fIn1 != sys.stdin:
    fIn1.close()
  if fOut1 != sys.stdout:
    fOut1.close()
  sys.stderr.write('Reads in %s: %d\n' % (args[1], total))
  sys.stderr.write('  printed: %d\n' % count)

  # check reads from 2nd input file
  total = count = 0
  flag = True
  while flag:
    match = False
    for i in range(4):
      line = fIn2.readline()
      if not line:
        flag = False
        break
      if i == 0:
        if line.split(' ')[0][1:].rstrip() in d:
          match = True
      if match:
        fOut2.write(line)
    if flag:
      total += 1
      if match:
        count += 1
  if fIn2 != sys.stdin:
    fIn2.close()
  if fOut2 != sys.stdout:
    fOut2.close()
  sys.stderr.write('Reads in %s: %d\n' % (args[2], total))
  sys.stderr.write('  printed: %d\n' % count)

if __name__ == '__main__':
  main()
