#!/usr/bin/python

# JMG 6/2017

# Merge and split two FASTQ files:
#   - keep only reads found in both files
#   - write them out to separate outputs

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

def main():
  args = sys.argv[1:]
  if len(args) < 4:
    sys.stderr.write('Usage: python mergeSplit.py  ' \
      + '<inFQ1> <inFQ2>  <outFQ1> <outFQ2>\n')
    sys.exit(-1)

  # open SAM file
  fIn1 = openRead(args[0])
  fIn2 = openRead(args[1])
  fOut1 = openWrite(args[2])
  fOut2 = openWrite(args[3])

  # load reads from 1st input file
  d = dict()
  flag = True
  while flag:
    head = seq = qual = ''
    for i in range(4):
      line = fIn1.readline().rstrip()
      if not line:
        flag = False
        break
      if i == 0:
        head = line
      elif i == 1:
        seq = line
      elif i == 3:
        qual = line
    if flag:
      key = head.split(' ')[0][1:]
      d[key] = (head, seq, '+', qual)
  if fIn1 != sys.stdin:
    fIn1.close()

  # check reads in 2nd input, write out both
  #   if there's a match
  flag = True
  while flag:
    head = seq = qual = ''
    for i in range(4):
      line = fIn2.readline().rstrip()
      if not line:
        flag = False
        break
      if i == 0:
        head = line
      elif i == 1:
        seq = line
      elif i == 3:
        qual = line
    if flag:
      key = head.split(' ')[0][1:]
      if key in d:
        fOut1.write('\n'.join(d[key]) + '\n')
        fOut2.write('\n'.join([head, seq, '+', qual]) + '\n')
        #sys.exit(0)

if __name__ == '__main__':
  main()
