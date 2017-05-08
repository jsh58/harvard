#!/usr/bin/python

# JMG 5/8/17

# make fixedStep wiggle from bedGraph (4-column)

import sys

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
    sys.stderr.write('Error! Cannot read input file %s\n' % filename)
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
    sys.stderr.write('Error! Cannot write to output file %s\n' % filename)
    sys.exit(-1)
  return f

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python %s  <bedGraph>  <wig>\n' % sys.argv[0])
    exit(-1)

  # open files
  fIn = openRead(args[0])
  fOut = openWrite(args[1])
  fOut.write('track type=wiggle_0')
  # print basename to header
  #name = '.'.join(args[0].split('.')[:-1])
  #fOut.write(' name=' + name)
  fOut.write('\n')

  # process files
  chrom = ''
  step = span = 0
  line = fIn.readline()
  while line:

    # load data
    newReg = False
    spl = line.rstrip().split('\t')
    if spl[0] != chrom:
      chrom = spl[0]
      newReg = True
    start = int(spl[1])
    length = int(spl[2]) - start
    dist = length

    # calculate distance to next step
    nextLine = fIn.readline()
    if nextLine and nextLine.split('\t')[0] == chrom:
      dist = int(nextLine.split('\t')[1]) - start
    else:
      dist = step  # avoids writing new header for last interval of chromosomes

    # initialize new region, if necessary
    #if newReg or dist != step or length != span:
    if newReg or dist != step:
      chrom = spl[0]
      step = dist
      span = length
      fOut.write('fixedStep chrom=' + chrom \
        + ' start=' + str(start + 1) \
        + ' step=' + str(step) \
        + ' span=' + str(span) \
        + '\n')

    # print output
    fOut.write(spl[3] + '\n')
    line = nextLine

if __name__ == '__main__':
  main()
