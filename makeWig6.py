#!/usr/bin/python

# JMG 2/16/18

# transform bedgraph -- average values within given window size

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

def loadChromSizes(f):
  '''Load chromosome sizes.'''
  chrom = dict()
  for line in f:
    spl = line.rstrip().split('\t')
    if len(spl) < 2:
      sys.stderr.write('Error! "chrom.size" file improperly formatted\n')
      sys.exit(-1)
    chrom[spl[0]] = int(spl[1])
  if f != sys.stdin:
    f.close()
  return chrom

def parseFiles(fIn, fOut, d, window):
  '''
  Create windows from input bedGraph.
  Average values within windows.
  '''
  chrom = ''
  start = end = pos = 0
  value = 0.0

  for line in fIn:
    spl = line.rstrip().split('\t')
    if len(spl) < 4:
      sys.stderr.write('Error! Input bedGraph improperly formatted\n')
      sys.exit(-1)

    if spl[0] != chrom:
      # process last interval(s) from previous chrom
      while end and start < d[chrom]:
        if not pos:
          value = 'NA'
        else:
          value /= pos
        fOut.write('\t'.join(map(str, [chrom, start, end, value])) + '\n')
        start = end
        end = min(start + window, d[chrom])
        pos = 0
      if chrom in d:
        del d[chrom]

      # start new chromosome
      chrom = spl[0]
      pos = 0
      value = 0.0
      if chrom not in d:
        sys.stderr.write('Error! Chromosome %s not in chrom.size file\n' % chrom)
        sys.stderr.write('  (or, input bedGraph not sorted)\n')
        sys.exit(-1)

      # print empty intervals at beginning
      start = 0
      end = min(start + window, d[chrom])
      while start < d[chrom] and end < int(spl[1]):
        fOut.write('\t'.join(map(str, [chrom, start, end, 'NA'])) + '\n')
        start = end
        end = min(start + window, d[chrom])

    # load values from bed interval
    for i in range(int(spl[1]), int(spl[2])):

      # process interval
      while i >= end:
        if not pos:
          value = 'NA'
        else:
          value /= pos
        fOut.write('\t'.join(map(str, [chrom, start, end, value])) + '\n')

        # start new region
        start = end
        end = min(start + window, d[chrom])
        pos = 0
        value = 0.0

      pos += 1
      value += float(spl[3])

  # process last interval
  while end and start < d[chrom]:
    if not pos:
      value = 'NA'
    else:
      value /= pos
    fOut.write('\t'.join(map(str, [chrom, start, end, value])) + '\n')
    start = end
    end = min(start + window, d[chrom])
    pos = 0

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 4:
    sys.stderr.write('Usage: python %s  <BGin>  ' % sys.argv[0] \
      + '<chrom.size>  <windowSize>  <BGout>\n' \
      + '  <chrom.size>  File listing chromosome names and lengths\n' \
      + '                  (one per line, tab-separated)\n')
    sys.exit(-1)

  # load chromosome sizes
  fChr = openRead(args[1])
  chrom = loadChromSizes(fChr)

  window = int(args[2])
  if window <= 0:
    sys.stderr.write('Error! Window size must be > 0\n')
    sys.exit(-1)
  fIn = openRead(args[0])
  fOut = openWrite(args[3])

  # parse input, produce output
  parseFiles(fIn, fOut, chrom, window)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
