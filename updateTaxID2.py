#!/usr/bin/python

# JMG 12/2017
# Update merged and deleted taxonomic IDs.

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
    sys.stderr.write('Usage: python updateTaxID2.py  <mergedIDs>  ' \
      + '<deletedIDs>  <out>  [<in>]+\n')
    sys.exit(-1)

  # load merged taxIDs to dict
  d = dict()
  f = openRead(args[0])
  for line in f:
    spl = line.rstrip().split('|')
    if len(spl) < 2:
      sys.stderr.write('Error! Poorly formatted record in merged file\n')
      sys.exit(-1)
    d[spl[0].strip()] = spl[1].strip()
  if f != sys.stdin:
    f.close()

  # load deleted taxIDs to dict
  f2 = openRead(args[1])
  for line in f2:
    spl = line.rstrip().split('|')
    d[spl[0].strip()] = '0'  # assign deleted to tax ID '0'
  if f != sys.stdin:
    f2.close()

  # open output file
  merge = printed = 0
  fOut = openWrite(args[2])

  # parse input files, write output on the fly
  for arg in args[3:]:
    fIn = openRead(arg)

    # parse header
    accIdx = taxIdx = -1
    spl = fIn.readline().rstrip().split('\t')
    try:
      accIdx = spl.index('accession.version')
      taxIdx = spl.index('taxid')
    except ValueError:
      sys.stderr.write('Error! Cannot find header value '
        + '(\'accession.version\' or \'taxid\')')
      sys.exit(-1)

    # parse input file, produce output
    for line in fIn:
      spl = line.rstrip().split('\t')
      if spl[taxIdx] in d:
        spl[taxIdx] = d[spl[taxIdx]]
        merge += 1
      fOut.write(spl[accIdx] + '\t' + spl[taxIdx] + '\n')
      printed += 1
    if fIn != sys.stdin:
      fIn.close()

  fOut.close()
  sys.stderr.write('Records written: %d\n' % printed)
  sys.stderr.write('  Updated: %d\n' % merge)

if __name__ == '__main__':
  main()
