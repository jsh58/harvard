#!/usr/bin/python

# JMG 12/2017

# Filter a fasta file (e.g. the NCBI nt database):
#   - remove sequences of pure Ns
#   - remove sequences shorter than a specified
#       length (def. 25bp)
#   - remove sequences whose headers are in a
#       given list

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

def parseFasta(fIn, fOut, minLen, headers):
  '''
  Parse fasta file, write output on the fly.
  '''
  count = short = pureNs = xReads = total = 0
  head = ''    # header (1st space-delim token)
  read = ''    # full read (header + sequence)
  nseq = True  # sequence is pure Ns
  length = 0   # length of sequence

  # analyze fasta reads
  for line in fIn:
    if line[0] == '>':

      # process previous read
      if read:
        count += 1
        if length < minLen:
          short += 1
        elif nseq:
          pureNs += 1
        elif head in headers:
          xReads += 1
        else:
          fOut.write(read)
          total += 1

      # start new read
      read = line
      head = line.split(' ')[0][1:]
      length = 0
      nseq = True

    elif read:
      # save sequence
      read += line
      length += len(line) - 1
      if line.rstrip() != 'N' * (len(line) - 1):
        nseq = False

  if fIn != sys.stdin:
    fIn.close()

  # process last read
  if read:
    count += 1
    if length < minLen:
      short += 1
    elif nseq:
      pureNs += 1
    elif head in headers:
      xReads += 1
    else:
      fOut.write(read)
      total += 1

  if fOut != sys.stdout:
    fOut.close()

  return count, short, pureNs, xReads, total

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python filterNT.py  <input>  ' \
      + '<output>  [<minLen>]  [<headers]\n')
    sys.stderr.write('  <minLen>    Minimum sequence length (def. 25bp)\n')
    sys.stderr.write('  <headers>   File listing headers of sequences to exclude\n')
    sys.exit(-1)

  # get CL args
  fIn = openRead(args[0])
  fOut = openWrite(args[1])
  minLen = 25
  if len(args) > 2:
    minLen = int(args[2])

  # load headers of seqs to exclude
  headers = []
  if len(args) > 3:
    fRead = openRead(args[3])
    for line in fRead:
      headers.append(line.rstrip())

  # parse fasta
  count, short, pureNs, xReads, total \
    = parseFasta(fIn, fOut, minLen, headers)

  sys.stderr.write('Total fasta sequences in %s: %d\n' % (args[0], count))
  sys.stderr.write('  Shorter than %dbp: %d\n' % (minLen, short))
  sys.stderr.write('  Pure Ns: %d\n' % pureNs)
  sys.stderr.write('  Excluded: %d\n' % xReads)
  sys.stderr.write('  Written to %s: %d\n' % (args[1], total))

if __name__ == '__main__':
  main()
