#!/usr/bin/python

# JMG 7/2018

# Exclude accessions from a fasta file that are
#   assigned to a given taxon's subtree.

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

def parseFasta(fIn, fOut, taxa, acc2tax):
  '''
  Parse fasta file, write output on the fly.
  '''
  count = short = xReads = total = 0
  head = ''    # header (1st space-delim token)
  read = ''    # full read (header + sequence)

  # analyze fasta reads
  for line in fIn:
    if line[0] == '>':

      # process previous read
      if read:
        count += 1
        if head in acc2tax and acc2tax[head] in taxa:
          xReads += 1
        else:
          fOut.write(read)
          total += 1
        if head not in acc2tax:
          short += 1

      # start new read
      read = line
      head = line.rstrip().split(' ')[0][1:]

    else:
      # save sequence
      read += line

  # process last read
  if read:
    count += 1
    if head in acc2tax and acc2tax[head] in taxa:
      xReads += 1
    else:
      fOut.write(read)
      total += 1
    if head not in acc2tax:
      short += 1

  return count, short, xReads, total

def loadAcc(f):
  '''
  Load accession -> taxID info.
  '''
  d = {}
  for line in f:
    spl = line.rstrip().split('\t')
    if len(spl) < 4:
      sys.stderr.write('Error! Improperly formatted acc2taxid file\n')
      sys.exit(-1)
    d[spl[1]] = spl[2]
  return d

def findTaxa(f, taxon):
  '''
  Load all taxa in given taxon's subtree.
  '''
  # load immediate parent taxa and counts to dict
  temp = {}
  for line in f:
    spl = line.split('|')
    if len(spl) < 3:
      sys.stderr.write('Error! Improperly formatted tree file\n')
      sys.exit(-1)
    temp[spl[0].strip()] = spl[1].strip()

  # save taxa from subtree
  d = {}
  d[taxon] = 1
  for tax in temp:
    parent = temp[tax]
    while parent != 'None':
      if parent in d:
        d[tax] = 1
        break
      parent = temp[parent]

  return d

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 3:
    sys.stderr.write('Usage: python %s  ' % sys.argv[0] \
      + '<taxon>  <taxTree>  <acc2taxid>  <in>  <out>\n')
    sys.exit(-1)

  # find all taxa in given taxon's subtree
  taxon = args[0]
  fTax = openRead(args[1])
  taxa = findTaxa(fTax, taxon)
  if fTax != sys.stdin:
    fTax.close()
  sys.stderr.write('Taxa in subtree of %s: %d\n' % (taxon, len(taxa)))

  # load accession -> taxID info
  fAcc = openRead(args[2])
  acc2tax = loadAcc(fAcc)
  if fAcc != sys.stdin:
    fAcc.close()

  # print output
  fIn = openRead(args[3])
  fOut = openWrite(args[4])
  count, short, xReads, total = parseFasta(fIn, fOut, taxa, acc2tax)
  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()
  sys.stderr.write('Total fasta sequences in %s: %d\n' % (args[3], count))
  sys.stderr.write('  Excluded due to %s subtree: %d\n' % (taxon, xReads))
  sys.stderr.write('  (taxonomy unknown: %d)\n' % short)
  sys.stderr.write('  Written to %s: %d\n' % (args[4], total))

if __name__ == '__main__':
  main()
