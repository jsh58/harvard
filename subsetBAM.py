#!/usr/bin/python

# JMG 12/2018

# output a FASTQ of reads in a SAM whose barcodes match
#   those in a given file

import sys

def revComp(dna):
  '''
  Reverse-complement the given DNA sequence.
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

def getTag(tag, lis):
  '''
  Get optional tag from a SAM record.
  '''
  for t in lis:
    spl = t.split(':')
    if spl[0] == tag:
      return spl[-1]
  return None
  #sys.stderr.write('Error! Cannot find %s in SAM record\n' % tag)
  #sys.exit(-1)

def loadBarcodes(filename):
  '''Load barcodes from a given file.'''
  d = dict()
  f = open(filename, 'rU')
  for line in f:
    spl = line.rstrip().split('_')
    if len(spl) > 1:
      d[spl[1]] = 1
  f.close()
  return d

def main():
  '''Main.'''
  # check CL args
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: samtools view <BAM> | ')
    sys.stderr.write('python %s <barcode.csv> <FASTQ>\n' % sys.argv[0])
    sys.exit(-1)

  # load barcodes
  d = loadBarcodes(args[0])
  sys.stderr.write('Barcodes loaded: ' + str(len(d)) + '\n')

  # parse SAM, write output
  count = 0
  p = dict()  # for reads already printed
  f = sys.stdin
  fOut = open(args[1], 'w')
  for line in f:
    if line[0] == '@': continue
    spl = line.rstrip().split('\t')

    # retrieve barcode
    barcode = getTag('CB', spl[11:])
    if not barcode:
      continue
    barcode = barcode.split('-')[0]

    # check if barcode is in dict
    if barcode in d:
      # check if read has been printed already
      if spl[0] in p:
        continue
      p[spl[0]] = 1

      # print fastq record
      fOut.write('@' + spl[0] + '\n')
      if int(spl[1]) & 0x10:
        fOut.write(revComp(spl[9]) + '\n+\n' + spl[10][::-1] + '\n')
      else:
        fOut.write(spl[9] + '\n+\n' + spl[10] + '\n')
      count += 1

  fOut.close()
  sys.stderr.write('Reads written: ' + str(count) + '\n')

if __name__ == '__main__':
  main()
