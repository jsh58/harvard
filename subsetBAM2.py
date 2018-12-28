#!/usr/bin/python

# JMG 12/2018

# output a FASTQ of reads in a SAM whose alignments are
#   to hg19 (better than mm10)

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
    sys.stderr.write('python %s hg19 <FASTQ>\n' % sys.argv[0])
    sys.exit(-1)

  # load name of genome of interest
  gen = args[0]

  # parse SAM, save reads whose AS's are best and alns are to 'gen'
  score = dict()  # dict of alignment scores
  reads = dict()  # dict of reads (seq/qual)
  f = sys.stdin
  for line in f:
    if line[0] == '@': continue
    spl = line.rstrip().split('\t')
    flag = int(spl[1])
    if flag & 0x4: continue

    # retrieve alignment score
    aScore = getTag('AS', spl[11:])
    if not aScore:
      sys.stderr.write('missing AS: ' + spl[0] + '\n')
      continue
    aScore = int(aScore)

    # check scores dict
    if spl[0] in score:

      if aScore > score[spl[0]]:
        # better AS
        score[spl[0]] = aScore
        if spl[2][:len(gen)] == gen:
          # save read
          if flag & 0x10:
            reads[spl[0]] = (revComp(spl[9]), spl[10][::-1])
          else:
            reads[spl[0]] = (spl[9], spl[10])
        elif spl[0] in reads:
          # better aln not to 'gen': delete
          del reads[spl[0]]

      elif aScore == score[spl[0]]:
        # equal AS
        if spl[2][:len(gen)] == gen:
          # save read
          if flag & 0x10:
            reads[spl[0]] = (revComp(spl[9]), spl[10][::-1])
          else:
            reads[spl[0]] = (spl[9], spl[10])

    else:
      # not seen before: save AS
      score[spl[0]] = aScore
      if spl[2][:len(gen)] == gen:
        if flag & 0x10:
          reads[spl[0]] = (revComp(spl[9]), spl[10][::-1])
        else:
          reads[spl[0]] = (spl[9], spl[10])

  f.close()

  # print fastq
  count = 0
  fOut = open(args[1], 'w')
  for r in reads:
    fOut.write('@' + r + '\n' + reads[r][0] + '\n+\n' + reads[r][1] + '\n')
    count += 1

  fOut.close()
  sys.stderr.write('Reads written: ' + str(count) + '\n')

if __name__ == '__main__':
  main()
