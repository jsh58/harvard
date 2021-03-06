#!/usr/bin/python

# JMG 5/2017
# Count errors in a SAM based on quality scores.

import sys
import gzip
import re

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

def parseCigar(cigar):
  '''
  Return string representation of CIGAR.
  '''
  ops = re.findall(r'(\d+)([IDM])', cigar)
  cigar = ''
  for op in ops:
    cigar += int(op[0]) * op[1]
  return cigar

def getTag(lis, tag):
  '''
  Get optional tag from a SAM record.
  '''
  for t in lis:
    spl = t.split(':')
    if spl[0] == tag:
      return spl[-1]
  sys.stderr.write('Error! Cannot find %s in SAM record\n' % tag)
  sys.exit(-1)

def findDiffs(cigar, md):
  '''
  Find positions of differences using CIGAR and MD.
  '''
  loc = 0  # location on cigar
  parts = re.findall(r'(\d+|\D+|\^\D+)', md)
  for part in parts:
    try:
      # sequence match
      val = int(part)

      # adjust for insertions (which do not consume MD parts)
      while val:
        if cigar[loc] != 'I':
          val -= 1
        loc += 1

    except ValueError:

      if part[0] == '^':
        # deletion: remove from cigar
        i = 0
        while cigar[loc + i] == 'D':
          i += 1
        cigar = cigar[:loc] + cigar[loc + i:]
        if not i:
          sys.stderr.write('Error! MD deletion not in CIGAR\n' \
            + 'CIGAR: ' + cigar \
            + '\nMD: ' + md + '\n')
          sys.exit(0)
        if len(part) > i + 1:
          sys.stderr.write('Error! MD deletion not separated ' \
            + 'from substitution: ' + md + '\n')
          sys.exit(0)

      else:

        # substitution -- put 'X' in cigar
        for i in range(len(part)):
          # skip inserted bases
          while cigar[loc] == 'I':
            loc += 1
          cigar = cigar[:loc + i] + 'X' + cigar[loc + i + 1:]
          loc += 1

  return cigar

def countBases(res, diff, seq, qual):
  '''
  Count matches/mismatches/insertions/Ns.
  '''
  idx = {'M': 0, 'X': 1, 'I': 2}  # values for CIGAR characters
  for i in range(len(diff)):
    q = ord(qual[i]) - 33  # assume Sanger scale
    val = idx[diff[i]]
    if seq[i] == 'N':
      val = 3
    res[ q ][ val ] += 1

def printOutput(fOut, res):
  '''
  Produce output.
  '''
  fOut.write('\t'.join(['qual', 'match', 'sub', 'ins', 'N',
    'subRate', 'insRate', 'NRate']) + '\n')
  for i in range(len(res)):
    fOut.write('\t'.join(map(str, [i] + res[i])))
    total = float(sum(res[i]))
    for j in range(1, 4):
      if total:
        fOut.write('\t%.9f' % (res[i][j] / total))
      #else:
      #  fOut.write('\tNA')
    fOut.write('\n')

def processSAM(fIn, fOut):
  '''
  Process the SAM file. Count errors.
  '''
  res = [[0, 0, 0, 0] for i in range(41)]  # for collecting results
  for line in fIn:
    if line[0] == '@': continue
    spl = line.rstrip().split('\t')
    if len(spl) < 11:
      sys.stderr.write('Error! Poorly formatted SAM record:\n' + line)
      sys.exit(-1)
    flag = int(spl[1])
    if flag & 0x904: continue  # skip unmapped, sec/supp

    # determine positions of differences with CIGAR and MD
    cigar = parseCigar(spl[5])
    md = getTag(spl[11:], 'MD')
    diff = findDiffs(cigar, md)

    # count bases
    countBases(res, diff, spl[9], spl[10])

  printOutput(fOut, res)

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python countErrors.py  <inSAM>  <out>\n')
    sys.exit(-1)

  # open SAM file
  fIn = openRead(args[0])
  fOut = openWrite(args[1])

  # process SAM file
  res = processSAM(fIn, fOut)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
