#!/usr/bin/python

# JMG 5/2017
# Count errors in a SAM based on quality scores.
# Version 2: increased efficiency
# Version 3: option to limit analysis to overlap region
# Version 4: separating stitch matches and mismatches

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

def parseCigar(cigar, diff):
  '''
  Save positions of inserted bases.
  '''
  ops = re.findall(r'(\d+)([IM])', cigar)
  pos = 0
  for op in ops:
    if op[1] == 'I':
      for i in range(pos, pos + int(op[0])):
        diff[i] = 2
    pos += int(op[0])

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

def loadMismatch(filename, d):
  '''
  Load stitching mismatches to given dict.
  '''
  f = openRead(filename)
  for line in f:
    spl = line.rstrip().split('\t')
    if spl[0] in d:
      d[spl[0]].append(spl[1:])
    else:
      d[spl[0]] = [spl[1:]]

  if f != sys.stdin:
    f.close()

def findDiffs(diff, md):
  '''
  Find positions of substitutions using MD.
  '''
  loc = 0  # location on read
  parts = re.findall(r'(\d+|\D+|\^\D+)', md)
  for part in parts:
    try:
      # sequence match
      val = int(part)

      # adjust for insertions (which do not consume MD parts)
      while val:
        if diff[loc] != 2:
          val -= 1
        loc += 1

    except ValueError:

      # skip deletion
      if part[0] == '^':
        pass

      else:

        # substitution
        for i in range(len(part)):
          # skip inserted bases
          while diff[loc] == 2:
            loc += 1
          diff[loc + i] = 1
          loc += 1

def countBases(res, res2, diff, qual, start, end, pos):
  '''
  Count matches/mismatches/insertions/Ns.
  'pos' has list of positions of stitch mismatches --
    save the results to res2.
  '''
#  if len(pos) > 0: print res
  for i in range(start, end):
    q = ord(qual[i]) - 33  # assume Sanger scale
    if q > 40:
      sys.stderr.write('Error! Quality score outside of Sanger range\n')
      sys.exit(-1)
    if i in pos:
      res2[ q ][ diff[i] ] += 1
    else:
      res[ q ][ diff[i] ] += 1
#    if len(pos) > 0:
#      sys.stdout.write(qual[i])
#  if len(pos) > 0:
#    print
#    print res
#    raw_input()

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

def processSAM(fIn, fOut, length, mismatch):
  '''
  Process the SAM file. Count errors.
  '''
  res = [[0, 0, 0, 0] for i in range(41)]  # for collecting results
  res2 = [[0, 0, 0, 0] for i in range(41)]  # for results of stitch mismatches
  count = 0
  for line in fIn:
    if line[0] == '@': continue
    spl = line.rstrip().split('\t')
    if len(spl) < 11:
      sys.stderr.write('Error! Poorly formatted SAM record:\n' + line)
      sys.exit(-1)
    flag = int(spl[1])
    if flag & 0x904: continue  # skip unmapped, sec/supp

    # determine positions of differences using CIGAR and MD
    diff = [0] * len(spl[9])  # assume matches at every position (value=0)
    parseCigar(spl[5], diff)  # add insertions (value=2)
    findDiffs(diff, getTag(spl[11:], 'MD'))  # add substitutions (value=1)
    for i in range(len(spl[9])):  # add Ns (value=3)
      if spl[9][i] == 'N':
        diff[i] = 3

    # determine positions to analyze
    start = 0
    end = len(spl[10])
    if length:
      start = max(len(spl[10]) - length, 0)
      end = min(length, len(spl[10]))

    # find positions of stitch mismatches
    pos = []
    if spl[0] in mismatch:
      for t in mismatch[spl[0]]:
        if flag & 0x10:
          pos.append(len(spl[10])-int(t[0])-1) # adjust if rc
        else:
          pos.append(int(t[0]))

    # count bases by quality score
#    print spl[0]
    countBases(res, res2, diff, spl[10], start, end, pos)
    count += 1

  printOutput(fOut, res)
  if mismatch:
    fOut.write('\nStitch mismatches:\n')
    printOutput(fOut, res2)
  sys.stderr.write('Reads analyzed: %d\n' % count)

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python countErrors4.py  <inSAM>  ' \
      + '<out>  [<length>]  [<alnFile>]\n')
    sys.exit(-1)

  # open SAM file
  fIn = openRead(args[0])
  fOut = openWrite(args[1])
  length = 0   # length of original reads
  if len(args) > 2:
    length = int(args[2])
  mismatch = dict()
  if len(args) > 3:
    loadMismatch(args[3], mismatch)

  # process SAM file
  processSAM(fIn, fOut, length, mismatch)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
