#!/usr/bin/python

# JMG 5/2017
# Count errors in a SAM based on quality scores.

import sys
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

def findDiffs(cigar, md, pr):
  '''
  Find positions of differences using CIGAR and MD.
  '''
  if pr:
    print cigar
  pos = 0  # position on md string
  loc = 0  # location on cigar
  parts = re.findall(r'(\d+|\D+|\^\D+)', md)
  if pr:
    print parts
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
        if len(part) > i + 1:
          sys.stderr.write('Error! MD deletion not separated ' \
            + 'from substitution: ' + md + '\n')
          sys.exit(0)

      else:

        # substitution -- put 'X' in cigar
        #ins = 0
        for i in range(len(part)):
          # skip inserted bases
          while cigar[loc] == 'I':
            loc += 1
          cigar = cigar[:loc + i] + 'X' + cigar[loc + i + 1:]
          #cigar = cigar[:loc + ins + i] + 'X' + cigar[loc + ins + i + 1:]
          loc += 1
        #loc += len(part) + ins

  if pr:
    print cigar, '\n'
#  sys.exit(0)


def processSAM(f):
  '''
  Process the SAM file. Count errors.
  '''
  for line in f:
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

    ###
    # print certain reads' results
    pr = False
    head = ['HSQ-7001360:311:HYFH5BCXX:1:1103:7011:49414', 'HSQ-7001360:311:HYFH5BCXX:1:1101:2558:2213', 'HSQ-7001360:311:HYFH5BCXX:1:1107:13423:77493']
    if spl[0] in head:
      pr = True
    ###

    diff = findDiffs(cigar, md, pr)

    qual = spl[10]

    pos = 0  # position in seq/qual

#    for op in cigar:
#      print op[0], op[1]
#    print md + '\n' + seq + '\n' + qual + '\n'
#    sys.exit(0)

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python countErrors.py  <inSAM>  <out>\n')
    sys.exit(-1)

  # open SAM file
  fIn = openRead(args[0])
  fOut = openWrite(args[1])

  # process SAM file
  processSAM(fIn)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
