#!/usr/bin/python

# JMG 3/2017

# find reads that are better aligned in one file than another

import sys

def getTag(tag, lis):
  for l in lis:
    spl = l.split(':')
    if spl[0] == tag:
      return spl[-1]
  print 'Error! Cannot find %s tag\n' % tag
  sys.exit(-1)

def main():
  # load AS's from CL arg
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python %s  <SAM>  <out>\n' % sys.argv[0])
    sys.stderr.write('  NOTE: the 2nd SAM must be piped in\n')
    sys.exit(-1)

  # load AS's from SAM
  d = dict()
  f = open(args[0], 'rU')
  for line in f:
    if line[0] == '@': continue
    spl = line.rstrip().split('\t')
    aScore = int(getTag('AS', spl[11:]))
    if spl[0] in d:
      if aScore != d[spl[0]]:
        sys.stderr.write('Warning: not a match:' + spl[0] \
          + ' ' + str(d[spl[0]]) + ' ' + str(aScore) + '\n')
    else:
      d[spl[0]] = aScore
  f.close()
  print 'Trimmed/aligned AS\'s loaded:', len(d)

  # compare to AS's from piped-in SAM
  count = worse = better = equal = novel = notTrim = 0
  f = sys.stdin
  fOut = open(args[1], 'w')
  for line in f:
    if line[0] == '@': continue
    spl = line.rstrip().split('\t')
    if int(spl[1]) & 0x100: continue  # skip secondary
    aScore = int(getTag('AS', spl[11:]))
    if spl[0] in d:
      if aScore > d[spl[0]]:
        #print 'Worse:', spl[0], aScore, d[spl[0]]
        worse += 1
      elif aScore < d[spl[0]]:
        #print 'Better:', spl[0], aScore, d[spl[0]]
        fOut.write(spl[0] + '\tbetter\n')
        better += 1
      else:
        #print 'Equal:', spl[0], aScore, d[spl[0]]
        equal += 1
      del d[spl[0]]
    else:
      notTrim += 1
    count += 1
    if count % 10000000 == 0:
      print count

  for read in d:
    fOut.write(read + '\tnovel\n')
    #print 'Novel:', read, d[read]
    novel += 1
  fOut.close()

  print 'Reads in old:', count
  print '  Not trimmed:', notTrim
  print '  Worse AS:', worse
  print '  Same AS:', equal
  print '  Better AS:', better
  print 'New alignments:', novel

if __name__ == '__main__':
  main()
