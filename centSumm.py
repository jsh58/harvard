#!/usr/bin/python

# JMG 1/2018

# Producing a summary from centifuge's kraken-style report.

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

class Node:
  '''
  Node: contains taxon name, rank (DKPCOFGS),
    and score (percent of the sample).
  '''
  def __init__(self, taxon, rank, score):
    self.child = []
    self.taxon = taxon
    self.rank = rank
    self.score = float(score)

  #def children(self):
  #  ret = []
  #  for c in self.child:
  #    ret.append(c.taxon)
  #  return ret

def printLevel(f, n, level, cutoff):
  '''
  Print results (recursively).
  '''
  if n.score >= cutoff:
    f.write('%6.2f\t' % n.score + '  ' * level \
      + n.taxon + '\n')
  for m in n.child:
    printLevel(f, m, level + 1, cutoff)

def printOutput(f, unclass, root, cutoff):
  '''
  Begin printing results.
  '''
  f.write('percent\ttaxon\n')
  f.write('%6.2f\tunclassified\n' % unclass)
  for n in root.child:
    printLevel(f, n, 0, cutoff)

def findCutoff(score, n):
  '''
  Determine threshold for top n scores.
  '''
  if n >= len(score):
    return min(score)
  res = sorted(score, reverse=True)
  return res[n-1]

def loadScores(f):
  '''
  Save scores and taxonomic hierarchies from a
    Centrifuge report.
  '''
  rank = 'DKPCOFGS'
  unclass = 0.0  # 'unclassified' score
  root = Node('root', 'X', -1)
  temp = [root]  # temp copy of hierarchy
  score = []     # list of scores

  for line in f:
    spl = line.split('\t')

    # skip non-canonical levels
    if spl[3] == '-':
      if spl[5] in ['  other sequences\n', \
          '  unclassified sequences\n']:
        # add special node (for 'other sequences' and
        #   'unclassified sequences'); reset 'temp'
        n = Node(spl[5].strip(), 'X', spl[0])
        root.child.append(n)
        temp = temp[:1] + [n]
      continue
    elif spl[3] == 'U':
      # save unclassified value automatically
      unclass = float(spl[0])
      continue
    elif spl[3] not in rank:
      sys.stderr.write('Warning! Unknown taxonomic rank:' \
        + ' %s\n' % spl[3] + '  ' + line)
      continue

    # find placement of taxon in 'temp' hierarchy
    if spl[3] == 'D':
      temp = temp[:1]  # automatically reset for 'D'
    else:
      found = False
      for i in range(1, len(temp)):
        # find same rank already in temp
        if spl[3] == temp[i].rank:
          temp = temp[:i]
          found = True
          break
      if not found:
        # if same rank not found, look for parent rank
        par = rank[rank.find(spl[3]) - 1]
        for i in range(1, len(temp)):
          if par == temp[i].rank:
            temp = temp[:i + 1]
            break

    # save to tree
    taxon = spl[5].strip()
    n = Node(taxon, spl[3], spl[0])
    parent = temp[-1]
    parent.child.append(n)
    temp.append(n)

    # save score
    score.append(float(spl[0]))

  return unclass, root, score

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python %s  ' % sys.argv[0] \
      + '<kreport>  <out>  [<num>]\n' \
      + '  <num>    Number of taxa to print (def. 20)\n')
    sys.exit(-1)

  # load scores and taxonomic tree
  f = openRead(args[0])
  unclass, root, score = loadScores(f)
  if f != sys.stdin:
    f.close()

  # find cutoff score for top n taxa
  num = 20
  if len(args) > 2:
    num = int(args[2])
  cutoff = findCutoff(score, num)

  # print output
  fOut = openWrite(args[1])
  printOutput(fOut, unclass, root, cutoff)
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
