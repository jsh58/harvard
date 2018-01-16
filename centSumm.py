#!/usr/bin/python

# JMG 1/2018

# Producing a summary from centifuge's kraken-style report.

import sys

def printOutput(f, score, hier):
  '''
  Print results.
  '''
  # sort results
  res = sorted(score, key=score.get, reverse=True)

  # for sorting same-scored taxa by depth
  def lenTup(tup):
    return len(tup[2])

  i = 0
  while i < 20:
    # save taxa with identical scores
    tRes = []
    tScore = score[res[i]]
    while tScore == score[res[i]]:
      tRes.append((res[i], score[res[i]], hier[res[i]]))
      i += 1

    # print results
    for s in sorted(tRes, key=lenTup):
      f.write('%.2f' % s[1])
      for j in range(len(s[2])):
        f.write('\t' + ' ' * j + s[2][j] + '\n')


def loadScores(f):
  '''
  Save scores and taxonomic hierarchies from a
    Centrifuge report.
  '''
  rank = {'D':0, 'P':1, 'C':2, 'O':3, 'F':4, 'G':5, 'S':6}
  score = {} # dict of scores (% abundance)
  hier = {}  # dict of taxonomic hierarchies
  temp = []  # temp copy of hierarchy

  for line in f:
    spl = line.split('\t')
    if spl[3] == '-': # and spl[5][0] == ' ':
      # skip non-canonical levels
      continue
    taxon = spl[5].strip()
    if spl[3] not in rank:
      score[taxon] = float(spl[0])
      hier[taxon] = [taxon]
      temp = []   # reset
      continue

    # save taxon to temp, clear higher values
    level = rank[spl[3]]
    if level + 1 <= len(temp):
      temp = temp[:level]
    temp.append(taxon)

    score[taxon] = float(spl[0])
    hier[taxon] = temp[:]

  return score, hier

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python %s  <kreport>  <out>\n' % sys.argv[0])
    sys.exit(-1)

  f = open(args[0], 'rU')
  score, hier = loadScores(f)

  fOut = open(args[1], 'w')
  printOutput(fOut, score, hier)

if __name__ == '__main__':
  main()
