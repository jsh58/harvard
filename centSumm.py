#!/usr/bin/python

# JMG 1/2018

# Producing a summary from centifuge's kraken-style report.

import sys

def loadScores(f):
  rank = {'D':0, 'P':1, 'C':2, 'O':3, 'F':4, 'G':5, 'S':6}
  score = {} # dict of scores (% abundance)
  hier = {}  # dict of taxonomic hierarchies
  temp = []  # temp copy of hierarchy

  for line in f:
    spl = line.split('\t')
    if spl[3] == '-' and spl[5][0] == ' ':
      # skip non-canonical levels
      continue
    taxon = spl[5].strip()
    if spl[3] not in rank:
      score[taxon] = spl[0]
      hier[taxon] = [taxon]
      temp = []   # reset
      continue

    # save taxon to temp, clear higher values
    level = rank[spl[3]]
    if level + 1 >= len(temp):
      temp = temp[:level]
    temp.append(taxon)

    score[taxon] = spl[0]
    hier[taxon] = temp[:]

  return score, hier

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python %s  <kreport>  <out>\n' % sys.argv[0])
    sys.exit(-1)

  f = open(args[0], 'rU')
  score, hier = loadScores(f)

  for s in sorted(score, key=score.get, reverse=True):
    print s, score[s], hier[s]
    raw_input()

if __name__ == '__main__':
  main()
