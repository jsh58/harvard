# JMG 2/21/18
# finding adapters in a sequence

import sys

def revComp(dna):
  '''
  Reverse-complements the given DNA sequence.
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

def createSeeds(d):
  '''Create 20-mer seeds of given dict of sequences'''
  seed = {}
  seedRC = {}
  for r in d:
    seed[r] = []
    seedRC[r] = []
    seq = d[r]
    for i in range(len(seq) - 19):
      seed[r].append(seq[i:i+20])
      seedRC[r].append(revComp(seq[i:i+20]))
  return seed, seedRC

def loadSeqs(f):
  '''Return dict of fasta sequences.'''
  d = {}
  seq = ''
  head = ''
  for line in f:
    if line[0] == '>':
      if seq:
        d[head] = seq
        seq = ''
      head = line.rstrip().split(' ')[0][1:]
    else:
      seq += line.rstrip()
  if seq:
    d[head] = seq
  return d

def findMatches(query, seed, seedRC):
  '''Find adapter matches.'''
  ret = {}
  for q in query:
    seq = query[q]
    match = {}
    for s in seed:
      for i in range(len(seed[s])):
        # save list of seed matches
        pos = seq.find(seed[s][i])
        if pos != -1:
          #print seed[s][i], pos
          if (s, False) not in match:
            match[(s, False)] = [i]
          else:
            match[(s, False)].append(i)
        pos = seq.find(seedRC[s][i])
        if pos != -1:
          #print seedRC[s][i], '(rc)', pos
          if (s, True) not in match:
            match[(s, True)] = [i]
          else:
            match[(s, True)].append(i)

    # find longest set of seed matches
    long = {}
    for m in match:
      length = 0
      longest = [0, 0]
      val = match[m][0]
      i = 1
      while i < len(match[m]):
        if match[m][i] != match[m][i-1] + 1:
          if match[m][i-1] - val > length:
            length = match[m][i-1] - val
            longest = [val, match[m][i-1]]
          val = match[m][i]
        i += 1
      if i == 1:
        length = 1
        longest = [match[m][0], match[m][0]]
      elif match[m][i-1] - val > length:
        length = match[m][i-1] - val
        longest = [val, match[m][i-1]]
      long[m] = [length, longest[0], longest[1] + 20]

    maxVal = 0
    for m in sorted(long, key=lambda k: long[k][0], reverse=True):
      if not maxVal:
        maxVal = long[m][0]
      if long[m][0] != maxVal:
        break
      if q not in ret:
        ret[q] = []
      ret[q].append((m[0], m[1], long[m][1], long[m][2]))

  return ret

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python findSeq.py  <seq.fa>  <adapters.fa>\n')
    sys.exit(-1)

  f1 = open(args[1], 'rU')
  adapters = loadSeqs(f1)
  f1.close()
  seed, seedRC = createSeeds(adapters)

  f2 = open(args[0], 'rU')
  query = loadSeqs(f2)

  match = findMatches(query, seed, seedRC)
  for q in match:
    head = '>'
    for i in range(len(match[q])):
      m = match[q][i]
      if i:
        head += ';'
      else:
        seq = adapters[m[0]][m[2]:m[3]]
        if m[1]:
          seq = revComp(seq)

      head += '%s_%d-%d' % (m[0], m[2], m[3])
      if m[1]:
        head += 'rc'

    sys.stdout.write(head + '\n' + seq + '\n')

    pos = query[q].find(seq)
    sys.stdout.write('>%s_%d-%d\n' % (q, pos, pos + len(seq)))
    sys.stdout.write(query[q][pos:pos+len(seq)] + '\n')

if __name__ == '__main__':
  main()
