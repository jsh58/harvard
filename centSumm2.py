#!/usr/bin/python

# JMG 1/2018

# Producing a summary from centifuge's kraken-style report.
# Version 2: html output

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
  Node: contains taxon name, taxon number,
    list of child nodes, parent node,
    score (percent of the sample), and
    count (number of reads).
  '''
  def __init__(self, parent, name, taxon, score, count):
    self.child = []
    self.parent = parent
    self.name = name
    self.taxon = taxon
    self.score = float(score)
    self.count = int(count)

def printFooter(f, num, version, date):
  '''
  Print centrifuge/nt info, etc. for the output.
  '''
  f.write('<p>\nTop %d taxa, produced by ' % num \
    + '<a href="http://www.ccb.jhu.edu/software/centrifuge/manual.shtml">' \
    + 'Centrifuge</a>')
  if version:
    f.write(' (version %s)' % version)
  f.write(', querying the NCBI ' \
    + '<a href="https://www.ncbi.nlm.nih.gov/nucleotide">nt</a> database')
  if date:
    f.write(' (downloaded %s)' % date)
  f.write('.</p>\n')

def printLevel(f, n, level, cutoff):
  '''
  Print results for a node (if its count meets cutoff).
    Continue printing for children nodes (recursively).
  '''
  if n.count >= cutoff:  # or level == 0: # to include all children of root
    f.write('  <tr>\n' \
      + '    <td align="right">%.2f&emsp;</td>\n' % n.score \
      + '    <td>%s%s</td>\n' % (level * 2 * '&emsp;', n.name) \
      + '  </tr>\n')
  for m in n.child:
    printLevel(f, m, level + 1, cutoff)

def printOutput(f, unclass, root, num, cutoff, version, date):
  '''
  Begin printing results (header and unclassified).
    Start recursive tree printing.
  '''
  f.write('''<h2>Taxonomy Analysis</h2>
<table style="width:100%;border:1px solid;">
  <tr>
    <th align="right" width=10%>Percent&emsp;</th>
    <th align="left">Taxon</th>
  </tr>
''')
  f.write('  <tr>\n' \
    + '    <td align="right">%.2f&emsp;</td>\n' % unclass \
    + '    <td>unclassified</td>\n' \
    + '  </tr>\n')

  # print tree
  for n in root.child:
    printLevel(f, n, 0, cutoff)
  f.write('</table>\n')

  printFooter(f, num, version, date)

def findCutoff(score, x):
  '''
  Determine threshold for top x scores.
  '''
  if x >= len(score):
    return 0
  res = sorted(score, reverse=True)
  return res[x-1]

def checkTree(node, taxon):
  '''
  Search node and its children for given taxon.
  '''
  if node.taxon == taxon:
    return node
  for n in node.child:
    parent = checkTree(n, taxon)
    if parent:
      return parent
  return None

def loadScores(f, d):
  '''
  Create taxonomic tree (including scores) from a
    Centrifuge report.
  '''
  rank = 'DKPCOFGS'
  unclass = 0.0  # 'unclassified' score
  root = Node(None, 'root', '1', -1, -1) # root of tree
  temp = root    # pointer to previous node
  score = []     # list of scores (read counts)

  for line in f:
    spl = line.split('\t')
    if len(spl) < 6:
      sys.stderr.write('Error! Improperly formatted ' \
        + 'centrifuge-kreport file\n')
      sys.exit(-1)

    # skip non-canonical levels
    if spl[3] == '-':
      if spl[4] not in ['12908', '28384']:
        # make exception for these '-' taxa
        #   (see loadTax(), below)
        continue
    elif spl[3] == 'U':
      # save unclassified value automatically
      unclass = float(spl[0])
      continue
    elif spl[3] not in rank:
      sys.stderr.write('Warning! Unknown taxonomic rank:' \
        + ' %s\n' % spl[3] + '  ' + line)
      continue

    # find parent node in hierarchy
    if spl[4] not in d:
      sys.stderr.write('Warning! Unknown taxon: %s\n' % spl[4])
      continue
    parent = None
    # check current branch first
    while temp != None:
      if temp.taxon == d[spl[4]]:
        parent = temp
        break
      temp = temp.parent
    # if not found, check whole tree from root
    if parent == None:
      parent = checkTree(root, d[spl[4]])
    if parent == None:
      sys.stderr.write('Warning! Cannot find parent for ' \
        + 'taxon %s\n' % spl[4])
      continue

    # save to tree
    n = Node(parent, spl[5].strip(), spl[4], spl[0], spl[1])
    parent.child.append(n)
    temp = n

    # save score
    score.append(int(spl[1]))

  return unclass, root, score

def findParent(d, taxon):
  '''
  Find parent taxon (recursively).
  '''
  if taxon not in d or d[taxon][0] not in d:
    return None
  if d[ d[taxon][0] ][1]:
    return d[taxon][0]
  return findParent(d, d[taxon][0])

def loadTax(f):
  '''
  Load parents of each taxon, keeping only
    canonical (DKPCOFGS) ones.
    (include 1 -> 'root'
         12908 -> 'unclassified sequences'
         28384 -> 'other sequences')
  '''
  # load immediate parent taxa to dict
  temp = {}
  for line in f:
    spl = line.split('|')
    if len(spl) < 3:
      sys.stderr.write('Error! Improperly formatted tree file\n')
      sys.exit(-1)

    level = False
    if spl[2].strip() in ['superkingdom', 'kingdom', 'phylum', \
        'class', 'order', 'family', 'genus', 'species'] \
        or spl[0].strip() in ['1', '12908', '28384']:
      level = True
    temp[spl[0].strip()] = (spl[1].strip(), level)

  # save canonical parent taxa to dict
  d = {}
  for taxon in temp:
    if temp[taxon][1]:
      parent = findParent(temp, taxon)
      if parent:
        d[taxon] = parent
      else:
        sys.stderr.write('Warning! Cannot find parent ' \
          + 'of taxon %s\n' % taxon)

  return d

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 3:
    sys.stderr.write('Usage: python %s  ' % sys.argv[0] \
      + '<kreport>  <taxTree>  <out> \ \n' \
      + '    [<num>]  [<version>  <date>]\n' \
      + '  <num>      Number of taxa to print (def. 20)\n' \
      + '  <version>  Version of centrifuge\n' \
      + '  <date>     Date of nt download\n')
    sys.exit(-1)

  # load tax tree
  fTax = openRead(args[1])
  d = loadTax(fTax)
  if fTax != sys.stdin:
    fTax.close()

  # load scores and create taxonomic tree
  fIn = openRead(args[0])
  unclass, root, score = loadScores(fIn, d)
  if fIn != sys.stdin:
    fIn.close()

  # find cutoff score for top n taxa
  num = 20
  if len(args) > 3:
    num = int(args[3])
  cutoff = findCutoff(score, num)

  # load centrifuge version, date of nt download
  version = date = ''
  if len(args) > 5:
    version = args[4]
    date = args[5]

  # print output
  fOut = openWrite(args[2])
  printOutput(fOut, unclass, root, num, cutoff, version, date)
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
