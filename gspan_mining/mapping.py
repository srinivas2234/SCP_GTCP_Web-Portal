
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys

import codecs
import collections
import copy
import itertools
import time
mapping=collections.defaultdict(dict)
edge_map=collections.defaultdict(dict)
s=sys.argv[5]
s=s.split('/')

s=s[2].split('.')
s=s[0]
def read_si():
    with codecs.open(str(s)+'_mapping.txt', 'r', 'utf-8') as f:
        global mapping
        lines = [line.strip() for line in f.readlines()]
        for i, line in enumerate(lines):
            cols=line.split(' ')    
            mapping[cols[1]]=cols[0]
    
    with codecs.open(str(s)+'_edge_mapping.txt', 'r', 'utf-8') as g:
        global edge_map
        lines1 = [linea.strip() for linea in g.readlines()]
        for i,linea in enumerate(lines1):
            colsa=linea.split(' ')
            edge_map[colsa[0]]=colsa[1:]
    # print("ed",edge_map)
    return mapping,edge_map
read_si()
# print(edge  _map)
