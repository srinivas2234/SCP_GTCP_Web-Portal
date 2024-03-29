"""The main program that runs gSpan."""
# -*- coding=utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys
import codecs
import collections
import copy
import itertools
import time
import pickle
s=sys.argv[3]
s=s.split('/')
p=sys.argv[2]
s=s[2].split('.')
s=s[0]
si=sys.argv[4]
mins=sys.argv[2]
fk=open(str(si)+"_"+str(p)+"_"+str(s)+"_sid_stats.txt",'a+')

from .gspan import nsg,flis1,flis2,flis3,timez

from .config import parser
from .gspan import gSpan
graph_cnt=0
# print(min_no_vertices)
def read_graphs(FLAGS=None):
    global graph_cnt
    if FLAGS is None:
        FLAGS, _ = parser.parse_known_args(args=sys.argv[1:])

    database_file_name=FLAGS.database_file_name
    with codecs.open(database_file_name, 'r', 'utf-8') as f:
        lines = [line.strip() for line in f.readlines()]
        for i, line in enumerate(lines):
            cols = line.split(' ')
            if cols[0] == 't':
                graph_cnt += 1
    return graph_cnt
graph_cnt=read_graphs()
print(graph_cnt)
min_sup=int(graph_cnt)*float(mins)
def read_si():

    with codecs.open(si, 'r', 'utf-8') as f:
        lines = [line.strip() for line in f.readlines()]
        si_count=0
        si_vert=collections.defaultdict(dict)
        si_degree=collections.defaultdict(list)
        si_ne=collections.defaultdict()
        no_of_edges=collections.Counter()

        for i, line in enumerate(lines):
            # si_degree=collections.defaultdict(list)
            cols=line.split(' ')
            # print(cols)
           
            if cols[0] == 't':
                if(si_count > 0):
                    si_ne[si_count] = si_degree.copy()
                    si_degree.clear()
                si_count += 1
            elif cols[-1] == '-1':
                break    
            elif cols[0] == 'v':
                si_vert[si_count][int(cols[1])]=int(cols[2])
            elif cols[0] == 'e':
                no_of_edges[si_count] += 1
                si_degree[int(cols[1])].append(si_vert[si_count][int(cols[2])])
                si_degree[int(cols[2])].append(si_vert[si_count][int(cols[1])])
  
    return si_vert,si_ne,no_of_edges
    
def main(FLAGS=None):
    """Run gSpan."""
    print(sys.argv[0],sys.argv[1])
    if FLAGS is None:
        FLAGS, _ = parser.parse_known_args(args=sys.argv[1:])

    if not os.path.exists(FLAGS.database_file_name):
        print('{} does not exist.'.format(FLAGS.database_file_name))
        sys.exit()
    si_vert,si_ne,no_of_edges=read_si()
    for i in range(1,(len(si_vert)+1)):
        # print("kjjjjas",i)
        # print(len(si_vert[i]))
        # print("asa",min_sup)
        min_no=len(si_vert[i])
        max_no=len(si_vert[i])
        gs = gSpan(
            si_id=i,
            noe=no_of_edges[i],
            vert=si_vert[i],
            ne=si_ne[i],
            database_file_name=FLAGS.database_file_name,
            min_support=int(min_sup),
            # min_num_vertices=FLAGS.lower_bound_of_num_vertices,
            min_num_vertices=min_no,
            # max_num_vertices=FLAGS.upper_bound_of_num_vertices,
            max_num_vertices=max_no,
            max_ngraphs=FLAGS.num_graphs,
            is_undirected=(not FLAGS.directed),
            verbose=FLAGS.verbose,
            visualize=FLAGS.plot,
            where=FLAGS.where
        )

        gs.run()
        gs.time_stats()
    
    fp=sorted(flis3)
    f1=open(str(si)+"_"+str(p)+"_"+str(s)+"_results.txt",'w')
    sum=0
    
    
    for i in fp:
        sum+=len(flis3[i])
    print("grpah",graph_cnt)
    for i in range(graph_cnt):
        # print(graph_cnt,i,flis3[i])
        if(flis3[i]==[]):
            pass
        for j in range(len(flis3[i])):
            f1.write(str(flis3[i][j]))
            if(j != len(flis3[i])-1):
                f1.write(" ")
        f1.write("\n")
        
    f1.close()
    from .gspan import list_fs
    # f=open(str(si)+"_"+str(p)+"_"+str(s)+"_stats.txt",'w')
    
    # print(nsg,file=f)
    # df=sorted(flis2)
    # fg=open(str(si)+"_"+str(p)+"_"+str(s)+"si_subgraphsId.txt",'w')
    # for i in df:
    #     fg.write(str(i))
    #     fg.write(": ")
    #     fg.write(str(flis2[i])[1:-1])
    #     fg.write("\n")
    # fg.close()
    avg=sum/graph_cnt
    fk.write("avg size of transaction: ")
    fk.write(str(avg))
    fk.write("\n")
    fg=open("gSpan_FSM_"+str(s)+"_stats.txt",'w')
    fg.write(" minrf : "+str(min_sup)+'\n')
    fg.write(" no_of_subgraphs : "+str(sum)+'\n')
    fg.write(" total execution time : "+str(0)+'\n')
    fg.write(" Avg size of flat_trans : "+str(avg)+"\n")
    fg.close()
    print("result.txt",len(list_fs))
    with open('result.txt','wb') as fp:
        pickle.dump(list_fs,fp)
    return gs


if __name__ == '__main__':
    main()
