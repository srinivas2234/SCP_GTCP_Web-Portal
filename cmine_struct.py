import time
import sys
import os
import copy
import shutil
from PIL import Image
from werkzeug import datastructures
import numpy as np
global One_freqItems
global f
global scp
class cmine():
    def __init__(self, minRF, minCS, maxOR, inpfile, outfile,datasetname):
        global f
        global scp
        self.minRF = float(minRF)
        self.minCS = float(minCS)
        self.minCS1=minCS
        self.maxOR = float(maxOR)
        self.inpfile = inpfile
        self.outfile = outfile
        self.nots = self.getlines(inpfile)
        '''print(self.nots)
        print(self.minCS)
        print(self.minCS1)'''
        self.dataset=datasetname
        #self.fout = open(outfile,'w')
        self.NOk = []
        self.frequencies=[]
        self.noofSCP=0
        [self.items, self.database] = self.dbscan(inpfile)
        sorteditems = sorted(self.items.items(), key = lambda a: (-a[1],a[0]))
        '''print("sorted")
        print(sorteditems)'''
        mintracs = float(self.minRF) * 1.0 * float(self.nots)
        print("self.nots")
        print(self.nots)
        freqitems=[]
        f=open("./"+datasetname+"_"+str(self.minCS1)+"_"+str(self.minRF)+"_"+str(self.maxOR),'w')
        f.write("Coverage"+" "+"1"+'\n')
        f.write('\n')

        for i in sorteditems:
            '''print(i)
            print(i[1]>=mintracs)'''
            if(i[1]>=mintracs):
                
                freqitems.append(i)
        
        '''freqitems = filter(lambda x: (x[1] >= mintracs), sorteditems)'''
        '''one_size_coverage = filter(lambda x: (x[1] >= minCS*self.nots),freqitems)'''
        one_size_coverage=[]
        c=[]
        count=0
        scp=0
        for i in freqitems:
            '''print(i)
            print(i[1],float(self.minCS))'''
            print(i[1],float(minCS)*float(self.nots))
            if(i[1]>=float(minCS)*float(self.nots)):
                scp=scp+1
                one_size_coverage.append(i[0])
                count=count+1
                f.write("i"+" "+str(scp)+'\n')
                f.write("cs"+" "+str(float(i[1])/self.nots)+'\n')
                f.write("or"+" "+str(0)+'\n')
                f.write('\n')
                c.append(i[0])
        self.frequencies.append(c)
        f.write('\n')
        f.write('\n')
        '''print("one_size_coverage")'''

        print(one_size_coverage)
        self.freqitems = map(lambda x: x[0], freqitems)  
        #self.noofFreqItems=copy.copy(len(list(self.freqitems)))
        
        for l in self.freqitems:
            '''print("fe")'''
            self.NOk.append([l])
        '''print(self.NOk)'''
        outputfile="./scp_images_"+datasetname+"_"+str(self.minCS1)
        if(os.path.exists(outputfile)):
            shutil.rmtree(outputfile)
        os.mkdir(outputfile)
        if(os.path.exists(outputfile+"/Coverage_1")):
            shutil.rmtree(outputfile+"/Coverage_1")
        os.mkdir(outputfile+"/Coverage_1")
        
        
        self.source_dir="./"+datasetname+"data_2_"+str(self.minCS1)
        self.dest_dir="./scp_images_"+datasetname+"_"+str(self.minCS1)
        for i in one_size_coverage:
            '''print("hi")'''
            self.noofSCP=self.noofSCP+1
            '''print("k")'''
            shutil.copy(self.source_dir+"/"+str(i)+".svg",self.dest_dir+"/Coverage_1/")

            #self.fout.write("['"+str(i[0])+"']\n")
       
       
    def get_overlapratio_cs(self, pattern):
        ovr_nume_1=set()
        for i in pattern[:-1]:
            for j in range(len(self.database)):
                if i in self.database[j]:
                    ovr_nume_1.add(j)
        ovr_deno = set()
        for j in range(len(self.database)):
            if pattern[-1] in self.database[j]:
                ovr_deno.add(j)
        cs_nume = ovr_nume_1.union(ovr_deno)
        ovr_nume = ovr_nume_1.intersection(ovr_deno)
        # print ovr_nume,ovr_deno,ovr_nume_1
        return len(ovr_nume)*1.0/len(ovr_deno),len(cs_nume)*1.0/self.nots

    def expand(self):
        global f
        global scp
        cnt = 0
        cnt1 = 0
        length = 1
        coverage=2
        while len(self.NOk)>0:
            os.mkdir(self.dest_dir+"/Coverage_"+str(coverage))
            f.write("Coverage"+" "+str(coverage)+'\n')
            f.write('\n')
            temp_NOk = self.NOk
            
            self.NOk = []
            count=1
           
            for i in range(len(temp_NOk)):
                c=[]
                for j in range(i+1, len(temp_NOk)):
                    cnt += 1
                    if temp_NOk[i][:-1] == temp_NOk[j][:-1]:
                        cnt1 += 1
                        newpattern = temp_NOk[i] + [temp_NOk[j][-1]]
                        '''print("new")
                        print(newpattern)'''
                        # print "newpattern",newpattern
                        overlapratio,cs = self.get_overlapratio_cs(newpattern)
                        # print newpattern,overlapratio,cs
                        if float(overlapratio) <= float(self.maxOR):
                            self.NOk.append(newpattern)
                           
                            if float(cs) >= float(self.minCS):
                                scp=scp+1
                                print("coverage spport")
                                print(float(cs))
                                f.write("i"+" "+str(scp)+'\n')
                                f.write("cs"+" "+str(round(float(cs),2))+'\n')
                                f.write("or"+" "+str(round(overlapratio,2))+'\n')
                                f.write('\n')
                                opened_images=[]
                                for k in newpattern:
                                    '''print("true")'''
                                    opened_images.append(Image.open(self.source_dir+"/"+str(k)+".png"))
                                widths,heights=zip(*(k.size for k in opened_images))
                                total_width=sum(widths)
                                max_height=max(heights)
                                new_img = np.zeros(shape=(max_height, total_width, 3))
                                new_image=Image.new('RGB',(total_width,max_height))
                                x_offset=0
                                for k in opened_images:
                                    new_image.paste(k,(x_offset,0))
                                    x_offset=x_offset+k.size[0]
                                new_image.save(self.dest_dir+"/Coverage_"+str(coverage)+"/"+str(count)+".png")
                               
                                self.noofSCP=self.noofSCP+1
                                count=count+1
                                # print "This is coverage pattern",newpattern,cs,overlapratio
                                #self.fout.write(str(newpattern)+"\n")
                    else:
                        break
                self.frequencies.append(c)
            length += 1
            f.write('\n')
            f.write('\n')
            coverage=coverage+1
        '''print("total read, total calculate_or_cs",cnt,cnt1)'''
        #self.fout.close()
        f.close()
        return cnt1,self.noofSCP,length



    def dbscan(self,inputfile):
        f=open(inputfile,'r')
        a = {}
        database = []
        for row in f:
            row = row.rstrip('\n')
            row = row.split(" ")
            if len(row[-1]) == 0:
                row.pop()
            database.append(row)
            for j in row:
                if j in a:
                    a[j] += 1
                else:
                    a[j] = 1
        return [a,database]

    def getlines(self,filename):
        with open(filename,"r") as f:
            ''''return sum(1 for _ in f)'''
            a=0
            for i in f:
                if i.strip():
                    a=a+1
            return a
t1 = time.time()
print("came_t_struct")
minRF = float(sys.argv[1])
minCS =sys.argv[2]
maxOR = float(sys.argv[3])
structure=sys.argv[4]
datasetname = sys.argv[5]
inpfile = structure+".txt_"+str(minRF)+"_"+datasetname+"_results.txt"
#datasetname = "graph5_test2.txt"
#filepath = "C:\\Users\\user\\Documents\\CoverageGraph\\dataset1\\GraphData\\"
outfile = "Results.txt"
'''print(outfile)'''
obj=cmine(minRF, minCS, maxOR, inpfile, outfile,datasetname)
t3=time.time()
print( "data read",str(t3-t1))
candidate_patterns,SCPs,max_graphs = obj.expand()
t2 = time.time()
print( "process done",str(t2-t1))
f = open("./Results.txt",'w')
f.write("execution_time:"+str(t2-t1)+'\n')
f.write("Number of Candidate Patterns:"+str(candidate_patterns)+'\n')
f.write("Number of Coverage Patterns:"+str(SCPs)+'\n')
f.write("Max_graphs:"+str(max_graphs)+'\n')
f.close()
