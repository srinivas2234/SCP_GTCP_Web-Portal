import time
import sys
import os
import copy
import shutil
from PIL import Image
from werkzeug import datastructures
import numpy as np
import cv2
import pickle
import base64
import matplotlib.pyplot as plt
global One_freqItems
global f
global scp
global cp_data
class cmine():
    def __init__(self, minRF, minCS, maxOR, inpfile, outfile,datasetname,fsubgraphs_data):
        global f
        global scp
        global cp_data
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
        self.fsubgraphs_data=fsubgraphs_data
        #self.fout = open(outfile,'w')
        self.NOk = []
        self.frequencies=[]
        self.noofSCP=0
        self.total=0
        [self.items, self.database] = self.dbscan(inpfile)
        sorteditems = sorted(self.items.items(), key = lambda a: (-a[1],a[0]))
        
        mintracs = float(self.minRF) * 1.0 * float(self.nots)
        freqitems = filter(lambda x: (x[1] >= mintracs), sorteditems)
       
        '''for i in freqitems:
            print(i[1],self.minCS*self.nots)
            print((i[1]==self.minCS*self.nots))'''
        
        one_size_coverage = filter(lambda x: (x[1] >= self.minCS*self.nots),freqitems)
        f=open("./"+datasetname+"_"+str(self.minCS1)+"_"+str(self.minRF)+"_"+str(self.maxOR),'w')
        f.write("Coverage"+" "+"1"+'\n')
        f.write('\n')

        '''for i in sorteditems:
            print(i)
            print(i[1]>=mintracs)
            if(i[1]>=mintracs):
                
                freqitems.append(i)
        
        freqitems = filter(lambda x: (x[1] >= mintracs), sorteditems)
        print(freqitems)
        one_size_coverage = filter(lambda x: (x[1] >= minCS*self.nots),freqitems)
        one_size_coverage=[]
        c=[]'''
        count=1
        for i in range(0,len(one_size_coverage)):
            self.noofSCP=self.noofSCP+1
            f.write("i"+" "+str(self.noofSCP)+'\n')
            f.write("cs"+" "+str(round(float(one_size_coverage[i][1])/float(self.nots)))+'\n')
            f.write("or"+" "+str(self.maxOR)+'\n')
            count=count+1
            f.write('\n')
        f.write('\n')
        f.write('\n')
        '''print("one_size_coverage")'''
        self.source_dir="./"+datasetname+"data_2_"+str(float(self.minRF))
        self.dest_dir="./scp_images_"+datasetname+"_"+str(self.minCS1)
        self.freqitems = map(lambda x: x[0], freqitems) 
        self.source_dir_graph_smiles="./smiles_graphs_"+datasetname+"_"+str(self.minRF) 
        self.dest_dir_graph_smiles="./scp_smiles_graphs_"+datasetname+"_"+str(self.minCS1)
        if(os.path.exists(self.dest_dir_graph_smiles)):
            shutil.rmtree(self.dest_dir_graph_smiles)
        os.mkdir(self.dest_dir_graph_smiles)
        #self.noofFreqItems=copy.copy(len(list(self.freqitems)))
        
        for l in self.freqitems:
            '''print("fe")'''
            self.NOk.append([l])
        '''print(self.NOk)'''
        outputfile="./scp_images_"+datasetname+"_"+str(self.minCS1)
        if(os.path.exists(outputfile)):
            shutil.rmtree(outputfile)
        os.mkdir(outputfile)
        if(len(one_size_coverage)!=0):
            if(os.path.exists(outputfile+"/Coverage_1")):
                shutil.rmtree(outputfile+"/Coverage_1")
            os.mkdir(outputfile+"/Coverage_1")
            if(os.path.exists(self.dest_dir_graph_smiles+"/Coverage_1")):
                shutil.rmtree(self.dest_dir_graph_smiles+"/Coverage_1")
            os.mkdir(self.dest_dir_graph_smiles+"/Coverage_1")
        
        cp_data=[]
        '''print("frequent")'''
        
        for i in range(0,len(one_size_coverage)):
            h={"id":str(self.total),"or":str(self.maxOR),"cs":str(round((float(one_size_coverage[i][1])/float(self.nots)),2)),"size":1,"graphs":[],"string":''}
            h["graphs"].append(fsubgraphs_data[int(one_size_coverage[i][0])-1])
            h["string"]=h["string"]+fsubgraphs_data[int(one_size_coverage[i][0])-1]["string"]
            shutil.copy(self.source_dir_graph_smiles+"/"+str(one_size_coverage[i][0])+".svg",self.dest_dir_graph_smiles+"/Coverage_1/")
            self.total=self.total+1
            #print(str(self.source_dir_graph_smiles+"/"+str(one_size_coverage[i][0])+".png"))
            print(one_size_coverage[i][0])
            with open(self.source_dir_graph_smiles+"/"+str(one_size_coverage[i][0])+".svg","rb") as img:
                k = str(base64.b64encode(img.read()))
               
                if(i==0):
                    print(k)
                #print(count)
                #print(no_of_coverage)
                h["smile_graph_image"]=[k]
            
            cp_data.append(h)
            #shutil.copy(self.source_dir+"/"+str(one_size_coverage[i][0])+".png",self.dest_dir+"/Coverage_1/")
        

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
        global cp_data
        cnt = 0
        cnt1 = 0
        length = 1
        coverage=2
        self.source_dir="./"+datasetname+"data_2_"+str(float(minRF))
        self.dest_dir="./scp_images_"+datasetname+"_"+str(self.minCS1)
        while len(self.NOk)>0:
            if(os.path.exists(self.dest_dir+"/Coverage_"+str(coverage))):
                shutil.rmtree(self.dest_dir+"/Coverage_"+str(coverage))
            if(os.path.exists(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage))):
                shutil.rmtree(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage))
            f.write("Coverage"+" "+str(coverage)+'\n')
            f.write('\n')
            temp_NOk = self.NOk
            
            self.NOk = []
            count=1
            if(len(temp_NOk)!=0):
                coverage_n={"coverage":coverage,"graphs":[]}
            for i in range(len(temp_NOk)):
                c=[]
                for j in range(i+1, len(temp_NOk)):
                    cnt += 1
                    if temp_NOk[i][:-1] == temp_NOk[j][:-1]:
                        cnt1 += 1
                        newpattern = temp_NOk[i] + [temp_NOk[j][-1]]
                       
                        # print "newpattern",newpattern
                        overlapratio,cs = self.get_overlapratio_cs(newpattern)
                        # print newpattern,overlapratio,cs
                        if float(overlapratio) <= float(self.maxOR):
                            self.NOk.append(newpattern)
                           
                            if float(cs) >= float(self.minCS):
                                if(not os.path.exists(self.dest_dir+"/Coverage_"+str(coverage))):
                                    os.mkdir(self.dest_dir+"/Coverage_"+str(coverage))
                                if(not os.path.exists(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage))):
                                    os.mkdir(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage))
                                self.noofSCP=self.noofSCP+1
                                
                                f.write("i"+" "+str(self.noofSCP)+'\n')
                                f.write("cs"+" "+str(round(float(cs),2))+'\n')
                                f.write("or"+" "+str(overlapratio)+'\n')
                                f.write('\n')
                                opened_images=[]
                                opened_smiles_images=[]
                                h={"id":str(self.total),"or":str(overlapratio),"cs":str(round(float(cs),2)),"size":len(newpattern),"graphs":[],"string":''}
                                '''print("pj")'''
                                h["smile_graph_image"]=[]
                                self.total=self.total+1
                                for k in newpattern:
                                    '''print("true")'''
                                    '''print(k)'''
                                    h["graphs"].append(fsubgraphs_data[int(k)-1])
                                    h["string"]=h["string"]+fsubgraphs_data[int(k)-1]["string"]
                                    opened_images.append(cv2.imread(self.source_dir+"/"+str(k)+".png"))
                                    #if(os.path.exists(self.source_dir_graph_smiles+"/"+str(k)+".svg")):
                                        #print("yes")
                                    with open(self.source_dir_graph_smiles+"/"+str(k)+".svg","rb") as img:
                                        g = str(base64.b64encode(img.read()))
                                        
                                        h["smile_graph_image"].append(g)
                                    # c=os.path.join(self.source_dir_graph_smiles,str(k)+".svg")
                                    # opened_smiles_images.append(Image.open("./smiles_graphs_papertoydata_0.3/1.svg"))
                                # h_min = min(img.shape[0] for img in opened_images)
                                # print(k)
                                # print(opened_smiles_images)
                                # h_min_smiles_images=min(img.shape[0] for img in opened_smiles_images)
                                '''widths,heights=zip(*(k.size for k in opened_images))
                                total_width=sum(widths)
                                max_height=max(heights)
                                new_img = np.zeros(shape=(max_height, total_width, 3))
                                new_image=Image.new('RGB',(total_width,max_height))
                                x_offset=0
                                for k in opened_images:
                                    new_image.paste(k,(x_offset,0))
                                    x_offset=x_offset+k.size[0]'''
                                # im_list_resize = [cv2.resize(img,(int(img.shape[1] * h_min / img.shape[0]),h_min), interpolation=cv2.INTER_CUBIC) for img in opened_images]
                                # im_list_resize_smiles=[cv2.resize(img,(int(img.shape[1] * h_min_smiles_images / img.shape[0]),h_min_smiles_images), interpolation=cv2.INTER_CUBIC) for img in opened_smiles_images]
                                # new_image=cv2.hconcat(im_list_resize)
                                # new_image_smiles=cv2.hconcat(im_list_resize_smiles)
                                #cv2.imwrite(self.dest_dir+"/Coverage_"+str(coverage)+"/"+str(count)+".png",new_image)
                                # cv2.imwrite(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage)+"/"+str(count)+".svg",new_image_smiles)
                                
                                '''new_image.save(self.dest_dir+"/Coverage_"+str(coverage)+"/"+str(count)+".png")'''
                                # with open(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage)+"/"+str(count)+".svg","rb") as img:
                                #     k = str(base64.b64encode(img.read()))
                                    
                                #     h["smile_graph_image"]=k
                                cp_data.append(h)

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
minRF = float(sys.argv[1])
minCS =sys.argv[2]
maxOR = float(sys.argv[3])
datasetname = sys.argv[4]
with open('result.txt','rb') as fp:
    fsubgraphs_data=pickle.load(fp)
    #print(type(fsubgraphs_data))
#print(fsubgraphs_data)
inpfile = "./2_"+str(minRF)+"_"+datasetname+"Flat_tra.txt"
#datasetname = "graph5_test2.txt"
#filepath = "C:\\Users\\user\\Documents\\CoverageGraph\\dataset1\\GraphData\\"
outfile = "Results.txt"
'''print(outfile)'''
obj=cmine(minRF, minCS, maxOR, inpfile, outfile,datasetname,fsubgraphs_data)
#print(fsubgraphs_data)
t3=time.time()
print( "data read",str(t3-t1))
candidate_patterns,SCPs,max_graphs = obj.expand()
t2 = time.time()
#print("cjjfjfjf")
#print(cp_data)
#print(len(cp_data[-1]["smile_graph_image"]))
with open("cp_result.txt","wb") as f:
    pickle.dump(cp_data,f)
#print(cp_data)
f.close()
#print( "process done",str(t2-t1))
f = open("./Results.txt",'w')
f.write("execution_time:"+str(t2-t1)+'\n')
f.write("Number of Candidate Patterns:"+str(candidate_patterns)+'\n')
f.write("Number of Coverage Patterns:"+str(SCPs)+'\n')
f.write("Max-graphs:"+str(max_graphs)+'\n')
f.close()
