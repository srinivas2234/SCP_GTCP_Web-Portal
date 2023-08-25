
''' This code generates Set cover problem patterns and overlap ratio is calculated 
w.r.t transactions as items. Prune transactions based on minRF'''

import time
import sys
import os
import pickle
import numpy as np
import copy
import itertools
from bitarray import bitarray
import shutil
from PIL import Image
import base64
import cv2
tot_CPs=[]
freqList=[]
#freqEle=[]
global gtcp_data
class cmine():
    def __init__(self, minRF, minCS, maxOR, inpfile, outfile,datasetname,fsubgraphs_data):
        global gtcp_data
        self.minRF = minRF
        self.minCS = minCS
        self.maxOR = maxOR
        self.inpfile = inpfile
        self.outfile = outfile
        self.nofs = self.getlines(inpfile)
        #self.fout = open(outfile,'w')
        self.noofCTP=0
        self.noof_Candi_CTP=0
        self.NOk = []
        self.Candi_CTP=[]
        self.final_CTPs=[]
        self.nofs,self.bitpattern = self.dbscanSCP(inpfile)
        #print("pattern")
        #print(self.nofs,self.bitpattern)

        self.source_dir_graph_smiles="./"+"gtcp_smiles_images_"+datasetname
        self.dest_dir_graph_smiles="./gtcp_coverage_smiles_"+datasetname+"_"+str(self.minCS)
        if(os.path.exists(self.dest_dir_graph_smiles)):
            shutil.rmtree(self.dest_dir_graph_smiles)
        os.mkdir(self.dest_dir_graph_smiles)
        #[self.items, self.bitpattern,self.TidKey_Dict] = self.dbscan(inpfile)
        self.total=0
        "Checking feature coverage of a graph greaterthan minFC or not"
        #print("No.of 1 size candidates:",len(self.bitpattern))  
        if(os.path.exists(self.dest_dir_graph_smiles+"/Coverage_1")):
            shutil.rmtree(self.dest_dir_graph_smiles+"/Coverage_1")
        os.mkdir(self.dest_dir_graph_smiles+"/Coverage_1")  
        #print("hi")  
        gtcp_data=[]
        #print(len(fsubgraphs_data))
        for key, value in self.bitpattern.items():
            self.noof_Candi_CTP=self.noof_Candi_CTP+1
            if (1.0*value.count()/self.nofs) >= minCS:
                #print(key)
                tt=[[key],1.0*value.count()/self.nofs]
                #print(tt)
                tt.append(0)
                #print(fsubgraphs_data[key])
                
                h={"id":str(self.total),"or":str(self.maxOR),"cs":str(round(float(1.0*value.count()/self.nofs),2)),"size":str(1),"graphs":[],"string":''}
                h["graphs"].append(fsubgraphs_data[int(key)])
                #print(fsubgraphs_data[key])
                #print(key)
                self.total=self.total+1
                h["string"]=h["string"]+fsubgraphs_data[int(key)]["string"]
                with open(self.source_dir_graph_smiles+"/"+str(key)+".svg","rb") as img:
                    k = str(base64.b64encode(img.read()))
                    k = k[2:]
                    #print(count)
                    #print(no_of_coverage)
                    h["smile_graph_image"]=[k[:-1]]
                
                
                tot_CPs.append(tt)
                shutil.copy(self.source_dir_graph_smiles+"/"+str(key)+".svg",self.dest_dir_graph_smiles+"/Coverage_1")
                gtcp_data.append(h)
            elif (1.0*value.count()/self.nofs) >= minRF:
                temp=[1.0*value.count()/self.nofs,key]
                self.NOk.append(temp)
        #print(k)
        c=base64.b64decode(k)
        image=open("new1.png",'wb')
        image.write(c)
        '''if(os.path.exists(self.dest_dir)):
            shutil.rmtree(self.dest_dir)
        os.mkdir(self.dest_dir)       '''
        print("NO.of 1-size candidate patterns, 1-size GCP=",self.noof_Candi_CTP,len(tot_CPs))
        #self.NOk.sort(reverse=True)
        sorteditems = sorted(self.NOk, key = lambda a: (-a[0],a[1]))
        #print(len(sorteditems))        
        self.NOk=[]
        for i in sorteditems:
            self.NOk.append([i[1]])
            
            
           
        print("1-size nonOverlap transactions",len(self.NOk))
        freqItemcnt=[]
        for key, value in self.bitpattern.items():
            #print(key,value)
            freqItemcnt.append(value.count())
        #for ii in freqItemcnt:
         #   print(ii)    
        one_size_coverage=[]       
        
      
        
    def get_overlapratio_cs(self, new_tp):  
        cov_set=self.nofs*bitarray('0')
        last_tra_cov=self.nofs*bitarray('0')
        for tid in new_tp[:-1]:
            cov_set = cov_set | self.bitpattern[tid]
        #for tid in new_tp[-1]:
        last_tra_cov = self.bitpattern[new_tp[-1]]
            
        cov_sup = 1.0*(cov_set | last_tra_cov).count()/self.nofs
        
        ov_ratio= 1.0*(cov_set & last_tra_cov).count()/(last_tra_cov.count())
        return cov_sup,ov_ratio
   


    def expand(self):
        global gtcp_data
        cnt = 0
        cnt1 = 0
        length = 2
        coverage=2
        #print("went")
        while len(self.NOk)>0:
            if(os.path.exists(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage))):
                shutil.rmtree(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage))
            #print("length",length,len(self.NOk))
            l_size_GCP=0
            l_size_Candi_pat=0
            temp_NOk = self.NOk
            self.NOk = []
            #print(temp_NOk)
            count=1
            for i in range(len(temp_NOk)):
                for j in range(i+1, len(temp_NOk)):
                    #print("went_2")
                    cnt += 1
                    if temp_NOk[i][:-1] == temp_NOk[j][:-1]:
                        cnt1 += 1
                        newpattern = temp_NOk[i] + [temp_NOk[j][-1]]
                        cs,ov_ra=self.get_overlapratio_cs(newpattern)
                        l_size_Candi_pat=l_size_Candi_pat+1
                        self.noof_Candi_CTP=self.noof_Candi_CTP+1
                        #print(ov_ra,maxOR,type(ov_ra),type(maxOR))
                        #print(cs,minCS,type(cs),type(minCS))
                        if ov_ra <= maxOR:                            
                            if cs >= minCS:
                                #print("yes")
                                if(not os.path.exists(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage))):
                                    os.mkdir(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage))
                                opened_images=[]
                                opened_smiles_images=[]

                                h={"id":str(self.total),"or":str(ov_ra),"cs":str(round(float(cs),2)),"size":len(newpattern),"graphs":[],"string":''}
                                self.total=self.total+1
                                for k in newpattern:
                                    '''print("true")'''
                                    '''print(k)'''
                                    h["graphs"].append(fsubgraphs_data[int(k)-1])
                                    h["string"]=h["string"]+fsubgraphs_data[int(k)-1]["string"]
                                    #opened_images.append(cv2.imread(self.source_dir+"/"+str(k)+".png"))
                                    opened_smiles_images.append(cv2.imread(self.source_dir_graph_smiles+"/"+str(k)+".svg"))
                                #h_min = min(img.shape[0] for img in opened_images)
                                #h_min_smiles_images=min(img.shape[0] for img in opened_smiles_images)
                                '''widths,heights=zip(*(k.size for k in opened_images))
                                total_width=sum(widths)
                                max_height=max(heights)
                                new_img = np.zeros(shape=(max_height, total_width, 3))
                                new_image=Image.new('RGB',(total_width,max_height))
                                x_offset=0
                                for k in opened_images:
                                    new_image.paste(k,(x_offset,0))
                                    x_offset=x_offset+k.size[0]'''
                                #im_list_resize = [cv2.resize(img,(int(img.shape[1] * h_min / img.shape[0]),h_min), interpolation=cv2.INTER_CUBIC) for img in opened_images]
                                #im_list_resize_smiles=[cv2.resize(img,(int(img.shape[1] * h_min_smiles_images / img.shape[0]),h_min_smiles_images), interpolation=cv2.INTER_CUBIC) for img in opened_smiles_images]
                                #new_image=cv2.hconcat(im_list_resize)
                                #new_image_smiles=cv2.hconcat(im_list_resize_smiles)
                                #cv2.imwrite(self.dest_dir+"/Coverage_"+str(coverage)+"/"+str(count)+".png",new_image)
                                #cv2.imwrite(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage)+"/"+str(count)+".svg",new_image_smiles)
                                
                                '''new_image.save(self.dest_dir+"/Coverage_"+str(coverage)+"/"+str(count)+".png")'''
                                #k = str(base64.b64encode(new_image_smiles.read()))
                                #k = k[2:]
                                #print(count)
                                #print(no_of_coverage)
                                #h["smile_graph_image"]="one"
                                h["smile_graph_image"]=[]
                                with open(self.dest_dir_graph_smiles+"/Coverage_"+str(coverage)+"/"+str(count)+".svg","rb") as img:
                                    k = str(base64.b64encode(img.read()))
                                    k = k[2:]
                                    if(count==0):
                                        print(k)
                                    #print(count)
                                    #print(no_of_coverage)
                                    h["smile_graph_image"].append(k[:-1])
                                gtcp_data.append(h)

                                count=count+1
                                self.noofCTP=self.noofCTP+1
                                temp=[newpattern,cs,ov_ra]
                                l_size_GCP=l_size_GCP+1
                                tot_CPs.append(temp)
                            else:
                                self.NOk.append(newpattern)
                            
                                
                    else:
                        break
            coverage=coverage+1
            length += 1
          
        return self.noof_Candi_CTP,self.noofCTP,length
    
   
    def getlines(self,filename):
        with open(filename,"r") as f:
            return sum(1 for _ in f)
            
    def dbscanSCP(self,db):
        # x=[]  
        noofbits = 0
        Tran_feat_Dict={}
        fidlist=[]
        g= open(db, 'r')
        for row in g:
            #if row[0]=='x':
            #print(row)
            row = row.rstrip('\n')
            r = row.split(' ')
            for fid in r:
                if fid not in fidlist and fid != '':
                    fidlist.append(fid)
        self.nofs=len(fidlist)
        #print("No.of features",len(fidlist))            
        f = open(db, 'r')
        Tidbitpat ={}        
        fid=0
        count=0
        for row in f:
            #if row[0]=='x':
            row = row.rstrip('\n')
            r = row.split(' ')
            bit_tra=self.nofs*bitarray('0')
            
            for fid in r:
                if fid !='':
                    fid=int(fid)
                    bit_tra[int(fid)-1] = 1
            Tran_feat_Dict[count]=bit_tra
                   
            count=count+1
        return self.nofs,Tran_feat_Dict
              
              

start_time = time.time()
minRF = float(sys.argv[1])
minCS = float(sys.argv[2])
maxOR = float(sys.argv[3])
datasetname = sys.argv[4]
writePatterns = sys.argv[5]
with open('./all_graphs.txt','rb') as fp:
    fsubgraphs_data=pickle.load(fp)
#print(fsubgraphs_data)
#print(len(fsubgraphs_data))
if(int(sys.argv[6])==0):
    inpfile = "./2_"+str(minRF)+"_"+datasetname+"Flat_tra.txt"
else:
    inpfile="./structure.data_"+str(minRF)+"_"+datasetname+"_results.txt"
outfile = "./GTCP"+str(datasetname)+"SetCover_Results.txt"
obj=cmine(minRF, minCS, maxOR, inpfile, outfile,datasetname,fsubgraphs_data)
candidate_patterns,CTPs,length = obj.expand()
CTP_time = time.time()
with open("gtcp_result.txt","wb") as f:
    pickle.dump(gtcp_data,f)
f.close()

print(str(datasetname)+", TC="+str(minRF)+", TPC="+str(minCS)+",overlap ratio"+str(maxOR)+", Exex Time="+str(CTP_time-start_time)+", No.of Candidate Patterns="+str(candidate_patterns)+",No.of GTCP="+str(len(tot_CPs)))
outStr=str(datasetname)+", TC="+str(minRF)+", TPC="+str(minCS)+",overlap ratio"+str(maxOR)+", Exex Time="+str(CTP_time-start_time)+", No.of Candidate Patterns="+str(candidate_patterns)+",No.of GTCP="+str(len(tot_CPs))
#outStr=str(datasetname)+","+str(minRF)+","+str(minCS)+","+str(maxOR)+","+str(len(tot_CPs))+","+str(candidate_patterns)+","+str(CTP_time-start_time)+"\n"
f=open("./Results.txt","w")
f.write("execution_time:"+str(CTP_time-start_time)+'\n')
f.write("Number of Candidate Patterns:"+str(candidate_patterns)+'\n')
f.write("Number of Coverage Patterns:"+str(len(tot_CPs))+'\n')
f.write("maz size of gtcp:"+str(length)+'\n')



print(tot_CPs)
with open("./"+str(datasetname)+"_GTCPs.txt", 'a', newline = '') as q:
    for pat in tot_CPs:
        q.write(str(pat)+"\n")    
    q.close()

