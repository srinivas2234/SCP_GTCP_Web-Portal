from genericpath import isdir
from flask import Flask, render_template, request, redirect, send_from_directory
from flask import jsonify
from flask.wrappers import Response
from flask_cors import CORS,cross_origin
from flask import request
from werkzeug.utils import secure_filename
import networkx as nx
import matplotlib.pyplot as plt
import sys
import json
import time
import os
from os import listdir
from os.path import isfile, join
import base64
import graph
import subprocess
import pickle
from graph import AUTO_EDGE_ID
from graph import Graph
from graph import VACANT_GRAPH_ID
from graph import VACANT_VERTEX_LABEL
from graph import min_no_vertices
import collections
import shutil

app = Flask(__name__)
# app.config['SECRET_KEY'] = 'the quick brown fox jumps over the lazy   dog'
# app.config['CORS_HEADERS'] = 'Content-Type'
cors = CORS(app, resources={r"/*": {"origins": "http://localhost:3000"}})
app.config['CORS_HEADERS']='Content-Type'
# , resources={r"/*": {"origins": "*"}}
HOST = '0.0.0.0'
PORT = 5000

@app.route('/', methods = ['GET','POST'])
@cross_origin()

def check():
    print("check")
    if(request.form['data']=='1'):
        if(request.form['submit']!='1'):
            reading=request.form['inputdata']
            reading1=request.form['inputdata']
            reading=reading.split('\n')
            #print(reading)
            plt.clf()
            g = nx.Graph()
            vertexset=set()
            error_messages=''
            for line in reading:
                #print(vertexset)
                #print(line)
                line=line.strip().split()
                if(len(line)!=0):
                    if line[0]=='v':
                        #print("went in")
                        vertexset.add(line[1])
                    if line[0]=='e':
                        #print("wentinto edge")
                        if line[1] not in vertexset or line[2] not in vertexset:
                            #print("error line1")
                            if(line[1] not in vertexset):
                                error_messages=error_messages+line[1]+" is not a vertex"+" "
                            if(line[2] not in vertexset):
                                error_messages=error_messages+line[2]+" is not a vertex"+" "
                        else:
                            g.add_edge(int(line[1]),int(line[2]))
            return_img=''
            #print("erro meessag")
            #print(error_messages)
            nx.draw_planar(g, with_labels = True)
            plt.savefig("filename.png")
            with open("filename.png","rb") as img:
                k=str(base64.b64encode(img.read()))
                k=k[2:]
                return_img=k[:-1]
            # os.remove("filename.png")
            #print(error_messages)
            if(request.form['dosave']=='1'):
                #print("went into savw")
                cnt=request.form['no_of_si']
                #print("cout",type(cnt))
                if(int(cnt)==1):
                    #print("went into if")
                    if(os.path.exists('./final.txt')):
                        os.remove('./final.txt')
                    f=open('final.txt','w')
                else:
                    f=open('final.txt','a')
                trscn='t # '+str(cnt)+' '
                f.write(trscn)
                f.write(reading1)
                f.write('\n')
                f.write('\n')
                f.close()
            res=jsonify({"return_img":return_img,"error_messages":error_messages}),200
            return res
        else:
            if(os.path.exists('./final.txt')):
                f=open('./final.txt','a')
                f.write('t # -1')
            res=jsonify({"written":True})
            return res
    elif(request.form['data']=='0'):
        #print("data.0")
        maxor=request.form['maxor']
        mincs=request.form['mincs']
        minrf=request.form['minrf']
        analysis_type=request.form['analysistype']
        data_set_name=request.form['selected_data']
        #print(request.form["structure_of_interest"])
        structure_of_interest=int(request.form["structure_of_interest"])
        if(structure_of_interest):
            structure_content=request.form["structure_content_file"]

            c=structure_content.split(",")
            #print("Structure")
            
            structure_file=open("./structure.data",'w')
            for i in c:
                structure_file.write(i+'\n')
            structure_file.write("t # -1"+'\n')
            structure_file.close()
        
        data_set_file=data_set_name+".txt"

        r=request.form['file_content'][23:]
        r=base64.b64decode(r)
        r.decode('ascii')
        f=open('./graphdata/'+request.form['selected_dataset'],'wb')
        f.write(r)
        f.close()
        
        data_set_file=request.form['selected_dataset']
        data_set_name=data_set_file.split(".")[0]
        #print("structure")
        #print(structure_of_interest)
        
        if(analysis_type=='1'):
            start=time.time()

            if(structure_of_interest==0):
                print("no structure")
                s='python -m gspan_mining -s '+str(minrf)+" ./graphdata/"+data_set_file
            else:
                print("hi")
                s='python -m gspan_mining_struct -s '+str(minrf)+" "+"./graphdata/"+data_set_file+" "+'./structure.data'


            os.system(s)
            with open('result.txt','rb') as fp:
                d=pickle.load(fp)
                #print(type(d))
            for i in range(len(d)):
                #print(d[i])
                d[i]["vertices"]=len(d[i]["graph"]["nodes"])
                d[i]["edges"]=len(d[i]["graph"]["links"])
            
            img_folder = data_set_name+"data_2"+"_"+str(minrf)
            img_folder_path = str(img_folder+"/")
            img_file_names = [f for f in listdir(img_folder_path) if isfile(join(img_folder_path,f))]
            f=open("2_"+str(minrf)+"_"+data_set_name+"_output.txt",'r')
            fsg=[]
            supports=[]
            img=1
            edges=0
            vertices=0
            for i in f:
                i=i.strip('\n')
                j=i.split(" ")
                if(j[0]=='t'):
                    img=int(j[2])
                    temp={}
                    temp["image_name"]=img
                    temp["vertices"]=0
                    temp["edges"]=0
                    temp["nodes"]=[]
                    temp["links"]=[]
                    temp["image_src"]=""

                if(j[0] == "Support"):
                    temp["support"]=j[1]
                    if(j[1] not in supports):
                        supports.append(j[1])
                if(j[0] == "x"):
                    c=j[1].split(",")
                    temp["where"]=[]
                    for k in c:
                        temp["where"].append(k)
                    fsg.append(temp)
                    if(temp["vertices"]>vertices):
                        vertices=temp["vertices"]
                    if(temp["edges"]>edges):
                        edges=temp["edges"]
                if(j[0]=="v"):
                    temp["vertices"]=temp["vertices"]+1
                    temp["nodes"].append({"id":j[1],"label":j[2]})
                if(j[0]=="e"):
                    temp["edges"]=temp["edges"]+1
                    temp["links"].append({"source":j[1],"target":j[2],"label":j[3]})
                
            f = open("gSpan_FSM_"+data_set_name+"_stats.txt")
            li1 = []
            for i in f:
                i=i.strip('\n')
                i=i.split(":")
                li1.append(i[1])


            
            
            if(structure_of_interest==0):
                s='python cmine.py '+str(minrf)+" "+str(mincs)+" "+str(maxor)+" "+data_set_name+" "
            else:
                s='python cmine_struct1.py '+str(minrf)+" "+str(mincs)+" "+str(maxor)+" "+data_set_name
            end=time.time()
            os.system(s)
            no_of_coverage=0
            img_folder = "./scp_images_"+data_set_name+"_"+mincs
            img_folder_path = str(img_folder+"/")
            img_folder_names = [f for f in listdir(img_folder_path) if (os.path.isdir(join(img_folder_path,f)) and len(os.listdir(join(img_folder_path,f)))!=0)]
            img_folder_names=sorted(img_folder_names)

            smiles_graphs_img_folder="./scp_smiles_graphs_"+data_set_name+"_"+mincs
            smiles_graphs_img_folder_path=str(smiles_graphs_img_folder+"/")
            smiles_graphs_folder_names=[f for f in listdir(smiles_graphs_img_folder_path) if (os.path.isdir(join(smiles_graphs_img_folder_path,f)) and len(os.listdir(join(smiles_graphs_img_folder_path,f)))!=0)]
            smiles_graphs_folder_names=sorted(smiles_graphs_folder_names)

            coverages=0
            if(len(img_folder_names)!=0):
                coverages=img_folder_names[-1].split("_")[-1]
            coverage_patterns=[]
            f=open("./"+data_set_name+"_"+str(mincs)+"_"+str(minrf)+"_"+str(maxor))
            temp={}
            for i in f:
                i=i.strip('\n')
                j=i.split(" ")
                if(j[0]=="Coverage"):
                    if(temp!={}):
                        coverage_patterns.append(temp)
                    temp={}
                    temp["coverage"]=j[1]
                    temp["image_info"]=[]
                    temp["smile_graphs_image_info"]=[]
                elif(j[0]=='i'):
                    t={}
                    t['pattern_id']=j[1]
                elif(j[0]=='cs'):
                    t['cs']=j[1]
                elif(j[0]=='or'):
                    t['or']=j[1]
                    temp["image_info"].append(t)
            temp=coverage_patterns
            #print(temp)
            coverage_patterns=[]
            for i in temp:
                if(i['image_info']==[]):
                    pass
                else:
                    coverage_patterns.append(i)
            for i in img_folder_names:
                '''no_of_coverage=no_of_coverage+1'''
                #print(i)
                img_file_names=[f for f in listdir(img_folder_path+i) if isfile(join(img_folder_path+i+"/",f))]
                '''print("image_file_names")
                temp={}
                coverage=i.split("_")
                temp["coverage"]=coverage[1]
                
                temp["image_src"]=[]'''
                count=0
                for j in img_file_names: 

                    with open(str(img_folder_path+i+"/"+j),"rb") as img:
                        k = str(base64.b64encode(img.read()))
                        k = k[2:]
                        #print(count)
                        #print(no_of_coverage)
                        coverage_patterns[no_of_coverage]["image_info"][count]["image_src"]=k[:-1]
            
                    count=count+1
                scp_smile_graphs_image_file_names=[f for f in listdir(smiles_graphs_img_folder_path+i) if isfile(join(smiles_graphs_img_folder_path+i+"/",f))]
                for j in scp_smile_graphs_image_file_names:
                    with open(str(smiles_graphs_img_folder_path+i+"/"+j),"rb") as img:
                        k = str(base64.b64encode(img.read()))
                        k = k[2:]
                        #print(count)
                        #print(no_of_coverage)
                        coverage_patterns[no_of_coverage]["smile_graphs_image_info"].append(k[:-1])
                no_of_coverage=no_of_coverage+1
            with open('cp_result.txt','rb') as fcp:
                cp_data=pickle.load(fcp)
                #print(type(cp_data))
            f = open("./Results.txt")
            #print("coverage_patterns")
            #print(coverage_patterns)
            li2 = []
            for i in f:
                i=i.strip('\n')
                i=i.split(":")
                li2.append(i[1])
            #print(type(cp_data))
            #print(len(cp_data))
            cp1_data=cp_data
            sorted_by_coverage_support=cp_data.sort(key=lambda x:[x["cs"],x["or"],x["size"]])
            sorted_by_overlap_ratio=cp_data.sort(key=lambda x:[x["or"],x["cs"],x["size"]])
            sorted_by_size=cp_data.sort(key=lambda x:[x["size"],x["cs"],x["or"]])
            response= jsonify({"cp_data":cp1_data,"filter_cs":sorted_by_coverage_support,"filter_or":sorted_by_overlap_ratio,"filter_size":sorted_by_size,"data":d,"fsubgraphs":li1[1],"fexecution_time":li1[3],"avgtransactions":round(float(li1[3]),2),"image_info":fsg,"atype":1,"supports":supports,"vertices":vertices,"edges":edges,"coverage_patterns":coverage_patterns,"no_of_coverages":li2[-1],"number_of_candidate_patterns":li2[1],"number_of_scps":li2[2],"execution_time":str(round(float(li1[3])+float(li2[0]),2))}),200
            return response
        
        elif(analysis_type=='2'):
            #print("came to transac")
            os.system("python data_for_images.py ./graphdata/"+data_set_name)
            start=time.time()
            if(structure_of_interest==0):
                s='python -m gspan_mining -s '+str(minrf)+" ./graphdata/"+data_set_file
            else:
                #print("hi")
                s='python -m gspan_mining_struct -s '+str(minrf)+" ./graphdata/"+data_set_file+" "+"./structure.data"


            os.system(s)
            with open('result.txt','rb') as fp:
                d=pickle.load(fp)
                #print(type(d))
            for i in range(len(d)):
                #print(d[i])
                d[i]["vertices"]=len(d[i]["graph"]["nodes"])
                d[i]["edges"]=len(d[i]["graph"]["links"])
            
            #img_folder = data_set_name+"data_2"+"_"+str(minrf)
            #img_folder_path = str(img_folder+"/")
            #img_file_names = [f for f in listdir(img_folder_path) if isfile(join(img_folder_path,f))]
            f=open("2_"+str(minrf)+"_"+data_set_name+"_output.txt",'r')
            fsg=[]
            supports=[]
            img=1
            edges=0
            vertices=0
            for i in f:
                i=i.strip('\n')
                j=i.split(" ")
                if(j[0]=='t'):
                    img=int(j[2])
                    temp={}
                    temp["image_name"]=img
                    temp["vertices"]=0
                    temp["edges"]=0
                    temp["nodes"]=[]
                    temp["links"]=[]
                    #with open(str(img_folder_path+str(img)+".png"),"rb") as img:
                        #k = str(base64.b64encode(img.read()))
                        #k = k[2:]
                        #temp["image_src"]=k[:-1]

                if(j[0] == "Support"):
                    temp["support"]=j[1]
                    if(j[1] not in supports):
                        supports.append(j[1])
                if(j[0] == "x"):
                    c=j[1].split(",")
                    temp["where"]=[]
                    for k in c:
                        temp["where"].append(k)
                    fsg.append(temp)
                    if(temp["vertices"]>vertices):
                        vertices=temp["vertices"]
                    if(temp["edges"]>edges):
                        edges=temp["edges"]
                if(j[0]=="v"):
                    temp["vertices"]=temp["vertices"]+1
                    temp["nodes"].append({"id":j[1],"label":j[2]})
                if(j[0]=="e"):
                    temp["edges"]=temp["edges"]+1
                    temp["links"].append({"source":j[1],"target":j[2],"label":j[3]})
                
            f = open("gSpan_FSM_"+data_set_name+"_stats.txt")
            li1 = []
            for i in f:
                i=i.strip('\n')
                i=i.split(":")
                li1.append(i[1])


            
            
            if(structure_of_interest==0):
                s='python3 SetCoverProblem_gSpan.py '+str(minrf)+" "+str(mincs)+" "+str(maxor)+" "+data_set_name+" "+str(1)+" "+str(0)
            else:
                s='python3 SetCoverProblem_gSpan.py '+str(minrf)+" "+str(mincs)+" "+str(maxor)+" "+data_set_name+" "+str(1)+" "+str(1)
 
            os.system(s)
            end=time.time()
            no_of_coverage=0
            #img_folder = "./gtcp_images_"+data_set_name+"_"+mincs
            #img_folder_path = str(img_folder+"/")
            #img_folder_names = [f for f in listdir(img_folder_path) if (os.path.isdir(join(img_folder_path,f)) and len(os.listdir(join(img_folder_path,f)))!=0)]
            #img_folder_names=sorted(img_folder_names)
            #print(img_folder_names)
            #smiles_graphs_img_folder="./gtcp_smiles_graphs_"+data_set_name+"_"+mincs
            #smiles_graphs_img_folder_path=str(smiles_graphs_img_folder+"/")
            #smiles_graphs_folder_names=[f for f in listdir(smiles_graphs_img_folder_path) if (os.path.isdir(join(smiles_graphs_img_folder_path,f)) and len(os.listdir(join(smiles_graphs_img_folder_path,f)))!=0)]
            #smiles_graphs_folder_names=sorted(smiles_graphs_folder_names)

            #coverages=0
            #if(len(img_folder_names)!=0):
                #coverages=img_folder_names[-1].split("_")[-1]
            #coverage_patterns=[]
            '''f=open("./gtcp_"+data_set_name+"_"+str(mincs)+"_"+str(minrf)+"_"+str(maxor))
            temp={}
            for i in f:
                i=i.strip('\n')
                j=i.split(" ")
                if(j[0]=="Coverage"):
                    if(temp!={}):
                        coverage_patterns.append(temp)
                    temp={}
                    temp["coverage"]=j[1]
                    temp["image_info"]=[]
                    temp["smile_graphs_image_info"]=[]
                elif(j[0]=='i'):
                    t={}
                    t['pattern_id']=j[1]
                elif(j[0]=='cs'):
                    t['cs']=j[1]
                elif(j[0]=='or'):
                    t['or']=j[1]
                    temp["image_info"].append(t)
            temp=coverage_patterns
            #print(temp)
            coverage_patterns=[]
            for i in temp:
                if(i['image_info']==[]):
                    pass
                else:
                    coverage_patterns.append(i)
            for i in img_folder_names:
                no_of_coverage=no_of_coverage+1
                #print(i)
                img_file_names=[f for f in listdir(img_folder_path+i) if isfile(join(img_folder_path+i+"/",f))]
                #print("image_file_names")
                temp={}
                coverage=i.split("_")
                temp["coverage"]=coverage[1]
                
                #temp["image_src"]=[]
                count=0
                for j in img_file_names: 

                    with open(str(img_folder_path+i+"/"+j),"rb") as img:
                        k = str(base64.b64encode(img.read()))
                        k = k[2:]
                        #print(count)
                        #print(no_of_coverage)
                        coverage_patterns[no_of_coverage]["image_info"][count]["image_src"]=k[:-1]
            
                    count=count+1
                scp_smile_graphs_image_file_names=[f for f in listdir(smiles_graphs_img_folder_path+i) if isfile(join(smiles_graphs_img_folder_path+i+"/",f))]
                for j in scp_smile_graphs_image_file_names:
                    with open(str(smiles_graphs_img_folder_path+i+"/"+j),"rb") as img:
                        k = str(base64.b64encode(img.read()))
                        k = k[2:]
                        #print(count)
                        #print(no_of_coverage)
                        coverage_patterns[no_of_coverage]["smile_graphs_image_info"].append(k[:-1])
                no_of_coverage=no_of_coverage+1'''
            with open('gtcp_result.txt','rb') as fcp:
                cp_data=pickle.load(fcp)
                #print(type(cp_data))
            f = open("./Results.txt")
            #print("coverage_patterns")
            #print(coverage_patterns)
            li2 = []
            for i in f:
                i=i.strip('\n')
                i=i.split(":")
                li2.append(i[1])
            #print(len(cp_data),len(cp_data[0])) 
            #print(cp_data)
            #print(len(cp_data))
            cp_data.sort(key=lambda x:[x["cs"],x["or"]])
            cp_data=cp_data[:21]
            #print(type(cp_data))
            #print(len(cp_data))
            #print(li2[1])
            #print("no of coverages")
            print(li2)
            print(round(float(li1[3])+float(li2[0]),2))
            response= jsonify({"cp_data":cp_data,"data":d,"fexecution_time":li1[3],"fsubgraphs":li1[1],"avgtransactions":round(float(li1[3]),2),"image_info":fsg,"atype":1,"supports":supports,"vertices":vertices,"edges":edges,"coverage_patterns":[],"no_of_coverages":li2[-1],"number_of_candidate_patterns":li2[1],"number_of_scps":li2[2],"execution_time":str(round(float(li1[3])+float(li2[0]),2))}),200
            return response



            


            

        
if __name__ == "__main__":
	app.run(host=HOST, port=PORT, debug=True)