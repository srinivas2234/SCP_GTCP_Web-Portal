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
from graph import AUTO_EDGE_ID, Vertex
from graph import Graph
from graph import VACANT_GRAPH_ID
from graph import VACANT_VERTEX_LABEL
from graph import min_no_vertices
import subprocess
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
    print("CHECK")
    print(request.form['data'])
    if(request.form['data']=='1'):
        print("cameto structure")
        if(request.form['submit']!='1'):
            reading=request.form['inputdata']
            reading1=request.form['inputdata']
            reading=reading.split('\n')
            print(reading)
            plt.clf()
            g = nx.Graph()
            vertexset=set()
            error_messages=''
            for line in reading:
                print(vertexset)
                print(line)
                line=line.strip().split()
                if(len(line)!=0):
                    if line[0]=='v':
                        print("went in")
                        vertexset.add(line[1])
                    if line[0]=='e':
                        print("wentinto edge")
                        if line[1] not in vertexset or line[2] not in vertexset:
                            print("error line1")
                            if(line[1] not in vertexset):
                                error_messages=error_messages+line[1]+" is not a vertex"+" "
                            if(line[2] not in vertexset):
                                error_messages=error_messages+line[2]+" is not a vertex"+" "
                        else:
                            g.add_edge(int(line[1]),int(line[2]))
            return_img=''
            print("erro meessag")
            print(error_messages)
            nx.draw_planar(g, with_labels = True)
            plt.savefig("filename.png")
            with open("filename.png","rb") as img:
                k=str(base64.b64encode(img.read()))
                k=k[2:]
                return_img=k[:-1]
            # os.remove("filename.png")
            print(error_messages)
            if(request.form['dosave']=='1'):
                print("went into savw")
                cnt=request.form['no_of_si']
                print("cout",type(cnt))
                if(int(cnt)==1):
                    print("went into if")
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
        maxor=request.form['maxor']
        mincs=request.form['mincs']
        minrf=request.form['minrf']
        analysis_type=request.form['analysistype']
        print("i")
        print("j")
        data_set_name=request.form['selected_data']
        print(1)
        data_set_file=data_set_name+".txt"

        r=request.form['file_content'][23:]
        r=base64.b64decode(r)
        r.decode('ascii')
        f=open('./graphdata/'+request.form['selected_dataset'],'wb')
        f.write(r)
        f.close()
        data_set_file=request.form['selected_dataset']
        data_set_name=data_set_file.split(".")[0]
       
        print(2)
        print("analysis",analysis_type)
        if(analysis_type=='1'):


            s='python -m gspan_mining -s '+str(mincs)+" ./graphdata/"+data_set_file

            os.system(s)

            '''img_folder = 'papertoydatadata_2'+"_"+str(mincs)
            img_folder_path = str(img_folder+"/")
            img_file_names = [f for f in listdir(img_folder_path) if isfile(join(img_folder_path,f))]

            f=open("2_"+str(mincs)+"_"+"papertoydata_output.txt",'r')
            frequencies_of_images={}
            where_of_images={}
            file=0
            for i in f:
                i=i.strip('\n')
                j=i.split(" ")
                if(j[0]=='t'):
                    file=j[2]
                if(j[0] == "Support"):
                    frequencies_of_images[file]=j[1]
                if(j[0] == "x"):
                    where_of_images[file]=j[1]
            print(frequencies_of_images)
            d={}
            for i in img_file_names:
                i=i.split('.')
                d[i[0]] = frequencies_of_images[i[0]]
            print(img_file_names)
            print(where_of_images)
            images={"img_base64":[],"img_frequencies":[],"where":[]}
            for i in img_file_names:
                with open(str(img_folder_path+i),"rb") as img:
                    k = str(base64.b64encode(img.read()))
                    k = k[2:]
                    images["img_base64"].append(k[:-1])
                    i=i.strip('\n')
                    i=i.split(".")
                    images["img_frequencies"].append(d[i[0]])
                    images["where"].append(where_of_images[i[0]])

            f = open("gSpan_FSM_papertoydata_stats.txt")
            li = []
            for i in f:
                i=i.strip('\n')
                i=i.split(":")
                li.append(i[1])'''
            
            img_folder = data_set_name+"data_2"+"_"+str(mincs)
            img_folder_path = str(img_folder+"/")
            img_file_names = [f for f in listdir(img_folder_path) if isfile(join(img_folder_path,f))]
            f=open("2_"+str(mincs)+"_"+data_set_name+"_output.txt",'r')
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
                    with open(str(img_folder_path+str(img)+".png"),"rb") as img:
                        k = str(base64.b64encode(img.read()))
                        k = k[2:]
                        temp["image_src"]=k[:-1]

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
                if(j[0]=="e"):
                    temp["edges"]=temp["edges"]+1
                

            f = open("gSpan_FSM_"+data_set_name+"_stats.txt")
            li1 = []
            for i in f:
                i=i.strip('\n')
                i=i.split(":")
                li1.append(i[1])


            
            
        
            s='python cmine.py '+str(minrf)+" "+str(mincs)+" "+str(maxor)+" "+data_set_name
            os.system(s)
            print("hiiiiii")
            no_of_coverage=0
            img_folder = "./scp_images_"+data_set_name+"_"+mincs
            img_folder_path = str(img_folder+"/")
            img_folder_names = [f for f in listdir(img_folder_path) if (os.path.isdir(join(img_folder_path,f)) and len(os.listdir(join(img_folder_path,f)))!=0)]
            print("imagefolder")
            img_folder_names=sorted(img_folder_names)
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
                elif(j[0]=='i'):
                    t={}
                    t['pattern_id']=j[1]
                elif(j[0]=='cs'):
                    t['cs']=j[1]
                elif(j[0]=='or'):
                    t['or']=j[1]
                    temp["image_info"].append(t)

            for i in img_folder_names:
                '''no_of_coverage=no_of_coverage+1'''
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
                        coverage_patterns[no_of_coverage]["image_info"][count]["image_src"]=k[:-1]
                    count=count+1
                no_of_coverage=no_of_coverage+1
                
            f = open("./Results.txt")
            li2 = []
            for i in f:
                i=i.strip('\n')
                i=i.split(":")
                li2.append(i[1])
                
            response= jsonify({"fsubgraphs":li1[1],"avgtransactions":round(float(li1[3]),2),"image_info":fsg,"atype":1,"supports":supports,"vertices":vertices,"edges":edges,"coverage_patterns":coverage_patterns,"no_of_coverages":no_of_coverage,"number_of_candidate_patterns":li2[1],"number_of_scps":li2[2],"execution_time":str(round(float(li2[0]),2))}),200
            return response


            


            

        
if __name__ == "__main__":
	app.run(host=HOST, port=PORT, debug=True)