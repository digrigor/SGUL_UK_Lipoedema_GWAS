import os
import sys
import time
import glob
import subprocess
import re
import random
import string
import csv
import pandas as pd

from dependencies import *

def coordlist2vcf(coordlist, pgrange_df, col2get="file", chr="chr", start="start", end="end"):
    """Function that takes a list with genomic coordinates: ["chr5","3432423","3243342"]
    (if one of the positions miss you can just replace it with "") and returns a list
    of vcf files corresponding to this coordinates (read from a map text file)."""
    c1=coordlist[0]
    c2=coordlist[1]
    c3=coordlist[2]
    cindex=[]
    if c3=="":
        if c2!="":
            cindex+=pgrange_df.index[(pgrange_df[start]<int(c2)) & (pgrange_df[end]>int(c2)) & (pgrange_df[chr]==str(c1))].tolist()
            if len(cindex)==0:
                print("The coordinates you've entered are either not correct or there are not in the processed vcf files.")
                raise IndexError
            else:
                cfiles = pgrange_df.iloc[cindex][col2get]
            #cfiles=cfiles+pgrange_df[(pgrange_df[start]<int(c3)) & (pgrange_df[end]>int(c3)) & (pgrange_df[chr]==str(c1))][col2get].tolist()
        elif c2=="":
            cindex+=pgrange_df.index[(pgrange_df[chr]==str(c1))].tolist()
            if len(cindex)==0:
                print("The coordinates you've entered are either not correct or there are not in the processed vcf files.")
                raise IndexError
            else:
                cfiles = pgrange_df.iloc[cindex][col2get].tolist()
    elif c2!="" and c3!="":
        cfiles=[]
        cindex+=pgrange_df.index[(pgrange_df[start]<int(c2)) & (pgrange_df[end]>int(c2)) & (pgrange_df[chr]==str(c1))].tolist()
        cindex+=pgrange_df.index[(pgrange_df[start]<int(c3)) & (pgrange_df[end]>int(c3)) & (pgrange_df[chr]==str(c1))].tolist()
        try:
            cfiles+=[pgrange_df.iloc[cindex[0]][col2get]]
        except IndexError:
            print("The coordinates you've entered are not correct or there are not in the processed vcf files.")
        try:
            cfiles+=[pgrange_df.iloc[cindex[1]][col2get]]
        except IndexError:
            print("The coordinates you've entered are not correct or there are not in the processed vcf files.")
        try:
            cfiles+=pgrange_df.iloc[cindex[0]:cindex[0]][col2get].tolist()
        except IndexError:
            print("The coordinates you've entered are not correct or there are not in the processed vcf files.")
    cfiles = list(set(cfiles))
    return(cfiles)

def csvlines2list(csv_path):
    """Function that takes a path for a csv file, it opens it and returns each line as a
    python list element."""
    csv_in = csv_path.strip()
    with open(csv_in, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        outlist=[]
        for row in csv_reader:
            row = [x.strip() for x in row]
            row = [x.split("\xef\xbb\xbf")[1] if x.startswith("\xef\xbb\xbf") else x for x in row]
            outlist+=row
    return(outlist)


def getnum(text):
    retnnum = int(text.split("/")[-1].split("_")[5])
    return(retnnum)


def tidydirs(project, outdir, rundir):
    dirlist=[]
    dirlist.append(outdir+"/"+project) #project_out
    dirlist.append(rundir+"/"+project) #project_run
    dirlist.append(rundir+"/"+project+"/"+"jobscripts") #project_js
    dirlist.append(rundir+"/"+project+"/"+"tmp") #project_tmp
    dirlist.append(outdir+"/"+project+"/genomic_data") #genomic_data
    dirlist.append(outdir+"/"+project+"/interval_data") #interval_data
    for dir in dirlist:
        if not os.path.exists(dir):
            print(dir+" "+"has been created.")
            os.makedirs(dir)
    return(dirlist)


def randomString(stringLength=8):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def lsf_submit(command, jobprefix, to_wait_id="",
    wtime="24:00:00", nodes=1, cpu=1, mem="1gb", cwd="./", project_name="re_gecip_cardiovascular"):
    """Function which takes a pipeline command, runtime arguments, jobname, job dependencies
    as input, submits the job to the cluster and returns the job id"""
    if (type(to_wait_id) == str):
        to_wait_id = [to_wait_id]
    to_wait_id = [x for x in to_wait_id if x!=""]
    #print("debug:")
    #print(to_wait_id)
    #print(":".join(to_wait_id))
    if cwd.endswith("/") !=True:
        cwd = cwd + "/"
    out_filename = cwd + jobprefix + ".jobscript"
    to_wait_id_to_write = ["done("+x+")" for x in to_wait_id]
    whours = int(wtime.split(":")[0])
    mem = str(int(mem.split("gb")[0])*1000)
    if whours < 4: wqueue="short"
    elif whours >=4 and whours < 24: wqueue="medium"
    elif whours >= 24: wqueue="long"
    #print(out_filename)
    with open(out_filename, 'w+') as out_file:
        out_file.write('#!/bin/bash')
        out_file.write('\n')
        out_file.write('\n')
        out_file.write('#BSUB -q ' + wqueue)
        out_file.write('\n')
        out_file.write('#BSUB -P ' + project_name)
        out_file.write('\n')
        if len(to_wait_id_to_write)!=0:
            out_file.write('#BSUB -w "'+" && ".join(to_wait_id_to_write)+'"')
            out_file.write('\n')
        if whours >= 168:
            out_file.write('#BSUB -We 200:00')
            out_file.write('\n')
        #out_file.write('\n')
        out_file.write("#BSUB -n " + str(cpu))
        out_file.write('\n')
        #out_file.write("#BSUB -R rusage[mem=" + str(mem) + "]")
        out_file.write('\n')
        out_file.write("#BSUB -oo " + str(cwd + "/" + jobprefix + ".stdout"))
        out_file.write('\n')
        out_file.write("#BSUB -outdir " + cwd)
        out_file.write('\n')
        out_file.write("#BSUB -eo " + str(cwd + "/" + jobprefix + ".stderr"))
        out_file.write('\n')
        out_file.write("#BSUB -cwd " + str(cwd))
        out_file.write('\n')
        out_file.write('\n')
        out_file.write(command)
        out_file.write('\n')
        out_file.write('sleep 2')
        out_file.write('\n')
        out_file.write('exit')
    out_file.close()
    time.sleep(0.1)
    p = subprocess.Popen(["bsub"],stdout=subprocess.PIPE, stderr=subprocess.PIPE,stdin=open(out_filename, 'r'))
    out, err = p.communicate()
    jout = out
    jerr = err
    mout = re.findall(r"<(\d+)>",jout.decode("utf-8"))
    if mout[0]:
        mout = mout[0]
    else:
        print("No job identifier detected")
        raise()
    if "Job not submitted" in jerr.decode("utf-8"):
        print(jerr.decode("utf-8"))
        raise()
    return (mout)

def getnum(text):
    retnnum = int(text.split("/")[-1].split("_")[5])
    return(retnnum)
