import sys
import argparse
import subprocess
import os
import fnmatch
import sqlite3
import random
import numpy
import math
from scipy import stats

def db(args):

    conn = sqlite3.connect(args.prefix+".db")
    c = conn.cursor()
    A="CREATE TABLE SVDB (chr TEXT, start INT,end INT, id TEXT)"
    c.execute(A)

    input_tab=[]
    for line in open(args.tab):
        if line[0] == "#":
            continue
        content=line.strip().split()
        input_tab.append([ content[0].replace("chr",""), content[1] , content[2] ,content[3] ])
                
        if len(input_tab) > 1000000:
            c.executemany('INSERT INTO SVDB VALUES (?,?,?,?)',input_tab)          
            input_tab=[]
    if input_tab:
        c.executemany('INSERT INTO SVDB VALUES (?,?,?,?)',input_tab)
                    

    A="CREATE INDEX SNP ON SVDB (chr, start, end)"
    c.execute(A)
    conn.commit()
    conn.close()
    return()



def compute(args):    
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    
    for line in open(args.tab):
        if line[0] == "#" and line[1] == "#":
            print line.strip()
            continue

        content=line.strip().split("\t")
        content.append(str( distance(content[0].replace("chr",""),int(content[1]),c) ))
        content.append(str( distance(content[2].replace("chr",""),int(content[3]),c) ))

        content.append(str( alu_id(content[0].replace("chr",""),int(content[1]),c) ))
        content.append(str( alu_id(content[2].replace("chr",""),int(content[3]),c) ))

        print("\t".join(content))

    conn.close()
    return()

def distance(chr,pos,c):
    delta=10000
    i =0
    while True:
        i +=1
        A='SELECT start,end FROM SVDB WHERE chr == \'{}\' AND end > {} AND start < {} '.format(chr,int(pos)-delta,int(pos)+delta)
        d=[]
        for hit in c.execute(A):
            d.append( abs(pos- int( hit[0] )) )
            d.append( abs(pos- int( hit[1] )) )
            if pos >= int( hit[0]) and pos <= int( hit[1] ):
                d.append(0)
        if d:
            return( min(d) )
        delta = delta*10
        if i > 1000:
            return -1
def alu_id(chr,pos,c):
    delta=10000
    i =0
    while True:
        i +=1
        A='SELECT start,end,id FROM SVDB WHERE chr == \'{}\' AND end > {} AND start < {} '.format(chr,int(pos)-delta,int(pos)+delta)
        d=[]
        id_dict={}
        for hit in c.execute(A):

            d.append( abs(pos- int( hit[0] )) )
            id_dict[d[-1]]=hit[2]
            d.append( abs(pos- int( hit[1] )) )
            id_dict[d[-1]]=hit[2]
            if pos >= int( hit[0]) and pos <= int( hit[1] ):
                return(hit[2])
        if d:
            return( id_dict[min(d)] )

        delta = delta*10
        if i > 1000:
            return -1  

def bootstrap(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()

    reference,order=load_fasta(args.fa)
    contig_sizes={}
    for contig in reference:
        contig_sizes[contig.replace("chr","")]=len(reference[contig])
    del reference

    average_distance=[]
    distances=[]
    for line in open(args.tab):
        if line[0] == "#":
            continue
        content=line.strip().split("\t")

        if content[-1] == -1 or content[-2] == -1:
            continue

        distances.append([int(content[-4]),int(content[-3])])
        average_distance += [int(content[-4]),int(content[-3])]
    total_distance= sum(average_distance)
    n = 0
    f=open(args.out,"w")
    f.write("#{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("chromosomeA","posA","chromosomeB","posB","index","distA","distB"))
    p=0
    averages=[]
    for i in range(0,args.n):
        dist=0
        for j in range(0,len(distances)):
            chromosome=random.choice(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"])
            posA=random.randint(20000, contig_sizes[chromosome]-20000)
            posB=random.randint(20000, contig_sizes[chromosome]-20000)

            pos_a_dist=distance(chromosome,posB,c)
            pos_b_dist=distance(chromosome,posA,c)
            dist += pos_b_dist +  pos_a_dist
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chromosome,posA,chromosome,posB,i,pos_a_dist,pos_b_dist))
        averages.append(dist/len(distances))
        if dist <= total_distance:
            p+=1

    print("p_val,variant_avg_dist,simulated_avg_dist,simulated_sd")
    print("{},{},{},{}".format(p/float(args.n),numpy.average(distances),numpy.average(averages),numpy.std(averages) )) 
    conn.close()

def simulate(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()

    reference,order=load_fasta(args.fa)

    contig_sizes={}
    for contig in reference:
        contig_sizes[contig.replace("chr","")]=len(reference[contig])
    del reference

    f=open(args.out,"w")
    f.write("#{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("chromosomeA","posA","chromosomeB","posB","index","distA","distB"))
    chromosomes=numpy.random.choice(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],args.n)
    i=0
    for chromosome in chromosomes:
        posA=random.randint(20000, contig_sizes[chromosome]-20000)
        posB=random.randint(20000, contig_sizes[chromosome]-20000)

        A= distance(chromosome.replace("chr",""),posA,c)
        B= distance(chromosome.replace("chr",""),posB,c)
          
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chromosome,posA,chromosome,posB,i,A,B))
        i+=1
    f.close()
    distances=[]
    for line in open(args.out):
        if line[0] == "#":
            continue
        content=line.strip().split("\t")    
        distances += [int(content[-2]),int(content[-1])]

    print("n\taverage_dist\tstdev\t")
    print("{}\t{}\t{}".format(args.n,numpy.average(distances),numpy.std(distances))) 
    conn.close()

def homology(args):
    reference,order=load_fasta(args.fa)
    contig_sizes={}
    for contig in reference:
        contig_sizes[contig.replace("chr","")]=len(reference[contig])

    homology_len=[]
    for line in open(args.tab):
        if line[0] == "#":
            continue
        content=line.strip().split("\t")
        homology_len.append(int(content[0]))
    avg=[]
    p_list=[]
    for length in homology_len:
        dist=[]
        chromosomes=numpy.random.choice(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],args.n)
        p_var=0
        for chromosome in chromosomes:
            while True:
                posA=random.randint(20000, contig_sizes[chromosome]-20000)
                posB=random.randint(20000, contig_sizes[chromosome]-20000)
                d=compute_homology_len(reference[chromosome],posA,posB) 
                if d > -1:
                    dist.append(d)
                    break

            if dist[-1] >= length:
                p_var+=1
        p_list.append(p_var/float(args.n))

        avg.append(numpy.average(dist))
    
    chromosomes=numpy.random.choice(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],args.n*len(homology_len) )
    p=0
    avg_var=[]
    for n in range(0,args.n):
        dist=[]
        for length in homology_len:
            while True:
                posA=random.randint(20000, contig_sizes[chromosome]-20000)
                posB=random.randint(20000, contig_sizes[chromosome]-20000)
                d=compute_homology_len(reference[chromosome],posA,posB)
                if d > -1: 
                    dist.append(d)
                    break  
          
        if sum(dist) >= sum(homology_len):
            p += 1
        avg_var.append(numpy.average(dist))

    p=p/args.n
    print("#{},{},{}".format(float(p),numpy.average(avg_var),numpy.std(avg_var) ))
    for i in range(0,len(homology_len)):
        print "{}\t{}\t{}".format( homology_len[i],p_list[i],avg[i] )             
    
def compute_homology_len(ref_seq,posA,posB):
    hom_len = 0
    direction_A=random.randint(0, 1)
    direction_B=random.randint(0, 1)
    while True:
        base_A=ref_seq[posA]
        base_B=ref_seq[posB]
        if base_A == "N" or base_B == "N":
            hom_len = -1
            break

        elif base_A == base_B:
            hom_len += 1
            if direction_A:
                posA +=1
            else:
                posA += 1
            if direction_B:
                posB += 1 
            else:
                posB +=1
        else:
            break

    return hom_len

def load_fasta(fafile):
    sequence={}
    chromosome_order=[]
    #read the fast file 
    with open(fafile, 'r+') as f:
        reference = f.read()
    split_reference=reference.split(">")
    del reference
    del split_reference[0]
    #store the reference as a dictionary
    chromosome_len={}
    chromosomes=[]
    simulated_bases=0
    for chromosome in split_reference:
        content=chromosome.split("\n",1)
        sequence[content[0].split()[0]]=content[1].replace("\n","")
        chromosome_order.append(content[0].split()[0])
    del split_reference
    return(sequence,chromosome_order)



parser = argparse.ArgumentParser("""analyse genomic positions""")
parser.add_argument('--compute',action="store_true",help="find the distance between the points of an input file and the regions of a database")
parser.add_argument('--db',action="store_true",help="generate the databases from input bed files")
parser.add_argument('--bootstrap',action="store_true",help="compute p value of the distance to repeats")
parser.add_argument('--homology',action="store_true",help="compute p value of the homology of breakpoints")
parser.add_argument('--simulate',action="store_true",help="simulate a large number of points and compute the distances from objects in the db")
args, unknown = parser.parse_known_args()

if args.compute:
    parser = argparse.ArgumentParser("""compute the distance between the points of an input bed file and database""")
    parser.add_argument('--db',type=str,required=True,help="the input database")
    parser.add_argument('--tab',type=str,required=True,help="the input tab file")
    args, unknown = parser.parse_known_args()
    compute(args)
elif args.homology:
    parser = argparse.ArgumentParser("""compute the distance between the points of an input bed file and database""")
    parser.add_argument('--fa',type=str,required=True,help="the input reference file")
    parser.add_argument('--tab',type=str,required=True,help="the input tab file")
    parser.add_argument('--n',type=int,default = 1000,help="itterations, default = 100")
    args, unknown = parser.parse_known_args()
    homology(args)
elif args.bootstrap:
    parser = argparse.ArgumentParser("""compute the distance between the points of an input bed file and database""")
    parser.add_argument('--db',type=str,required=True,help="the input database")
    parser.add_argument('--tab',type=str,required=True,help="the input tab file")
    parser.add_argument('--out',type=str,default="dist.tab",help="output name of distance file")
    parser.add_argument('--fa',type=str,required=True,help="the input reference file")
    parser.add_argument('--n',type=int,default = 1000,help="itterations, default = 100")
    args, unknown = parser.parse_known_args()
    bootstrap(args)
elif args.simulate:
    parser = argparse.ArgumentParser("""compute the distance between the points of an input bed file and database""")
    parser.add_argument('--db',type=str,required=True,help="the input database")
    parser.add_argument('--out',type=str,default="simulate.tab",help="output name of distance file")
    parser.add_argument('--fa',type=str,required=True,help="the input reference file")
    parser.add_argument('--n',type=int,default = 1000,help="itterations, default = 1000")
    args, unknown = parser.parse_known_args()
    simulate(args)
elif args.db:
    parser = argparse.ArgumentParser("""generate a database of genomic regions""")
    parser.add_argument('--db',action="store_true",help="generate the databases from input bed files")
    parser.add_argument('--prefix',type=str,required=True,help="the prefix of the database")
    parser.add_argument('--tab',type=str,required=True,help="the input tab file")
    args, unknown = parser.parse_known_args()
    db(args)
else:
    print "type python distance_calculator.py --help"

