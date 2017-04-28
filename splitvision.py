import argparse
import FindTranslocations
import itertools
import os
import xlwt
import readVCF
import sqlite3
import gzip


def find_repeat(chr,pos,c):
    delta=10000
    i =0
    while True:
        i +=1
        A='SELECT start,end,id FROM SVDB WHERE chr == \'{}\' AND end > {} AND start < {} '.format(chr,int(pos)-delta,int(pos)+delta)
        d={}
        for hit in c.execute(A):
            d[ abs(pos- int( hit[0] )) ] = str(hit[2])
            d[ abs(pos- int( hit[1] )) ] = str(hit[2])
            if pos >= int( hit[0]) and pos <= int( hit[1] ):
                d[0]=str(hit[2])
        if d:
            return( str(min(d)),d[min(d)] )
        delta = delta*10
        if i > 1000:
            return("","")

def find_snps(chr,pos,c,dist):
    delta=10000
    i =0
    closest_snp=[]
    snp_distance=[]
    while True:
        i +=1
        A='SELECT pos,ref,alt FROM SVDB WHERE chr == \'{}\' AND pos > {} AND pos < {} '.format(chr,int(pos)-delta,int(pos)+delta)
        d={}
        for hit in c.execute(A):
            snp_distance.append( abs(pos- int( hit[0] ))  )
            if hit[0] <= dist:
                closest_snp.append( "{}:{},{}->{}".format(chr,hit[0],str(hit[1]),str(hit[2])) )
        if snp_distance:
            return( str(min(snp_distance)),"|".join(closest_snp) )
        delta = delta*10
        if i > 1000:
            return("","")

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

def read_cigar(cigar,contig_len):
    deletions=0
    insertions=0
    SC = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
    length=0
    first=True
    clip_after=True

    aligned_range=[]
    current_pos=1
    for i in range(0,len(SC)/2):
        if first and SC[i*2+1] == "M":
            first = False
        elif first and SC[i*2+1] == "S":
            first = False
            clip_after=False
        if SC[i*2+1] == "M":
            length += int( SC[i*2] )
            bases=range(0,int( SC[i*2] ))
            for j in range(0,len(bases)):
                bases[j] += current_pos

            aligned_range += bases
            current_pos += int( SC[i*2] )
       	elif SC[i*2+1] == "I":
       	    insertions+=1
            length += int( SC[i*2] )
            bases=range(0,int( SC[i*2] ))
            for j in range(0,len(bases)):
                bases[j] += current_pos
            aligned_range += bases

            current_pos += int( SC[i*2] )
       	elif SC[i*2+1] == "D":
       	    deletions +=1
        else:
            current_pos += int( SC[i*2] )

    return deletions,insertions,length,clip_after,aligned_range

def retrieve_pos(args,input_file):
    sucess=False
    contig=""
    homology_seq=""
    insertion_seq = ""

    deletions=0
    homology=0
    insertions=0

    for line in open(input_file):
        if line[0] == "@":
            continue
        content= line.strip().split("\t")
        

        contig=content[9]
        AB=True
        if not content[2] == args.chrA and not content[2] == args.chrB:
            continue
        pos = int(content[3])
        if pos >= args.posA-args.padding and pos <= args.posA+args.padding:
            pass
        elif pos >= args.posB-args.padding and pos <= args.posB+args.padding:
            AB=False
        else:
            continue
        
        if not "SA:Z:" in line:
            continue
        SA_line=line.strip().split("SA:Z:")[-1].split("\t")[0]
        SA_fields=SA_line.strip(";").split(";")
        found = False
        for SA in SA_fields:
            deletions=0
            insertions=0
            #print SA
            split_read=SA.split(",")
            #print split_read
            pos = int(split_read[1])

            
            cigar_del,cigar_ins,length,clip_after,range_secondary=read_cigar(split_read[3],len(contig))
            SA_orientation=split_read[2]
            SA_start=pos
            SA_end=length+pos-1
            deletions += cigar_del
            insertions += cigar_ins
            SA_clip_after= clip_after
            SA_len=length

            #print  "{} {} {}".format(pos,startB,stopB)
            if split_read[0] == args.chrB and pos >= args.posB-args.padding and pos <= args.posB+args.padding and AB:
                found=True
                break
            elif split_read[0] == args.chrA and pos >= args.posA-args.padding and pos <= args.posA+args.padding and not AB:
                found=True
                break

        if not found:
            continue
        flag="{0:012b}".format(int(content[1]))
        #the read is only accepted if it is a primary alignemnt, not a secondary alignment, and if it is only split in two parts
        if not int(flag[-9]) and not int(flag[0]) and len(SA_fields) == 1:
            orientationA="+"
            orientationB=SA_orientation
            if int(flag[-5]) :
                orientationA="-"
            cigar_del,cigar_ins,length,clip_after, range_primary=read_cigar(content[5], len(contig))
            deletions += cigar_del
            insertions += cigar_ins

            posA=int(content[3])+length-1
            posB=SA_start
            if clip_after:
               posA=int(content[3])+length-1
            else:
               posA=int(content[3])

            if SA_clip_after:
               posB=SA_end
            else:
               posB=SA_start


            if orientationB != orientationA:
                for i in range(0,len(range_secondary)):
                    range_secondary[i]=len(contig)+1-range_secondary[i]
                range_secondary=sorted(range_secondary)

            homology=len( set(range_primary).intersection(set(range_secondary))  )
            homology_seq=""
            homologous_pos=sorted(list(set(range_primary).intersection(set(range_secondary))))
            for i in range(0,len(homologous_pos)):
                if i == len(homologous_pos) -1:
                    homology_seq +=   contig[ homologous_pos[i]-1 ]
                else:
                    if homologous_pos[i] +1 == homologous_pos[i+1]:
                        homology_seq +=   contig[ homologous_pos[i]-1 ]
                    else:
                        homology_seq +=   contig[ homologous_pos[i]-1 ]  + ","                    



            insertion_seq=""
            insertion_range=sorted( list(set(range(1,len(contig)+1)).difference( set(range_secondary).union( set(range_primary) ) )))
            if insertion_range:
                insertions = len(insertion_range)
            for i in range(0,len(insertion_range)):
                if i == len(insertion_range) -1:
                    insertion_seq +=   contig[ insertion_range[i]-1 ]
                else:
                    if insertion_range[i] +1 == insertion_range[i+1]:
                        insertion_seq +=   contig[ insertion_range[i]-1 ]
                    else:
                        
                        insertion_seq +=   contig[ insertion_range[i]-1 ]  + ","                    
            


            contigA=""
            for i in range(0,len(range_primary)):
                contigA += contig[range_primary[i]-1]
            
            contigB=""
            reverse_comp={"A":"T","a":"t","T":"A","t":"a","G":"C","g":"c","C":"G","c":"g"}

            for i in range(0,len(range_secondary)):
                contigB += contig[range_secondary[i]-1]

            if orientationA == "-":
                tmpA=""
                for i in range(0,len(contigA)):
                    tmpA += reverse_comp[contigA[len(contigA) -i-1 ] ]
                contigA=tmpA
                tmpContig=""

                tmphomology=""
                if not "," in homology_seq:
                    for i in range(0,len(homology_seq)):
                        tmphomology += reverse_comp[homology_seq[len(homology_seq) -i-1 ] ]
                else:
                    homologous_sequences=homology_seq.split(",")
                    for seq in homologous_sequences:
                        for i in range(0,len(seq)):
                            if not i == 0:
                                tmphomology += ","
                            tmphomology += reverse_comp[seq[len(seq) -i-1 ] ]
                homology_seq=tmphomology

                for i in range(0,len(contig)):
                    tmpContig += reverse_comp[contig[len(contig)-i-1]]
                contig = tmpContig

            if orientationB == "-" and orientationA == "-":
                tmpB=""
                for i in range(0,len(contigB)):
                   tmpB += reverse_comp[contigB[len(contigB)-i-1]]
                contigB=tmpB             
 


            fontseq = xlwt.easyfont('')
            fontA= xlwt.easyfont('color_index green')
            fontHOM= xlwt.easyfont('color_index red')
            fontSEQ = xlwt.easyfont('color_index orange')
            fontB= xlwt.easyfont('color_index blue')
            seq_norm=""
            tupleB=[]
            tupleA=[]
            tupleCtg=[]
            tupleH=[]

            for i in range(0,len(contigA)):
                pos = range_primary[i]
                if orientationA == "-":
                    pos=range_primary[len(range_primary) -1 -i]

                if pos in homologous_pos:
                    tupleA.append( (contigA[i],fontHOM) )
                else:
                    tupleA.append( (contigA[i],fontA) )

            for i in range(0,len(contigB)):
                pos = range_secondary[i]
                if orientationB == "-" and orientationA == "-":
                    pos=range_secondary[len(range_secondary) -1 -i]

                if pos in homologous_pos:
                    tupleB.append( (contigB[i],fontHOM) )
                else:
                    tupleB.append( (contigB[i],fontB) )

            for i in range(0,len(contig)):
                pos= i +1
                if orientationA == "-":
                    pos=len(contig)-i
                if pos in homologous_pos:
                    tupleCtg.append( (contig[i],fontHOM) )
                elif pos in range_primary:
                    tupleCtg.append( (contig[i],fontA) )
                elif pos in range_secondary:
                    tupleCtg.append( (contig[i],fontB) )
                else:
                    tupleCtg.append( (contig[i],fontSEQ) )

            for i in range(0,len(homology_seq)):
                tupleH.append( (homology_seq[i],fontHOM) )
            
            sucess = True
            break

    if sucess:
        if AB:
            args.posA=posA
            args.orientationA=orientationA
            args.posB=posB
            args.orientationB=orientationB
            args.lengthA=length
            args.lengthB=SA_len
            args.regionA=contigA
            args.regionB=contigB
            args.regionAsegments= tuple(tupleA)
            args.regionBsegments= tuple(tupleB)
            args.contigSegments= tuple(tupleCtg)

        else:
            args.posB=posA
            args.orientationB=orientationA
            args.posA=posB
            args.orientationA=orientationB
            args.lengthA=SA_len
            args.lengthB=length
            args.regionA=contigB
            args.regionB=contigA
            args.regionAsegments= tuple(tupleB)
            args.regionBsegments= tuple(tupleA)
            args.contigSegments= tuple(tupleCtg)
        args.HomologySegments=tuple(tupleH)



    return (args,sucess,contig,homology,homology_seq,insertions,insertion_seq,deletions)


def extract_splits(args,ws0):
    if args.bed:
        input_file=args.bed
    elif args.vcf:
        input_file=args.vcf
    else:
        print "error: missing bed or vcf"
        quit()

    if not args.sample:
        args.sample=args.bam.split("/")[-1].split(".")[0]
    if not args.working_dir:
        args.working_dir=args.bam.split("/")[-1].split(".")[0]

    if args.repeatmask:
        conn = sqlite3.connect(args.repeatmask)
        c = conn.cursor()

    if args.snps:
        conn_snp=sqlite3.connect(":memory:")
        c_snp = conn_snp.cursor()

        A="CREATE TABLE SVDB (chr TEXT, pos INT,ref TEXT, alt TEXT)"
        c_snp.execute(A)

        input_vcf=[]

        if args.snps.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open


        for line in opener(args.snps):
            if line[0] == "#":
                continue
            content=line.strip().split()
            input_vcf.append([ content[0].replace("chr",""), content[1] , content[3] ,content[4] ])                
            if len(input_vcf) > 1000000:
                c_snp.executemany('INSERT INTO SVDB VALUES (?,?,?,?)',input_vcf)          
                input_vcf=[]
        if input_vcf:
            c_snp.executemany('INSERT INTO SVDB VALUES (?,?,?,?)',input_vcf)
                    

        A="CREATE INDEX SNP ON SVDB (chr, pos)"
        c_snp.execute(A)
        conn_snp.commit()

        

    row=1
    detected_splits={}    


    i = 0
    for line in open(input_file):
        if line[0] == "#":
            continue
        
        if args.bed:
            content=line.strip().split()
            args.chrA= content[0]
            args.posA= int(content[1])
            args.chrB= content[2]
            args.posB= int(content[3])
            args.type=content[4]
            args.orientationA=""
            args.orientationB=""
            args.id=str(i)
            var_id=str(i)
            args.lengthA=""
            args.lengthB=""
            args.regionA=""
            args.regionB=""
            insertion_seq = ""
            homology_seq = ""
            args.regionAsegments= ()
            args.regionBsegments= ()
            args.contigSegments= ()
            args.HomologySegments = ()
            args.repeatA= ""
            args.repeatB= ""
            i+=1
        elif args.vcf:
            chrA,posA,chrB,posB,event_type,INFO,FORMAT = readVCF.readVCFLine(line)
            args.chrA= chrA
            args.posA= int(posA)
            args.chrB= chrB
            args.posB= int(posB)
            args.type= event_type

            args.orientationA=""
            args.orientationB=""
            args.id=line.strip().split("\t")[2]
            var_id=line.strip().split("\t")[2]
            args.lengthA=""
            args.lengthB=""
            args.regionA=""
            args.regionB=""
            insertion_seq = ""
            homology_seq = ""
            args.regionAsegments= ()
            args.regionBsegments= ()
            args.contigSegments= ()
            args.HomologySegments = ()
            args.repeatA= ""
            args.repeatB= ""
            i+=1


        found=FindTranslocations.main(args)

        splits=0
        bp_homology=""
        insertions=""
        deletions=""
        sucess = False

        if args.skip_assembly and found:
            wd=os.path.join(args.working_dir,var_id)
            target = open(args.working_dir + "/" + var_id +"/" + "softclip.fa", 'w')
            for line in open(os.path.join(args.working_dir,var_id,"splits.sam")):
                content= line.strip().split("\t")
                flag="{0:012b}".format(int(content[1]))
                if not int(flag[-9]) and not int(flag[0]) and not int(flag[1]) and not int(flag[2]):
                   target.write( ">" + content[0] + "\n")
                   target.write(content[9] + "\n")
                   splits += 1
            target.close()

            try:
                args,sucess,contig,bp_homology,homology_seq,insertions,insertion_seq,deletions = retrieve_pos(args,os.path.join(args.working_dir,var_id,"splits.sam"))
            except:
                homology_seq="WARNING:unable to determine the breakpoint sequence"   

        elif found:
            wd=os.path.join(args.working_dir,var_id)
            target = open(args.working_dir + "/" + var_id +"/" + "softclip.fa", 'w')
            for line in open(os.path.join(args.working_dir,var_id,"splits.sam")):
                content= line.strip().split("\t")
                flag="{0:012b}".format(int(content[1]))
                if not int(flag[-9]) and not int(flag[0]) and not int(flag[1]) and not int(flag[2]):
                   target.write( ">" + content[0] + "\n")
                   target.write(content[9] + "\n")
                   splits += 1
            target.close()
            trials=[20,60,90]
            for k in trials:
                os.system("ABYSS -c {} -e {} -k {} -o {}_{}.fa {} > /dev/null 2>&1".format(1,0,k,os.path.join(args.working_dir,var_id,"abyss"),k,os.path.join(wd,"softclip.fa") ))
            os.system("cat {}_20.fa {}_60.fa {}_90.fa > {}".format(os.path.join(args.working_dir,var_id,"abyss"),os.path.join(args.working_dir,var_id,"abyss"),os.path.join(args.working_dir,var_id,"abyss")   ,os.path.join(args.working_dir,var_id,"abyss.fa")))

            if not os.stat( os.path.join(args.working_dir,var_id,"abyss.fa") ).st_size == 0:
                os.system("bwa mem {} {} > {}".format(args.fa,os.path.join(args.working_dir,var_id,"abyss.fa"),os.path.join(wd,"aligned_contig.sam")))
                try:
                    args,sucess,contig,bp_homology,homology_seq,insertions,insertion_seq,deletions = retrieve_pos(args,os.path.join(args.working_dir,var_id,"aligned_contig.sam"))
                except:
                    pass
            if not sucess:

                try:
                    args,sucess,contig,bp_homology,homology_seq,insertions,insertion_seq,deletions = retrieve_pos(args,os.path.join(args.working_dir,var_id,"splits.sam"))
                except:
                    homology_seq="WARNING:unable to determine the breakpoint sequence"
        else:
            wd=os.path.join(args.working_dir,var_id)
            os.system("samtools view -bh {} {}:{}-{} > {}/regionA.bam".format(args.bam,args.chrA,args.posA-args.padding,args.posA+args.padding,wd))
            os.system("samtools view -bh {} {}:{}-{} > {}/regionB.bam".format(args.bam,args.chrB,args.posB-args.padding,args.posB+args.padding,wd))
            os.system("samtools merge -f {}/region.bam {}/regionA.bam {}/regionB.bam ".format(wd,wd,wd))
            os.system("samtools bam2fq {}/region.bam > {}/region.fq".format(wd,wd))
            trials=[20,60,90]
            for k in trials:
                os.system("ABYSS -c {} -e {} -k {} -o {}_{}.fa {} > /dev/null 2>&1".format(1,10,k,os.path.join(args.working_dir,var_id,"abyss"),k,os.path.join(wd,"region.fq") ))
            os.system("cat {}_20.fa {}_60.fa {}_90.fa > {}".format(os.path.join(args.working_dir,var_id,"abyss"),os.path.join(args.working_dir,var_id,"abyss"),os.path.join(args.working_dir,var_id,"abyss")   ,os.path.join(args.working_dir,var_id,"abyss.fa")))

            if not os.stat( os.path.join(args.working_dir,var_id,"abyss.fa") ).st_size == 0:
                os.system("bwa mem {} {} > {}".format(args.fa,os.path.join(args.working_dir,var_id,"abyss.fa"),os.path.join(wd,"aligned_contig.sam")))
                try:
                    args,sucess,contig,bp_homology,homology_seq,insertions,insertion_seq,deletions = retrieve_pos(args,os.path.join(wd,"aligned_contig.sam"))
                except:
                    homology_seq="WARNING:unable to determine the breakpoint sequence"
            if not sucess:
                contig=""

        distanceA=""
        distanceB=""
        if args.repeatmask:
            distanceA,args.repeatA= find_repeat(args.chrA,args.posA,c)
            distanceB,args.repeatB= find_repeat(args.chrB,args.posB,c)

        snpDistanceA=""
        snpDistanceB=""
        snpsA=""
        snpsB=""
        if args.snps:
            snpDistanceA,snpsA= find_snps(args.chrA,args.posA,c_snp,args.snp_distance)
            snpDistanceB,snpsB= find_snps(args.chrB,args.posB,c_snp,args.snp_distance)

        row_content=[args.sample,var_id,args.type,splits,args.chrA,args.posA,args.orientationA,args.repeatA,distanceA,args.chrB,args.posB,args.orientationB,args.repeatB,distanceB,bp_homology,args.HomologySegments,insertions,insertion_seq,args.lengthA,args.lengthB,len(contig),args.regionAsegments,args.regionBsegments,args.contigSegments]
        j=0
        for item in row_content:

            if j in [19,25,26,27]:
                ws0.write_rich_text(row, j, item)
            else:
                ws0.write(row, j, item)
            j+=1
        row += 1



parser = argparse.ArgumentParser("""SplitVision - SV breakpoint analysis software""")
parser.add_argument('--analyse',action="store_true",help="analyse breakpoints")
parser.add_argument('--db',action="store_true",help="generate the repeatmask database")
args, unknown = parser.parse_known_args()

if args.analyse:

    parser = argparse.ArgumentParser("""SplitVision - SV breakpoint analysis software""")
    parser.add_argument('--analyse',action="store_true",help="analyse breakpoints")
    parser.add_argument('--vcf'        , type=str, help="input vcf file containing breakpoints of interest(use only bed or vcf at a time)")
    parser.add_argument('--bed', type=str, help="input bed file(tab separted) containing the sv breakpoints(format: chrA,posA,chrB,posB)")
    parser.add_argument('--bam', type=str,required=True ,help="the input bam file")
    parser.add_argument('--fa', type=str,required=True ,help="the reference fasta file")
    parser.add_argument('--skip_assembly',required=False, action="store_true",help="skip the assembly analysis")
    parser.add_argument('--working_dir', type=str ,help="working directory")
    parser.add_argument('--sample', type=str ,help="sample id")
    parser.add_argument('--snps', type=str ,help="a vcf file containing snps, the software will compute the distance to the closest snp, and report snps within the snp_distance")
    parser.add_argument('--snp_distance', type=int,default=100 ,help="report snps within this distance ")
    parser.add_argument('--repeatmask', type=str,help="database file generated from the uscs repeat mask")
    parser.add_argument('--padding', type=int,default=1000 ,help="search for reads mapped within this distance fromt the breakpoint position")
    args = parser.parse_args()
    
    wb =  xlwt.Workbook()
    ws0 = wb.add_sheet("SplitVision",cell_overwrite_ok=True)
    
    header=["sampleID","variant_id","variant_type","split_reads","ChrA","PosA","OrientationA","repeatA","repeat_distanceA","snps/indelsA","snp_distance_A","ChrB","PosB","OrientationB","repeatB","repeat_distanceB","snps/indelsB","snp_distance_B","breakpoint_microhomology(bp)","breakpoint_microhomology(sequence)","insertions(bp)","insertions(sequence)","lengthA","lengthB","contig_length","regionA_sequence","regionB_sequence","contig_sequence"]
    j=0
    for item in header:
        ws0.write(0, j, item)
        j += 1
    detected_splits=extract_splits(args,ws0)
    wb.save(os.path.join(args.working_dir,args.sample+".xls"))

elif args.db:
    parser = argparse.ArgumentParser("""generate a database of genomic regions""")
    parser.add_argument('--db',action="store_true",help="generate the databases from input bed files")
    parser.add_argument('--prefix',type=str,required=True,help="the prefix of the database")
    parser.add_argument('--tab',type=str,required=True,help="the input repeatmask file")
    args, unknown = parser.parse_known_args()
    db(args)
else:
    print("invalid option, use --analyse to analyse breakpoints, or --db to generate a repeatmasksk db file")



