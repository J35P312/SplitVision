import argparse
import FindTranslocations
import itertools
import os

parser = argparse.ArgumentParser("""SplitVision - SV breakpoint analysis software""")
parser.add_argument('--vcf'        , type=str, help="input vcf file containing breakpoints of interest(use only bed or vcf at a time)")
parser.add_argument('--bed', type=str, help="input bed file(tab separted) containing the sv breakpoints(format: chrA,posA,chrB,posB)")
parser.add_argument('--bam', type=str,required=True ,help="the input bam file")
parser.add_argument('--fa', type=str,required=True ,help="the reference fasta file")
parser.add_argument('--working_dir', type=str,default="out" ,help="working directory")
parser.add_argument('--sample', type=str,default="patient" ,help="sample id")
parser.add_argument('--coverage', type=int,default=20 ,help="library coverage, does not need to be exact +-15 or so will do")
parser.add_argument('--padding', type=int,default=1000 ,help="search for reads mapped within this distance fromt the breakpoint position")
args = parser.parse_args()

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
            current_pos += int( SC[i*2] )
        else:
            current_pos += int( SC[i*2] )

    return deletions,insertions,length,clip_after,aligned_range

def retrieve_pos(args,input_file):
    sucess=False
    contig=""

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
        if not int(flag[-9]):
            orientationA="+"
            orientationB=SA_orientation
            if int(flag[-5]) :
                orientationA="-"
            cigar_del,cigar_ins,length,clip_after, range_primary=read_cigar(content[5], len(contig))
            deletions += cigar_del
            insertions += cigar_ins

            posA=int(content[3])+length-1
            posB=SA_start
            if orientationA == "+":
                if clip_after:
                    posA=int(content[3])+length-1
                else:
                    posA=int(content[3])
            else:
                if clip_after:
                    posA=int(content[3])
                else:
                    posA=int(content[3])+length-1

            if orientationB == "+":
                if SA_clip_after:
                    posB=SA_end
                else:
                    posB=SA_start
            else:
                if SA_clip_after:
                    posB=SA_start
                else:
                    posB=SA_end

            if orientationB != orientationA:
                for i in range(0,len(range_secondary)):
                    range_secondary[i]=len(contig)+1-range_secondary[i]

            homology=len( set(range_primary).intersection(set(range_secondary))  )

            contigA=""
            for i in range(0,len(range_primary)):
                contigA += contig[range_primary[i]-1]
            
            contigB=""
            reverse_comp={"A":"T","a":"t","T":"A","t":"a","G":"C","g":"c","C":"G","c":"g"}
            for i in range(0,len(range_secondary)):
                contigB += contig[range_secondary[i]-1]       
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

        else:
            args.posB=posA
            args.orientationB=orientationA
            args.posA=posB
            args.orientationA=orientationB
            args.lengthA=SA_len
            args.lengthB=length
            args.regionA=contigB
            args.regionB=contigA
    return (args,sucess,contig,homology,insertions,deletions)


def extract_splits(args):
    detected_splits={}    

    if args.bed:
        input_file=args.bed
    elif args.vcf:
        input_file=args.vcf
    else:
        print "error: missing bed or vcf"
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
            args.orientationA=""
            args.orientationB=""
            args.id=str(i)
            var_id=str(i)
            args.lengthA=""
            args.lengthB=""
            args.regionA=""
            args.regionB=""
            i+=1
        elif args.vcf:
            print "to be implemented, please use bed for now"
            break

        found=FindTranslocations.main(args)

        splits=0
        bp_homology=""
        insertions=""
        deletions=""

        if found:
            target = open(args.working_dir + "/" + var_id +"/" + "softclip.fa", 'w')
            for line in open(os.path.join(args.working_dir,var_id,"splits.sam")):
                content= line.strip().split("\t")
                flag="{0:012b}".format(int(content[1]))
                if not int(flag[-9]):
                   target.write( ">" + content[0] + "\n")
                   target.write(content[9] + "\n")
                   splits += 1
            target.close()
            trials=[20,90]
            for k in trials:
                os.system("ABYSS -c 1 -e 0 -k {} -o {} {} > /dev/null 2>&1".format(k,os.path.join(args.working_dir,var_id,"abyss.fa"),os.path.join(args.working_dir,var_id,"softclip.fa") ))
                os.system("bwa mem {} {} > {}".format(args.fa,os.path.join(args.working_dir,var_id,"abyss.fa"),os.path.join(args.working_dir,var_id,"aligned_contig.sam")))
                args,sucess,contig,bp_homology,insertions,deletions = retrieve_pos(args,os.path.join(args.working_dir,var_id,"aligned_contig.sam"))
                if sucess:
                    break

            if not sucess:
                args,sucess,contig,bp_homology,insertions,deletions = retrieve_pos(args,os.path.join(args.working_dir,var_id,"splits.sam"))
        else:
            wd=os.path.join(args.working_dir,var_id)
            os.system("samtools view -bh {} {}:{}-{} > {}/regionA.bam".format(args.bam,args.chrA,args.posA-args.padding,args.posA+args.padding,wd))
            os.system("samtools view -bh {} {}:{}-{} > {}/regionB.bam".format(args.bam,args.chrB,args.posB-args.padding,args.posB+args.padding,wd))
            os.system("samtools merge -f {}/region.bam {}/regionA.bam {}/regionB.bam ".format(wd,wd,wd))
            os.system("samtools bam2fq {}/region.bam > {}/region.fq".format(wd,wd))
            trials=[20,40,60,90]
            for k in trials:
                os.system("ABYSS -c {} -e {} -k {} -o {} {} > /dev/null 2>&1".format(1,10,k,os.path.join(args.working_dir,var_id,"abyss.fa"),os.path.join(wd,"region.fq") ))
                os.system("bwa mem {} {} > {}".format(args.fa,os.path.join(args.working_dir,var_id,"abyss.fa"),os.path.join(wd,"aligned_contig.sam")))
                args,sucess,contig,bp_homology,insertions,deletions = retrieve_pos(args,os.path.join(wd,"aligned_contig.sam"))
                if sucess:
                    break
            if not sucess:
                contig=""

        print"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(args.sample,var_id,splits,args.chrA,args.posA,args.orientationA,args.chrB,args.posB,args.orientationB,bp_homology,insertions,deletions,args.lengthA,args.lengthB,len(contig),args.regionA,args.regionB,contig).strip()

header=["sampleID","variant_id","split_reads","ChrA","PosA","OrientationA","ChrB","PosB","OrientationB","breakpoint_homology(bp)","insertions","deletions","regionA_sequence","regionB_sequwnce","lengthA","lengthB","contig_length","contig_sequence"]
print "\t".join(header)      
detected_splits=extract_splits(args)
