import os
import sys
#import readVCF
import itertools

def get_sam_data(line):
    content=line.strip().split("\t")
    start=int(content[3])
    end = int(content[3])
    SC = ["".join(x) for _, x in itertools.groupby(content[5], key=str.isdigit)]



    length=0
    for i in range(0,len(SC)/2):
        if SC[i*2+1] == "M":
            length += int( SC[i*2] )
	end+=length
		    
    sam={"ID":content[0],"chr":content[2],"start":start,"end": end,"Q":int(content[4]),"CIGAR":content[5],"length":length}
    return(sam)

def find_variants(bam,chrA,startA,stopA,chrB,startB,stopB,working_dir,it):
    bam_prefix=bam.split("/")[-1]
    prefix=working_dir + "/" +  bam_prefix[0:-4]
    os.system("samtools view -h -F 256  {} {}:{}-{} | samtools view -Sh -F 512 - | samtools view -Sh -F 1024 - | samtools view -S -F 2048 - | grep S  > {}_test.sam".format(bam,chrA,startA,stopA,prefix))
    splits=[]        
    found = False
    for line in open("{}_test.sam".format(prefix)):
        content=line.strip().split()
        if "SA:Z" in line and not "H" in content[5]: 
            SA_line=line.strip().split("SA:Z:")[-1].split("\t")[0]
            SA_fields=SA_line.strip(";").split(";")
            for SA in SA_fields:
                #print SA
                split_read=SA.split(",")
                #print split_read
                pos = int(split_read[1])
                #print  "{} {} {}".format(pos,startB,stopB)
                if split_read[0] == chrB and pos >= startB and pos <= stopB:
                    found = True
                    splits.append(line.strip() +"\n")
                    continue
           
                SC = ["".join(x) for _, x in itertools.groupby(SA, key=str.isdigit)]
                length=0
                for i in range(0,len(SC)/2):
                    if SC[i*2+1] == "M":
                        length += int( SC[i*2] )

                if split_read[0] == chrB and pos+length >= startB and pos+length <= stopB:
                    found = True
                    splits.append(line.strip() +"\n")
            
    os.system("samtools view -H {} > {}_header.sam".format(bam,prefix))

    target = open(working_dir +"/" + "splits.sam", 'w')
    splits=set(splits)
    for line in splits:
        target.write(line)
    target.close()
    if it == 1:
       os.system("cat {}_header.sam {}/splits.sam | samtools view -Shb - > {}/splits_1.bam".format(prefix,working_dir,working_dir  ))
    elif it == 2:
       os.system("cat {}_header.sam {}/splits.sam | samtools view -Shb - > {}/splits_2.bam".format(prefix,working_dir,working_dir  ))
       os.system("samtools merge -f {}/splits.bam {}/splits_1.bam {}/splits_2.bam".format(working_dir, working_dir, working_dir))
       os.system("samtools view {}/splits.bam > {}/splits.sam".format(working_dir,working_dir))

       unique_alignments=[]
       for line in open("{}/splits.sam".format(working_dir)):
          unique_alignments.append(line.strip() )

       unique_alignments=set(unique_alignments)
       unique_aln_sam = open(working_dir +"/" + "splits.sam", 'w')
       for line in unique_alignments:
          unique_aln_sam .write(line + "\n")
       unique_aln_sam.close()

       #os.system("cat {}_header.sam {}/splits.sam | samtools view -Shb - > {}/splits.bam".format(prefix,working_dir,working_dir  ))
       #os.system("samtools bam2fq {}/splits.bam > {}/splits.fq".format(working_dir,working_dir))

    #os.system("samtools bam2fq {}/splits.bam > {}/splits.fastq".format(working_dir,working_dir) )
    return found

def main(args):
    if not os.path.isdir( os.path.join(args.working_dir, args.id ) ):
        os.makedirs( os.path.join(args.working_dir, args.id ) )

    padding=args.padding
    var_id=args.id
    chrA=args.chrA
    startA= args.posA -padding
    stopA= args.posA + padding
    
    chrB=args.chrB
    startB= args.posB - padding
    stopB= args.posB + padding
    working_dir=os.path.join(os.path.join(args.working_dir, args.id ))
    found=find_variants(args.bam,chrA,startA,stopA,chrB,startB,stopB,working_dir,1)
    if not found:
        #print 1
        found =find_variants(args.bam,chrB,startB,stopB,chrA,startA,stopA,working_dir,2)
    else:
       # print 
        tmp =find_variants(args.bam,chrB,startB,stopB,chrA,startA,stopA,working_dir,2)
        if tmp:
            found = True

    return found  
