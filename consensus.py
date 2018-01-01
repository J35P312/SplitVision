import sys
import os
clustal="clustalw2"
cdhit="cd-hit"

def cluster(cdhit,folder,softclipfa):
    #cluster the reads, and analyse the reads belonging to the largest cluster
    softclip_clusters=os.path.join(folder,"cdhit_softclip")
    os.system("{} -d 0 -i {} -o {} -c 0.85 -g 1 > /dev/null 2>&1".format(cdhit,softclipfa,softclip_clusters))
    clusters={}
    n_clusters=0
    for line in open(softclip_clusters+".clstr"):
       if line[0] == ">":
           n_clusters+=1
           clusters[n_clusters]=[]
           continue
       sequence_id=">"+line.strip().split(">")[-1].split("...")[0]
       clusters[n_clusters].append(sequence_id)
    softclip_dictionary={}
    for line in open(softclipfa):
        if line[0] == ">":
            current_id=line.strip()
            continue
        softclip_dictionary[current_id]=line.strip()

    largest_cluster_size=-1
    largest_cluster=""
    for cluster in clusters:
        cluster_size=len(clusters[cluster])
        if cluster_size > largest_cluster_size:
            largest_cluster_size=cluster_size
            largest_cluster=cluster

    new_fa=""
    for sequence in clusters[largest_cluster]:
        new_fa+="{}\n{}\n".format(sequence,softclip_dictionary[sequence])
    new_fa=new_fa.strip()

    softclipfa=os.path.join(folder,"softclip.clustered.fa")
    target = open( softclipfa, 'w')
    target.write(new_fa)
    target.close()
    return(largest_cluster_size)

def reverse_read(read):
    reversed=""
    reverse={"A":"T","T":"A","G":"C","C":"G","N":"N"}
    for i in range(0,len(read)):
        if read[-1-i].upper() in reverse:
           reversed+=reverse[read[-1-i].upper()]
        else:
           reversed+=N
    return(reversed)

def majority(seqs,i):
    nucleotides={"A":0,"T":0,"C":0,"G":0,"-":0}
    for seq in seqs:
        if seq[i] in nucleotides:
            nucleotides[seq[i]]+=1
    major="A"
    for nucleotide in nucleotides:
        if nucleotides[nucleotide] > nucleotides[major]:
            major=nucleotide
    return (major,nucleotides[major])


def main(folder,input_sam,name):
    splits=0
    softclipfa=os.path.join(folder,"softclip.fa")
    target = open( softclipfa, 'w')
    #extract split reads
    for line in open(input_sam):
        content= line.strip().split("\t")
       	flag="{0:012b}".format(int(content[1]))
        if not int(flag[0]) and not int(flag[1]) and not int(flag[2]) and not int(flag[3]):
            if "SA:Z" in line:
                split_pos=line.split("SA:Z:")[-1].split("\t")[0]
                chromosome_aln=content[2]
                #print line
                #print flag
                chromosome_split=split_pos.split(",")[0]

                pos_aln=int(content[3])
                pos_split=int(split_pos.split(",")[1])
                orientation_split= 0
                if split_pos.split(",")[2] == "-":
                     orientation_split = 1
                target.write( ">" + content[0]+":"+str(splits) + "\n")
                if int(flag[-5]) == orientation_split:
                    target.write(content[9] + "\n")
                else:
                    if chromosome_aln == chromosome_split and pos_aln < pos_split:
                       target.write(content[9] + "\n")
                    elif chromosome_aln < chromosome_split:
                       target.write(content[9] + "\n")
                    else:
                       target.write(reverse_read(content[9]) + "\n")
                splits += 1
    target.close()

    largest_cluster_size=cluster(cdhit,folder,softclipfa)
    softclipfa=os.path.join(folder,"softclip.clustered.fa")

    if largest_cluster_size == 1:
        sequence=""
        for line in open(softclipfa):
            if line[0] == ">":
               continue
            sequence= line.strip()
            break
        return(">{} {} {} {} {} {}\n{}".format( name,len(sequence),largest_cluster_size,splits,len(sequence), ":".join(["1"]*len(sequence)),sequence) )

    #generate the multiple sequence alignment
    consensusfa=os.path.join(folder,"consensus.fa")
    os.system("{} -INFILE={} -ALIGN > /dev/null 2>&1".format(clustal,softclipfa))
    #os.system("{} -INFILE={} -ALIGN".format(clustal,softclipfa))
    segments={}
    first=True
    aln_start=False
    it=0
    for line in open(softclipfa.replace(".fa",".aln")):
        if first:
            first=False
            it=1
            continue
        if line != "" and line[0] != "\n":
            
            if line[0] != " ":
                line=line.rstrip()
                if not it in segments:
                    segments[it]=[]
                segments[it].append(line.split()[1])
                block_end=len(line)
                if not aln_start:
                    for i in range(1,len(line)):
                        if line[-i] == " ":
                            aln_start=len(line)-i+1
                            break
            else:
                matches=[]
                for i in range(aln_start,len(line)-1):
                    if line[i] == "*":
                        matches.append("*")
                    else:
                        matches.append("-")
                segments[it].append("".join(matches))
                it+=1
    if not segments:
       return ("")

    full_alignments=[]
    for i in range(0,largest_cluster_size+1):
        full_alignments.append("")

    for i in range(1,len(segments)+1):
        for j in range(0,len(segments[i])):
            full_alignments[j]+=segments[i][j]

    for i in range(0,len(full_alignments)-1):
        len_before=len(full_alignments[i])
        full_alignments[i]=full_alignments[i].rstrip("-")
        len_after=len(full_alignments[i])
        full_alignments[i]+="#" * (len_before-len_after)

        full_alignments[i]=full_alignments[i].strip("-")
        len_after=len(full_alignments[i])
        full_alignments[i]="#" * (len_before-len_after) + full_alignments[i]

    sequence=[]
    scores=[]
    similarity=full_alignments[-1]
    identical_bases=0
    del full_alignments[-1]
    for i  in range(0,len(similarity)):
        if similarity[i] == "*":
            scores.append(str(len(full_alignments)))
            sequence.append( full_alignments[0][i] )
            identical_bases+=1
        else:
            #print segments[segment]
            n,s=majority(full_alignments,i)
            if not n == "-":
                scores.append(str(s))
                sequence.append(n)
    return ">{} {} {} {} {} {}".format( name,len(sequence),largest_cluster_size,splits,identical_bases, ":".join(scores) ) +"\n"+"".join(sequence)

print main(sys.argv[1],sys.argv[2],"consensus")

