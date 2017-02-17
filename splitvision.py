import argparse
import FindTranslocations

parser = argparse.ArgumentParser("""SplitVision - SV breakpoint analysis software""")
parser.add_argument('--vcf'        , type=str, help="input vcf file containing breakpoints of interest(use only bed or vcf at a time)")
parser.add_argument('--bed', type=str, help="input bed file(tab separted) containing the sv breakpoints(format: chrA,posA,chrB,posB)")
parser.add_argument('--bam', type=str,required=True ,help="the input bam file")
parser.add_argument('--working_dir', type=str,default="out" ,help="working directory")
parser.add_argument('--padding', type=int,default=1000 ,help="search for reads mapped within this distance fromt the breakpoint position")
args = parser.parse_args()

i=0
def extract_splits(args):
    detected_splits={}    

    if args.bed:
        input_file=args.bed
    elif args.vcf:
        input_file=args.vcf
    else:
        print "error: missing bed or vcf"

    for line in open(input_file):
        if line[0] == "#":
            continue
        
        if args.bed:
            content=line.strip().split()
            args.chrA= content[0]
            args.posA= content[1]
            args.chrB= content[2]
            args.posB= content[3]
            args.id=str(i)
            i += 1
        elif args.vcf:
            pass

        found=FindTranslocations.main(args)
        detected_splits["var_id"]={}
        detected_splits["var_id"]["splits_found"]=found





detected_splits=extract_splits(args)
