import sys


def reverse_read(read):
    reversed=""
    reverse={"A":"T","T":"A","G":"C","C":"G"}
    for i in range(0,len(read)):

        reversed+=reverse[read[-1-i].upper()]

    return(reversed)

reads={}
for line in open(sys.argv[1]):
    content=line.strip().split()
    flag="{0:012b}".format(int(content[1]))
    if not int(flag[0]) and not int(flag[1]) and not int(flag[2]) and not int(flag[3]) and not int(flag[4]):
       reads[content[0]]=content[9]
       if not int(flag[7]):
          reads[content[0]]=reverse_read(content[9])

       if not int(flag[5]):
          reads[content[0]]=reverse_read(content[9])

for read in reads:
    print ">{}".format(read)
    print reads[read]
