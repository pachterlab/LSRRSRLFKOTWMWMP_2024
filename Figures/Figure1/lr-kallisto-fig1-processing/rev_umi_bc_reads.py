import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--fastq', help='reverse the sequences and quality scores', type=argparse.FileType('r'))
parser.add_argument('--ofastq', help='output reversed fastq sequences here', type=argparse.FileType('w'))

args = parser.parse_args()
c = 0
writelines = ""
for l in args.fastq.readlines():
    c+=1
    if l.startswith("@") or l.startswith("+"):
        #args.ofastq.write(l)
        writelines+=l
    elif l != "\n":
        writelines+=l[16:24]+l[8:16]+l[0:8]+"\n"
        #args.ofastq.write(l[::-1])
    else:
        writelines+=l
    if c == 4:# and l != "\n":
        args.ofastq.write(writelines)
        c = 0
        writelines = ""
    '''
    elif c == 4:
        c = 0
        writelines = ""
    '''

