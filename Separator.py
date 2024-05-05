import os
import sys
import Levenshtein
import shutil
import multiprocessing
from tqdm import tqdm

#============================
# Functions for data reading
#============================
def Read_Primers(path):
    ''' '''
    return [(x.split('\t')[0],x.split('\t')[1],RevComp(x.split('\t')[2]))for x in open(path,'r').read().split('\n') if x!= '']

def FromFastqGenerateSeqs(path,Size=1000):
    ''' '''
    with open(path, 'r') as file:
        lines = []
        for line in file:
            lines.append(line.strip())
            if len(lines) == 4*Size:
                Seqs=[x[1] for x in split_list(lines,4)]
                yield Seqs
                lines = []
        if lines:
            Seqs=[x[1] for x in split_list(lines,4)]
            yield Seqs

def Read_File(filename, batch_size=None):
    if batch_size is not None:
        with open(filename, 'r') as file:
            lines = []
            for line in file:
                lines.append(line.strip())
                if len(lines) == batch_size:
                    yield lines
                    lines = []
            if lines:
                yield lines
    else:
        with open(filename, 'r') as file:
            for line in file:
                yield line.strip()

class Read_DIR():
    def __init__(self,Dirpath,kw='.'):
        self.FNames = [x for x in os.listdir(Dirpath) if kw in x]
        self.FPaths = [os.path.join(Dirpath,x) for x in self.FNames]
    def GetFileNamesAndPaths(self):
        return self.FNames,self.FPaths

#===================================
# Functions for sequence processing
#===================================
complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N': 'N','-':'-'}
def RevComp(dna_sequence):
    return ''.join([complement_dict[base] for base in dna_sequence[::-1]])

def SubSeqF(seq, k): # First half of the sequence
    return [seq[i:i+k] for i in range(len(seq)//2-k+1)]

def SubSeqR(seq, k): # Second half of the sequence
    return [seq[i:i+k] for i in range(len(seq)//2,len(seq)//2+len(seq)//2-k+1)]

def FindTagF(Seq,Tag,Threshold=2):   
    record,subseqlst = [],SubSeqF(Seq, len(Tag))
    for idx, subseq in enumerate(subseqlst):
        rst = Levenshtein.distance(subseq, Tag)
        if rst <= Threshold:
            record.append((idx, rst))
            if len(record) > 1 and rst <= record[-2][1]:
                return record[-2][0]
    if record:
        return record[-1][0]
    return False

def FindTagR(Seq,Tag,Threshold=2): 
    record,subseqlst = [],SubSeqR(Seq, len(Tag))
    for idx, subseq in enumerate(subseqlst):
        rst = Levenshtein.distance(subseq, Tag)
        if rst <= Threshold:
            record.append((idx-len(subseqlst), rst))
            if len(record) > 1 and rst <= record[-2][1]:
                return record[-2][0]
    if record:
        return record[-1][0]
    return False 

def LabelSeq(seq,primers,Threshold=1,trim=0):
    for primer in primers:
        label,Fprimer,Rprimer = primer[0],primer[1],primer[2]
        idxF=FindTagF(seq,Fprimer,Threshold)
        if idxF:
            idxR=FindTagR(seq,Rprimer,Threshold)
            if idxR:
                return label,seq[idxF+trim:idxR-trim]
        else:
            seq=RevComp(seq)
            idxF=FindTagF(seq,Fprimer,Threshold)
            if idxF:
                idxR=FindTagR(seq,Rprimer,Threshold)
                if idxR:
                    return label,seq[idxF+trim:idxR-trim]
    return None

#=================
# Other Functions
#=================
def split_list(lst, n=4):
    '''Merge elements of a list into specified number of groups'''
    return [lst[i:i+n] for i in range(0, len(lst), n)]


#===============
# Core Function
#===============
def CoreRun(Seqs,Primers,OutPath,MinLenth,MisMatch,TagTrim):
    results=[]
    for seq in tqdm(Seqs):
        results.append(LabelSeq(seq,Primers,MisMatch,TagTrim))
    f=open(OutPath,'w',encoding='utf-8')
    [print(x[0]+'\t'+x[1],file=f) for x in results if x != None and len(x[1])>=MinLenth]
    f.close()

#================
# Final Function
#================
def Final(Input_DIR,OTPT_DIR):
    File_LST=Read_DIR(Input_DIR).FPaths
    DICT={}
    for file in File_LST:
        for line in Read_File(file):
            gene_pop=line.split('\t')[0]
            seq=line.split('\t')[1]
            if gene_pop not in DICT.keys():
                DICT[gene_pop]=[]
            else:
                DICT[gene_pop].append(seq)      
    try:
        os.mkdir(OTPT_DIR)
    except FileExistsError:
        shutil.rmtree(OTPT_DIR)
        os.mkdir(OTPT_DIR)
        
    for gene_pop in tqdm(DICT.keys()):
        O_file=os.path.join(OTPT_DIR,gene_pop+'.fasta')
        SeqLst=DICT[gene_pop]
        with open(O_file,'w',encoding='utf-8') as f:
            [print(f">{gene_pop}_{str(i+1)}\n{SeqLst[i]}",file=f) for i in range(len(SeqLst))]    



if __name__=='__main__':
    
    # Default:
    MinLenth= 2000
    MisMatch= 1
    TagTrim=0
    CPUnum=3
    
    # Input
    for i in range(len(sys.argv)):
        j=i+1
        if '--Fastq' == sys.argv[i]:
            fastq=str(sys.argv[j])
        elif '--Primer' == sys.argv[i]:
            primer=str(sys.argv[j])
        elif '--OutDir' == sys.argv[i]:
            outdir=str(sys.argv[j])
        elif '--MinLenth' == sys.argv[i]:
            MinLenth=int(sys.argv[j])
        elif '--MisMatch' == sys.argv[i]:
            MisMatch=int(sys.argv[j])
        elif '--TagTrim' == sys.argv[i]:
            TagTrim=int(sys.argv[j])
        elif '--CPU' == sys.argv[i]:
            CPUnum=int(sys.argv[j])
    


    try:
        os.mkdir(outdir+"_temp")
    except FileExistsError:
        shutil.rmtree(outdir+"_temp")
        os.mkdir(outdir+"_temp")


    SEQs_BATCH,OUTs_PATH=[],[]
    primers=Read_Primers(primer)
    cnt=0
    ReadsNumber=0
    for Seqs in FromFastqGenerateSeqs(fastq):
        cnt+=1
        ReadsNumber+=len(Seqs)
        print(f"{ReadsNumber} reads have been loaded...")
        SEQs_BATCH.append(Seqs)
        OUTs_PATH.append(os.path.join(outdir+"_temp",f"LabelSeqSet_{str(cnt)}.txt"))
        if len(SEQs_BATCH) == CPUnum:
            print(f"Please be patient: the multiprocessing is currently in progress.")
            processes=[]
            for seqs,out in zip(SEQs_BATCH,OUTs_PATH):
                p=multiprocessing.Process(target=CoreRun,args=(seqs,primers,out,MinLenth,MisMatch,TagTrim))
                processes.append(p)
                p.start()
            for p in processes:
                p.join()
            SEQs_BATCH,OUTs_PATH=[],[]
    if SEQs_BATCH:
        print(f"Please be patient: the multiprocessing is currently in progress.")
        processes=[]
        for seqs,out in zip(SEQs_BATCH,OUTs_PATH):
            p=multiprocessing.Process(target=CoreRun,args=(seqs,primers,out,MinLenth,MisMatch,TagTrim))
            processes.append(p)
            p.start()
        for p in processes:
            p.join()
    print(f"Performing the final integration of all reads.")
    Final(outdir+"_temp",outdir)
    shutil.rmtree(outdir+"_temp")
    print(f"Done!")
    