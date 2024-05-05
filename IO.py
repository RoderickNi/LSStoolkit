# IO.py
import os
from Bio import SeqIO
from Other_Functions import split_list,RevComp
from statistics import mean


class Read_Fasta():
    def __init__(self,path):
        self.init = SeqIO.parse(path, "fasta")
        self.index = dict([
            (
            record.id,
            {
            'Seq':record.seq.__str__().upper(),
            'Len':len(record)
            }
            ) for record in self.init])
        self.path = path
    def GetNameSeq(self):
        return dict([(name,info['Seq'])for name,info in self.index.items()])
    def GetNameSeq_to_lst(self):
        return [[name,info['Seq']]for name,info in self.index.items()]
    def Get_NameLst_SeqLst(self):
        return [name for name in self.index.keys()],[info['Seq'] for info in self.index.values()]
    def Get_SeqLst(self):
        return [info['Seq'] for info in self.index.values()]
    def Get_NameLst(self):
        return [name for name in self.index.keys()]
    
        
class Read_Fastq():
    def __init__(self,path):
        self.init = SeqIO.parse(path, "fastq")
        self.index = dict([
            (
            record.id,
            {
            'Seq':record.seq.__str__().upper(),
            'Qual':record.letter_annotations["phred_quality"],
            'Len':len(record)
            }
            ) for record in self.init])
        self.path = path 
    def QC_stat(self):
        return {seq_name:{'Mean_Q':mean(info['Qual']),'Len':info['Len']} for seq_name,info in self.index.items()}
    
        
class Read_Aln():
    def __init__(self,path):
        self.init = split_list([x for x in open(path,'r',encoding='utf-8').read().split('\n') if x!=''],n=3)
        self.index = dict([
            (
            record[0].replace('@',''),
            {
            'Ref':record[1],
            'Seq':record[2]
            }
            ) for record in self.init])
        self.path = path
        
    

class Read_Seq():
    def __init__(self,path):
        self.index = [x for x in open(path,'r',encoding='utf-8').read().split('\n') if x!='']
        self.path = path
        
class Read_DIR():
    def __init__(self,Dirpath,kw='.'):
        self.FNames = [x for x in os.listdir(Dirpath) if kw in x]
        self.FPaths = [os.path.join(Dirpath,x) for x in self.FNames]
    def GetFileNamesAndPaths(self):
        return self.FNames,self.FPaths
        
def Read_Primers(path):
    return [(x.split('\t')[0],x.split('\t')[1],RevComp(x.split('\t')[2]))for x in open(path,'r').read().split('\n') if x!= '']

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