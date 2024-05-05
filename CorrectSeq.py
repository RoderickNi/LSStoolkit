# CorrectSeq.py
# need muscle.exe or 

import os
from IO import Read_Fasta,Read_DIR
from multiprocessing import Pool
from collections import Counter
import time

def get_max_key(dictionary):
    return max(dictionary, key=dictionary.get)

def Muscle(filename,filepath,OS):
    print(f'Perform MSA on {filename} by Muscle')
    if OS == 'win':
        out_aln=filepath.replace('_PreAlnTemp.fasta','_AlnedTemp.fasta')
        os.system(f'''
            MuscleWin.exe -in {filepath} -out {out_aln} -maxiters 100 -quiet
        ''')
    elif OS == 'linux':
        out_aln=filepath.replace('_PreAlnTemp.fasta','_AlnedTemp.fasta')
        os.system(f'''
            muscle -in {filepath} -out {out_aln} -maxiters 100 -quiet
        ''')
    print(f'{filename} Done...')


def AligmentCombine(filename,filepath):
    AlleleName=filename.split('_')[0]
    NameLst,SeqLst=Read_Fasta(filepath).Get_NameLst_SeqLst()  ########
    Total_Num=NameLst[0].split('::')[1]
    SEQUENCE=''
    for i in range(len(SeqLst[0])):
        Base_lst=[]
        for seq in SeqLst:
            Base_lst.append(seq[i])
        SEQUENCE+=get_max_key(Counter(Base_lst))
    return '>'+AlleleName+'::'+Total_Num+'\n'+SEQUENCE.replace('-', '')

def Correct(fasta_path,Output,OS='win',cpu_num=1):
    # 数据载入与准备
    Names_Seqs=Read_Fasta(fasta_path).GetNameSeq_to_lst()
    LST=[]
    for name_seq in Names_Seqs:
        filename=name_seq[0].split('@')[0]+'_PreAlnTemp.fasta'
        seqname=name_seq[0].split('@')[1]
        seq=name_seq[1]
        File_path=os.path.join(os.path.dirname(fasta_path),filename)
        with open(File_path,'a',encoding='utf-8') as f:
            print(f'>{seqname}\n{seq}',file = f)
    # MSA
    FileNames,FilePaths=Read_DIR(os.path.dirname(fasta_path),kw='_PreAlnTemp.fasta').GetFileNamesAndPaths()
    Inputs=[(filename,filepath,OS) for filename,filepath in zip(FileNames,FilePaths)]
    pool=Pool(processes=cpu_num)
    pool.starmap(Muscle,Inputs)
    pool.close()
    pool.join()
    time.sleep(1)
    # 序列矫正
    FileNames,FilePaths=Read_DIR(os.path.dirname(fasta_path),kw='_AlnedTemp.fasta').GetFileNamesAndPaths()
    Inputs=[(filename,filepath) for filename,filepath in zip(FileNames,FilePaths)]
    pool=Pool(processes=cpu_num)
    results=pool.starmap(AligmentCombine,Inputs)
    results=sorted(results,key=lambda x:int(x.split('\n')[0].split('-')[1].split('::')[0]))
    pool.close()
    pool.join()
    # 输出结果

    with open(Output,'w',encoding='utf-8') as f:
        for rst in results:
            print(rst,file=f)

    time.sleep(1)
    #删除临时文件
    if OS == 'win':            
        os.system(f'''
            del {os.path.join(os.path.dirname(fasta_path),'*Temp.fasta')}
        ''')
    elif OS == 'linux':
        os.system(f'''
            rm -rf {os.path.join(os.path.dirname(fasta_path),'*Temp.fasta')}
        ''')    

