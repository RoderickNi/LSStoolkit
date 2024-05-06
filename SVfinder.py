import os
import sys


if __name__=='__main__':
    
    # default:
    PF='map-ont'
    SVlen=50
    CPU=3
    # argv:
    for i in range(len(sys.argv)):
        j=i+1
        if '--Platform' == sys.argv[i]:  # map-pb or map-ont
            PF=str(sys.argv[j])
        elif '--CPU'  == sys.argv[i]:
            CPU=str(sys.argv[j])
        elif '--REF'  == sys.argv[i]:
            REF=str(sys.argv[j])
        elif '--Fasta' == sys.argv[i]:
            Fasta = str(sys.argv[j])
        elif '--SVlen' == sys.argv[i]:
            SVlen = str(sys.argv[j])
        elif '--OUT'  == sys.argv[i]:
            Out = str(sys.argv[j])
    
    
    print('=======================================================')
    print('        Map to Reference Sequence (by minimap2)        ')
    print('=======================================================')
    os.system(f'''
        minimap2 --MD -ax {PF} -t {CPU} {REF} {Fasta} > mapped.sam
              ''')
    
    os.system(f'''
        samtools view -bS mapped.sam | samtools sort -m 10G -o mapped_sorted.bam
        samtools index -b mapped_sorted.bam
              ''')
    print('=======================================================')
    print('                SV Calling (by sniffles)               ')
    print('=======================================================')
    os.system(f'''
        sniffles -n -1 -s 1 -t {CPU} -l {SVlen}  -m mapped_sorted.bam -v {Out} 
              ''')
    