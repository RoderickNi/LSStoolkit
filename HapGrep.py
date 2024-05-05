import os
import sys
from CorrectSeq import Correct

if __name__ == '__main__':
    
    # default:
    size=50
    otpkl=None
    kmin=1
    kmax=5
    radius=2
    num_closest=10
    CPU=3
    # argv:
    for i in range(len(sys.argv)):
        mark=sys.argv[i]
        if '--Fasta' == mark:
            path = str(sys.argv[i+1])
        elif '--Out' == mark:
            otpt = str(sys.argv[i+1])
        elif '--Size' == mark:
            size = str(sys.argv[i+1])
        elif '--Pkl' == mark:
            otpkl = str(sys.argv[i+1])
        elif '--RepReadsNum' == mark:
            num_closest = str(sys.argv[i+1])
        elif '--Radius' == mark:
            radius = str(sys.argv[i+1])
        elif '--KmerMin' == mark:
            kmin = str(sys.argv[i+1])
        elif '--KmerMax' == mark:
            kmax = str(sys.argv[i+1])
        elif '--CPU' == mark:
            CPU = str(sys.argv[i+1])

            
    print('=======================================================')
    print('Reads deduplication and representative reads extraction')
    print('=======================================================')
    os.system(f'''
            python ReduceAndRepSeq.py -fa {path} -out {otpt} -size {size} -pkl {otpkl} -otnum {num_closest} -radius {radius} -kmin {kmin} -kmax {kmax}
              ''')
    print('=======================================================')
    print('Reads deduplication and representative reads extraction')
    print('=======================================================')
    HapOut=otpt+'.Haplotypes.fasta'
    print(otpt)
    Correct(otpt,HapOut,OS='linux',cpu_num=int(CPU))
    
