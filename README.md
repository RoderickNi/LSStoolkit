# LSStoolkit
It is a multifunctional module toolkit  for processing long fragment single-molecule sequencing clean data of target genes or segments.

## Installation
- Conda environment    
```
conda install muscle=3.8
conda install python=3.10
conda install numpy scipy scikit-learn numba
pip install levenshtein
pip install biopython==1.81
pip install umap-learn
pip install tqdm
```
- Get LSStoolkit
```
git clone https://github.com/RoderickNi/LSStoolkit.git
```
## Usage
- Separator
```
 python Separator.py --Fastq TestFastq.fq --Primer TestPrimers.txt --OutDir TestRst --MinLenth 5000 --MisMatch 1 --TagTrim 6 --CPU 8
```
![image execution flow](https://github.com/RoderickNi/LSStoolkit/blob/main/Separator.png)
- HapGrep
```
 python HapGrep.py --Fasta ./TestFasta.fasta --Out ./TestRepReads.fasta --Pkl ./TestRst.pkl --CPU 8
```
![image execution flow](https://github.com/RoderickNi/LSStoolkit/blob/main/HapGrep.png)
- SVfinder
```
 python SVfinder.py ...
```

