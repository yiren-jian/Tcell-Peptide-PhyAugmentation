## **PhyAugmentation by HDOCK for T-cell Receptors and Peptides**

This repo covers implementation of running HDOCK to compute pseudo-labels used in **T-Cell Receptor-Peptide Interaction Prediction with Physical Model Augmented Pseudo-Labeling (KDD 2022)** by Yiren Jian, Erik Kruus and Martin Renqiang Min. 

You will also need standalone [hdock](http://hdock.phys.hust.edu.cn/) program in main directory.

TL,DR: `python main.py 10000` and results are stores in `results.txt`.


### Buidling  3D structure for a fasta sequence

#### Requirements

- Install Anaconda and create a environment
```bash
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh
```
Remember to answer `yes` for appending the PATH to `.bashrc`. Once Anaconda is installed, create a environment for the project and activate it.
```bash
conda create -n modeller python=3.7
conda activate modeller
```

- Install a standalone blast+ package and blastdb Database
The official guide can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz  
tar zxvpf ncbi-blast-2.12.0+-x64-linux.tar.gz
export PATH=$PATH:$HOME/ncbi-blast-2.12.0+/bin
```
Create a blast database
```bash
mkdir $HOME/blastdb
export BLASTDB=$HOME/blastdb
```
Download the PDB amino acids database
```bash
cd $HOME/blastdb
perl $HOME/ncbi-blast-2.12.0+/bin/update_blastdb.pl --passive --decompress pdbaa
```

- Install Biopython as the interface for blast+
```bash
pip install biopython
```

- Install MODELLER
```bash
conda config --add channels salilab
conda install modeller
```

#### Structure building
`build.py` is the script for building structures.

```python
import wget
import argparse
import os

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

from modeller import *
from modeller.automodel import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_sequence', type=str, default='TCR.fasta')
    args = parser.parse_args()

    # local blastp
    args.blast_xml = args.fasta_sequence[:-6] + '.xml'
    fasta_string = open(args.fasta_sequence).read()
    blastx_cline = NcbiblastpCommandline(query=args.fasta_sequence, db='pdbaa', outfmt=5, out=args.blast_xml)
    blastx_cline()

    result_handle = open(args.blast_xml)
    modeller_templates = []
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print('****Alignment****')
                print('sequence:', alignment.title[4:10])
                modeller_templates.append(alignment.title[4:10])


    modeller_template = modeller_templates[0]
    pdb_id = modeller_template[:4]
    pdb_chain = modeller_template[-1]

    # Download PDB file
    pdb_url = 'http://www.pdb.org/pdb/files/%s.pdb'%pdb_id
    pdb_name = wget.download(pdb_url)
    pdb_id = pdb_id.lower()

    # MODELLLER align2d.py
    with open(args.fasta_sequence.replace('fasta', 'ali'), 'w') as f:
        f.write('>P1;%s\n'%args.fasta_sequence[:-6])
        f.write('sequence:%s:::::::0.00: 0.00\n'%args.fasta_sequence[:-6])
        f.write(fasta_string)
        f.write('*\n')

    env = Environ()
    aln = Alignment(env)
    mdl = Model(env, file=pdb_id, model_segment=('FIRST:%s'%pdb_chain,'LAST:%s'%pdb_chain))
    aln.append_model(mdl, align_codes='%s%s'%(pdb_id, pdb_chain), atom_files='%s.pdb'%pdb_id)
    aln.append(file=args.fasta_sequence.replace('fasta', 'ali'), align_codes=args.fasta_sequence[:-6])
    aln.align2d(max_gap_length=50)
    aln.write(file='%s-%s.ali'%(args.fasta_sequence[:-6], pdb_id+pdb_chain), alignment_format='PIR')
    aln.write(file='%s-%s.pap'%(args.fasta_sequence[:-6], pdb_id+pdb_chain), alignment_format='PAP')

    # MODELLER model-single.py
    a = AutoModel(env, alnfile='%s-%s.ali'%(args.fasta_sequence[:-6], pdb_id+pdb_chain),
                  knowns=pdb_id+pdb_chain, sequence=args.fasta_sequence[:-6],
                  assess_methods=(assess.DOPE,
                                  assess.GA341))
    a.starting_model = 1
    a.ending_model = 1
    a.make()

    ### Clean the directory
    for item in os.listdir():
        if item.startswith(args.fasta_sequence[:-6]) and item.endswith('.pdb'):
            pass  # keep the built 3D structure
        elif item.endswith('.fasta') or item.endswith('.py'):
            pass  # keep the original fasta sequence
        else:
            os.remove(item)


if __name__ == '__main__':
    main()
```

`build.py` takes a fasta (for example, `TCR.fasta`) sequence as input. The fasta sequence looks like the following
```
CASSQEEGGGSWGNTIYF
```

create a new working directory by `mkdir example`, put `build.py` and `TCR.fasta` in it. And run the command

```bash
python build.py --fasta_sequence TCR.fasta
```

The only output is `TCR.B99990001.pdb`.
