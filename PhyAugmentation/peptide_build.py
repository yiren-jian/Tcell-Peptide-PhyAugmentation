import wget
import argparse
import os
import pickle

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

from modeller import *
from modeller.automodel import *
import random
import numpy as np
import shutil

def main():
    peps = pickle.load( open('pep_database.pickle', 'rb' ) )

    for pep in peps:
        amino_acids = [letter for letter in 'ARNDCEQGHILKMFPSTWYV']
        if len(pep) < 14:
            # pep = pep + ''.join(np.random.choice(amino_acids, size=10, replace=True))
            pass

        build_history = []
        build_history_str = []

        fasta_string = pep

        with open('my.fasta', 'w') as f:
            f.write(fasta_string)
        blastx_cline = NcbiblastpCommandline(query='my.fasta', task='blastp-short', word_size=2, db='pdbaa', outfmt=5, out='my.xml')
        blastx_cline()

        result_handle = open('my.xml')
        modeller_templates = []
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    print('****Alignment****')
                    print('sequence:', alignment.title[4:10])
                    modeller_templates.append(alignment.title[4:10])

        try:
            modeller_template = modeller_templates[0]
            pdb_id = modeller_template[:4]
            pdb_chain = modeller_template[-1]

            # Download PDB file
            pdb_url = 'http://www.pdb.org/pdb/files/%s.pdb'%pdb_id
            pdb_name = wget.download(pdb_url)
            pdb_id = pdb_id.lower()

            # MODELLLER align2d.py
            with open('my.ali', 'w') as f:
                f.write('>P1;my\n')
                f.write('sequence:my:::::::0.00: 0.00\n')
                f.write(fasta_string)
                f.write('*\n')

            env = Environ()
            aln = Alignment(env)
            mdl = Model(env, file=pdb_id, model_segment=('FIRST:%s'%pdb_chain,'LAST:%s'%pdb_chain))
            aln.append_model(mdl, align_codes='%s%s'%(pdb_id, pdb_chain), atom_files='%s.pdb'%pdb_id)
            aln.append(file='my.ali', align_codes='my')
            aln.align2d(max_gap_length=50)
            aln.write(file='my-%s.ali'%(pdb_id+pdb_chain), alignment_format='PIR')
            aln.write(file='my-%s.pap'%(pdb_id+pdb_chain), alignment_format='PAP')

            # MODELLER model-single.py
            a = AutoModel(env, alnfile='my-%s.ali'%(pdb_id+pdb_chain),
                          knowns=pdb_id+pdb_chain, sequence='my',
                          assess_methods=(assess.DOPE,
                                          assess.GA341))
            a.starting_model = 1
            a.ending_model = 1
            a.make()
            os.rename('my.B99990001.pdb', '%s.pdb'%pep)
            shutil.move('%s.pdb'%pep, '../PDB-for-Peptides/%s.pdb'%pep)
            build_history.append(True)
            build_history_str.append(pep)
            build_history_str.append('True')
        except:
            build_history.append(False)
            build_history_str.append(pep)
            build_history_str.append('False')


        ##### clean the directory
        skip_list = ['build.py', 'hdock', 'tcr_unknown.pickle', 'pep_database.pickle', 'main.py',
                     'clean.py', 'results.txt', 'peptide_build.py', 'hdock_example.py',
                     'README.md']
        for item in os.listdir():
            if item in skip_list:
                pass
            else:
                os.remove(item)

        with open('results.txt', 'a') as f:
            f.write(' '.join(build_history_str))
            f.write('\n')



if __name__ == '__main__':
    main()
