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

    modeller_template = modeller_templates[0]   ### 1O9J|A
    pdb_id = modeller_template[:4]    ###  1O9J
    pdb_chain = modeller_template[-1]   ###  A
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
    ### Clean the directory


if __name__ == '__main__':
    main()
