#!/usr/bin/python
import os
import sys
import re
from Bio import SeqIO

genome = sys.argv[1]
site = sys.argv[2]
output_dir = sys.argv[3]
primer = sys.argv[4]

with open(f'{output_dir}/primer.fasta', 'w') as primer_fasta:
    primer_fasta.write(f'>primer\n{primer}\n')

def check_restriction_sites(sequence, enzyme_site):
    enzyme_site_tmp = enzyme_site.split('_')
    enzyme_site = enzyme_site_tmp[0]+enzyme_site_tmp[1]
    enzyme_pattern = r"{0}".format(enzyme_site)
    sites = re.finditer(enzyme_pattern, sequence)
    return [match.start()+len(enzyme_site_tmp[0]) for match in sites]

def generate_restriction_fragments(sequence, founded_sites):
    founded_sites.sort()
    fragments = []
    for start in range(len(founded_sites) - 1):
        fragment_start = founded_sites[start] 
        fragment_end = founded_sites[start + 1]
        fragment = sequence[fragment_start:fragment_end]
        fragments.append(fragment)
    return fragments

def find_bends(fragments, blastn_out):
    with open(fragments, 'r') as fragm, \
    open(blastn_out, 'r') as bl_out, \
    open(f'{output_dir}/bends_size.txt', 'w') as bends:
        dict_primer_sites = {}
        list_bends = []
        for line in bl_out:
            line = line.split('\t')
            if line[1] not in dict_primer_sites:
                dict_primer_sites[line[1]] = []
            if line[1] in dict_primer_sites:
                dict_primer_sites[line[1]].append(line[8])
        for seq in SeqIO.parse(fragm, 'fasta'):
            if seq.id in dict_primer_sites:
                lenght = len(seq.seq)
                for primer_site in dict_primer_sites[seq.id]:
                    bent_lenght = lenght - int(primer_site)
                    list_bends.append(bent_lenght)
        list_bends.sort()
        for bend in list_bends:
            bends.write(f'{bend}\n')

with open(genome, 'r') as fasta, \
open(f'{output_dir}/restriction_fragments.fasta', 'w') as output:
    dict_sites = {}
    dict_fragments = {}
    counter = 0
    for seq in SeqIO.parse(fasta, 'fasta'):
        dict_sites[seq.id] = check_restriction_sites(str(seq.seq), site)
        dict_sites[seq.id].append(0)
        dict_sites[seq.id].append(len(seq.seq))
        dict_fragments[seq.id] = generate_restriction_fragments(str(seq.seq), dict_sites[seq.id])
    for chromosome in dict_fragments:
        for fragment in dict_fragments[chromosome]:
            counter += 1
            output.write(f'>restriction_fragment_{counter}\n{fragment}\n')

os.system('makeblastdb -in {0} -dbtype nucl'.format(f'{output_dir}/restriction_fragments.fasta'))
os.system('blastn -db {0} -max_target_seqs 30000 \
-query {1} \
-evalue 0.1 -task blastn -outfmt 6 -word_size 6 \
-out {2}'.format(f'{output_dir}/restriction_fragments.fasta', f'{output_dir}/primer.fasta', f'{output_dir}/results.xml'))

find_bends(f'{output_dir}/restriction_fragments.fasta', f'{output_dir}/results.xml')

#os.remove(f'{output_dir}/results.xml')
os.remove(f'{output_dir}/primer.fasta')



