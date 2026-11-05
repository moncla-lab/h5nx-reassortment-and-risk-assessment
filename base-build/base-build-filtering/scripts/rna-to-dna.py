import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--alignments', type=str, required=True, nargs='+', help='path to alignment fasta files for each segment')
parser.add_argument('--dna_alignments', type=str, required=True, nargs='+', help='path to output alignment fasta files for each segment, with rna converted to dna')
args = parser.parse_args()

for alignment, dna_alignment in zip(args.alignments, args.dna_alignments):
    with open(alignment) as f:
        t = f.read()
        lines = []
        for l in t.split('\n'):
            if not l.startswith('>'):
                # replace U bases with T
                l = l.replace('u', 't')
                l = l.replace('U', 'T')
            lines.append(l)
    with open(dna_alignment, 'w') as f:
        f.write('\n'.join(lines))