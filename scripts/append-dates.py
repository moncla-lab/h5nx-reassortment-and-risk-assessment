import argparse
import pandas as pd
from Bio import SeqIO

def extend_wildcards(f_list, wildcard, genes, verbose=None):
    f_list_copy = f_list[:]
    for f in f_list:
        if wildcard in f:
            f_list_copy.extend([f.replace(args.wildcard, gene) for gene in args.genes])
            f_list_copy.remove(f)
    if verbose and f_list != f_list_copy:
        print(f"Replaced wildcard(s) in {verbose} FASTA paths with: {', '.join(genes)}")
    return f_list_copy


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--metadata', action='store', type=str, dest='metadata',
                        required=True, help='Path to input metadata tsv')
    parser.add_argument('--fastas', action='store', type=str, nargs='+', dest='fastas',
                        required=True, help='Path to input FASTAs, separated by a comma')
    parser.add_argument('--output_metadata', action='store', type=str, dest='output_metadata',
                        required=True, help='Path to output metadata tsv with dates appended to strain names')
    parser.add_argument('--output_fastas', action='store', type=str, nargs='+', dest='output_fastas',
                        required=True, help='Path to output FASTAs with dates appended to strain names, separated by a comma')
    parser.add_argument('--wildcard', action='store', type=str, dest='wildcard', default='*',
                        required=False, help='Wildcard symbol used in place of gene names in FASTA paths')
    parser.add_argument('--genes', action='store', type=str, nargs='+', dest='genes',
                        default=['pb2', 'pb1', 'pa', 'ha', 'np', 'na', 'mp', 'ns'],
                        required=False, help='List of genes to fill in place of wildcard in FASTA paths')

    args = parser.parse_args()

    print(args.fastas)
    
    # expand any wildcard symbols in input and output FASTAs
    args.fastas = extend_wildcards(args.fastas, args.wildcard, args.genes, verbose='input')
    args.output_fastas = extend_wildcards(args.output_fastas, args.wildcard, args.genes, verbose='output')

    if len(args.fastas) != len(args.output_fastas):
        print(f"WARNING: Specified {len(args.fastas)} input FASTAs, but {len(args.output_fastas)} output FASTAs")
        

    df = pd.read_csv(args.metadata, sep='\t')
    orig_strains = list(df['strain'])
    df['strain'] = df['strain'] + '|' + df['date']

    df.to_csv(args.output_metadata, sep='\t', index=False)
    
    mapping = dict(zip(orig_strains, df['strain']))

    for fasta, output_fasta in zip(args.fastas, args.output_fastas):
        records = []
        for record in SeqIO.parse(fasta, 'fasta'):
            record.id = mapping[record.id]
            record.description = record.id
            records.append(record)
        SeqIO.write(records, output_fasta, 'fasta')