import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--ha_tree', type=str, required=True)
parser.add_argument('--ha_alignment', type=str, required=True)
parser.add_argument('--nonha_trees', type=str, required=True, nargs='+')
parser.add_argument('--nonha_alignments', type=str, required=True, nargs='+')
parser.add_argument('--segments', type=str, required=True, nargs='+')
parser.add_argument('--descriptor', type=str, required=True)

args = parser.parse_args()

with open(args.descriptor, 'w') as f:
    f.write(f"*HA,../../../{args.ha_alignment},../../../{args.ha_tree}\n")
    for seg,aln,tree in zip(args.segments, args.nonha_alignments, args.nonha_trees):
        f.write(f"{seg},../../../{aln},../../../{tree}\n")
