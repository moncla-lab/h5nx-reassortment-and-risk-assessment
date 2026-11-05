import re
import dendropy
import pandas as pd
import baltic as bt
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--treesort_tree', type=str, required=True, help='path to treesort .tre file')
parser.add_argument('--newick_tree', type=str, required=True, help='path to refined tree .nwk file')
parser.add_argument('--cleaned_treesort_tree', type=str, required=True, help='path to export a cleaned treesort .tre file')
parser.add_argument('--output', type=str, required=True, help='path to output node data json')
args = parser.parse_args()

with open(args.treesort_tree) as f:
    text = f.read()

# Function to replace all commas *inside double quotes* in [&...] blocks
def replace_commas_inside_quotes(match):
    block = match.group(0)
    # Replace all commas that are inside quotes
    # This keeps quotes and other commas outside intact
    def comma_replacer(m):
        return m.group(0).replace(',', ';')
    return re.sub(r'"[^"]*"', comma_replacer, block)

# Apply to all [&...] annotation blocks
fixed_text = re.sub(r'\[&[^\]]*\]', replace_commas_inside_quotes, text)

# Save cleaned file
with open(args.cleaned_treesort_tree, "w") as f:
    f.write(fixed_text)

# Load the tree
tree = dendropy.Tree.get(
    path=args.cleaned_treesort_tree,
    schema="nexus",
    preserve_underscores=True
)

records = []

# Iterate through all nodes (internal and terminal)
for node in tree:
    # Get descendant tip names
    descendants = [leaf.taxon.label for leaf in node.leaf_nodes() if leaf.taxon]

    # Extract metadata (annotations)
    traits = {ann.name: ann.value for ann in node.annotations}

    traits['is_reassorted'] = {'0':False, '1':True}[traits.get('is_reassorted', '0')]
    traits['rea'] = traits.get('rea', None)
    

    records.append({
        "node": node.label if node.label else None,
        "n_descendants": len(descendants),
        "descendants": descendants,
        **traits  # expand the metadata as columns
    })

# Convert to a DataFrame for easy inspection
df = pd.DataFrame(records)


tree = bt.loadNewick(args.newick_tree)
tip_mapping = {k.name:k for k in tree.getExternal()}

GENES = ["PB2", "PB1", "PA", "NP", "NA", "MP", "NS"] # no HA, can't be reassorted relative to itself

rea_df = df.query("is_reassorted")
desc = rea_df['descendants']
rea = rea_df['rea']

node_data = {}

for d,r in zip(desc, rea):
    if len(d) > 1:
        mrca = tree.commonAncestor([tip_mapping[n] for n in d])
        name = mrca.traits['label']
    else:
        name = d[0]
        tip = tip_mapping[name]
    r = r.replace(';', ', ')

    if name in node_data:
        # print(f'Adding to {name}')
        # print(node_data[name]['Reassortment Annotation'])
        # print(node_data[name]['Reassorted Segments'])
        # print(node_data[name]['Number of Reassorted Segments'])

        node_data[name]['Reassortment Annotation'] = node_data[name]['Reassortment Annotation'] + ', ' + r
        node_data[name]['Reassorted Segments'] = ', '.join([node_data[name]['Reassorted Segments']] + [gene for gene in GENES if gene in r])
        node_data[name]['Number of Reassorted Segments'] = node_data[name]['Number of Reassorted Segments'] + sum(gene in r for gene in GENES)

        # print(node_data[name]['Reassortment Annotation'])
        # print(node_data[name]['Reassorted Segments'])
        # print(node_data[name]['Number of Reassorted Segments'])

    else:
        node_data[name] = {
            'Reassorted': 'True',
            'Reassortment Annotation': r,
            'Reassorted Segments': ', '.join([gene for gene in GENES if gene in r]),
            'Number of Reassorted Segments': sum(gene in r for gene in GENES)
        }

    for gene in GENES:
        node_data[name][f'{gene} Reassorted'] = node_data[name].get(f'{gene} Reassorted') or str(gene in r)
    
        
non_rea_df = df.query("~is_reassorted")
non_rea_desc = non_rea_df['descendants']

for d in non_rea_desc:
    if len(d) > 1:
        mrca = tree.commonAncestor([tip_mapping[n] for n in d])
        name = mrca.traits['label']
    else:
        name = d[0]
        tip = tip_mapping[name]

    if name not in node_data:
        node_data[name] = {
            'Reassorted': 'False',
            'Reassortment Annotation': '',
            'Reassorted Segments': '',
            'Number of Reassorted Segments': 0
        }

        for gene in GENES:
            node_data[name][f'{gene} Reassorted'] = 'False'
    


# format for auspice node_date
out_dict = {'nodes': node_data}

# export as json
with open(args.output, 'w') as f:
    json.dump(out_dict, f)