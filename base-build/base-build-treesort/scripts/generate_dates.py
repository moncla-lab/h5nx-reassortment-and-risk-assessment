import argparse
from datetime import datetime

parser = argparse.ArgumentParser()

parser.add_argument('--alignment', type=str, required=True)
parser.add_argument('--dates', type=str, required=True)

args = parser.parse_args()

def dec_date(date):
    return date.year + ((date.month - 1) * 30 + date.day) / 365.0

with open(args.alignment) as f:
    lines = f.read().split('\n')
    strains = [l[1:] for l in lines if l.startswith('>')]
    dates = [s.split('|')[1] for s in strains]
    dates = [d.upper().replace('XX','01') for d in dates]
    dates = [datetime.strptime(d, '%Y-%m-%d') for d in dates]
    dates = [dec_date(d) for d in dates]

with open(args.dates, 'w') as f:
    f.write('name,date')
    for s,d in zip(strains,dates):
        f.write('\n')
        f.write(f'{s},{d}')