#!/usr/bin/env python

if __name__ == "__main__":
    import sys


    if len(sys.argv) < 2:
        print('USAGE: python nwk2mat.py TREE.nwk OUTDIR')
        sys.exit(1)
    
    import pandas as pd
    import itertools
    from Bio import Phylo
    from pathlib import Path
    import re
    import os

    ifile = sys.argv[1]
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path('.')
    
    t = Phylo.read(ifile, 'newick')

    d = {}
    for x, y in itertools.combinations(t.get_terminals(), 2):
        v = t.distance(x, y)
        d[x.name] = d.get(x.name, {})
        d[x.name][y.name] = v
        d[y.name] = d.get(y.name, {})
        d[y.name][x.name] = v
    for x in t.get_terminals():
        d[x.name][x.name] = 0

    m = pd.DataFrame(d)
    m.to_csv(str(output_dir / "tree.matrix.tsv"), sep='\t', index=False) 
    for col in m.columns:
        if re.search(r'str', col):
            series = m[col].sort_values()
            
            # Output column series to file
            with open(str(output_dir / (col + ".sorted_distance.txt")), 'w') as f:
                for index, value in series.items():
                    f.write(index + '\t' + str(value) + '\n')