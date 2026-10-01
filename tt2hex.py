#!/usr/bin/env python3
# Usage : python3 tt2hex.py entree.txt > liste.txt
# Entree : une fonction par ligne = 64 valeurs 0/1 (separees par des espaces, ou collees),
#          f(0) en premier. Un prefixe du type "TT=" est ignore. Lignes vides et # ignorees.
import sys, re
for n, line in enumerate(open(sys.argv[1]), 1):
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    bits = re.findall(r'[01]', line.split('=', 1)[-1])
    if len(bits) != 64:
        sys.stderr.write(f"ligne {n}: {len(bits)} valeurs au lieu de 64, ignoree\n")
        continue
    f = sum(int(b) << i for i, b in enumerate(bits))
    if bin(f).count('1') % 2 == 0:
        sys.stderr.write(f"ligne {n}: poids pair, ignoree\n")
        continue
    print('%016x' % f)
