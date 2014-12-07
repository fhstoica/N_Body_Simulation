#! /usr/bin/python

dim        = 3
tot_points = 16
vscale     = 4.0
pscale     = 4.0

import sys
import os
import random

OUTFILE      = 'initial_data.dat'
r_seed       = float(sys.argv[1])

random.seed(r_seed)
f  = open(OUTFILE , 'w')
for a in range(tot_points):
    for b in range(dim):
        f.write(str("%1.6f" % (pscale*(random.random() - 0.5)))+' ')
    for b in range(dim):
        f.write(str("%1.6f" % (vscale*(random.random() - 0.5)))+' ')
    f.write('\n')
    
f.close()

