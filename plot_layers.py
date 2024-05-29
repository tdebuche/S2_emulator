from geometry.plot import plot_uv
from geometry.plot import provide_events
import numpy as np
import awkward as ak
import uproot
import math
import yaml
import os
import argparse
import cppyy
import numpy as np

parser = argparse.ArgumentParser(description='Stage-2 Emulator Parameters')
parser.add_argument('-n',          type=int, default=1,         help='Provide the number of events')
parser.add_argument('--particles', type=str, default='photons', help='Choose the particle sample')
parser.add_argument('--pileup',    type=str, default='PU0',     help='Choose the pileup - PU0 or PU200')
parser.add_argument('--plot',        action='store_true', help='Create plots')
parser.add_argument('--col',         action='store_true', help='Create plots using column numbers')
parser.add_argument('--phi',         action='store_true', help='Create plots using phi coordinates')
parser.add_argument('--performance', action='store_true', help='Create performance plots: distance gen_particle/max_TC')
parser.add_argument('--thr_seed',    action='store_true', help='Create efficiency plots post seeding')
parser.add_argument('--cl_energy',   action='store_true', help='Create plot of gen_pt vs recontructed energy')
parser.add_argument('--Sector',      type=int, default=0, help='Sector of S2 Board')
parser.add_argument('--Edges',   default = 'no', help='20*24 or 20*28 bins')

args = parser.parse_args()

events = provide_events(args.n, args.particles, args.pileup)
for idx, event in enumerate(events):
    print(event[0]['good_tc_layer'])
    plot_uv(event)
  

