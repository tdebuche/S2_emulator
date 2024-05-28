import numpy as np
import awkward as ak
import uproot
import math
import yaml
import os
import argparse
import cppyy
import numpy as np
cppyy.add_include_path(os.path.abspath(''))
cppyy.load_library('lib_configuration.so')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h')

from cppyy.gbl import l1thgcfirmware
from data_handle.S1simulator import build_pTTsCEE
from data_handle.geometry import read_txt_pTTs
from data_handle.geometry import read_xml_pTTs
from data_handle.event_pTT import provide_events

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

S1pTTCEE,S1pTTCEH = read_txt_pTTs(args.Edges,args.Sector)
xml_allocation = read_xml_pTTs(args.Edges)

events = provide_events(args.n, args.particles, args.pileup)

events._pTT_packer(args, xml_allocation,S1pTTCEE,S1pTTCEH)

print(events.pTT_packer)