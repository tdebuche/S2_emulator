

import numpy as np
import awkward as ak
import uproot
import math
import yaml
import os
import argparse
import numpy as np
import xml.etree.ElementTree as ET
from collections import defaultdict
import matplotlib.pyplot as plt

            
#######################################################################################
############################### PROVIDE EVENTS ########################################
#######################################################################################

with open('config.yaml', "r") as afile:
    cfg_particles = yaml.safe_load(afile)["particles"]


def provide_events(n, particles, PU):
    base_path = cfg_particles['base_path']
    name_tree = cfg_particles[PU][particles]["tree"]
    filepath  = base_path + cfg_particles[PU][particles]["file"]

    branches_tc = [
        'good_tc_x', 'good_tc_y', 'good_tc_z',
        'good_tc_phi', 'good_tc_layer', 'good_tc_cellv',
        'good_tc_waferu', 'good_tc_waferv',
        'good_tc_pt', 'good_tc_subdet'
    ]

    tree = uproot.open(filepath)[name_tree]
    events_ds = []
    for ev in range(n):
      data = tree.arrays(branches_tc, entry_start=ev, entry_stop=ev+1, library='ak')
      events_ds.append(data)
    return events_ds

def plot_uv(event):
    for layer in range(47):
        if not ((layer)<27 & (layer%2 ==0)):
            plt.figure(figsize = (12,8))
            TCs = event[event['good_tc_layer']==layer][0]
            plt.scatter(TCs['good_tc_x'],TCs['good_tc_y'])
            for TC_idx in range(len(TCs)):
                if TCs[TC_idx]['good_tc_waferu'] :
                    plt.annotate('('+str(TCs[TC_idx]['good_tc_waferu'])+','+str(TCs[TC_idx]['good_tc_waferv'])+')',(TCs[TC_idx]['good_tc_x'],TCs[TC_idx]['good_tc_y']))
            plt.savefig('geometry/plot_geometry/'+ 'Layer' +str(layer)+'.png')
  

