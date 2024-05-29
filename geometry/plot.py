

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

modules = np.load('geometry/ModulesGeometry.npy')

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
    Modules = np.load('geometry/ModulesGeometry.npy')
    Modules = Modules.tolist()
    for layer in range(5):
        L = []
        if not ((layer<27) & (layer%2 ==0)):
            plt.figure(figsize = (12,8))
            TCs = event[event['good_tc_layer']==layer][0]
            #plt.scatter(TCs['good_tc_x']*10,TCs['good_tc_y']*10)
            modules = Modules[layer]
            for module_idx in range(len(modules)):
                module = modules[module_idx]
                a = 0 
                for i in range(6):
                    if module[0][i]!=0 or module[1][i]!=0:
                        a+=1
                if a != 0:
                    plt.plot(module[0][0:a] + [module[0][0]],module[1][0:a] + [module[1][0]],color = 'black')
            for TC_idx in range(len(TCs['good_tc_layer'])):
                u,v = TCs['good_tc_waferu'][TC_idx],TCs['good_tc_waferv'][TC_idx]
                #if TCs[TC_idx]['good_tc_waferu'] :
                if not (u,v) in L:
                    plt.annotate('('+str(TCs['good_tc_waferu'][TC_idx])+','+str(TCs['good_tc_waferv'][TC_idx])+')',(TCs['good_tc_x'][TC_idx]*10,TCs['good_tc_y'][TC_idx]*10))
                    L.append((u,v))
            plt.savefig('geometry/plot_geometry/'+ 'Layer' +str(layer)+'.png')
  

