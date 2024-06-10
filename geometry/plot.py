

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
    for layer in range(1,48):
        L = []
        if not ((layer<27) & (layer%2 ==0)):
            plt.figure(figsize = (12,8))
            TCs = event[event['good_tc_layer']==layer][0]
            #plt.scatter(TCs['good_tc_x']*10,TCs['good_tc_y']*10)
            modules = Modules[layer-1]
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
                    L.append((u,v))
                    u,v,sector = getuvsector(layer,u,v)
                    plt.annotate('('+str(u)+','+str(v)+')',(TCs['good_tc_x'][TC_idx]*10,TCs['good_tc_y'][TC_idx]*10))
                    plt.scatter(TCs['good_tc_x'][TC_idx]*10,TCs['good_tc_y'][TC_idx]*10)
            plt.savefig('geometry/plot_geometry/'+ 'Layer' +str(layer)+'.png')

def Sector0(layer,u,v):
        if (layer <34) and (layer != 30) and (layer != 32) and (layer != 28):
            if (v-u > 0) and (v >= 0):
                return(True)
        if (layer >= 34) and (layer%2 == 0):
            if (v-u > 0) and (v > 0):
                return(True)
        if (layer >= 34) and (layer%2 == 1):
            if (v-u >= 0) and (v >= 0):
                return(True)
        if (layer == 28) or (layer == 30) or (layer == 32):
            if (u - 2*v <0) and (u+v >= 0):
                return(True)
        return False

def getuvsector(layer,u,v):
        if u == -999:
            return (u,v,0)
        if Sector0(layer,u,v):
            if (layer != 28) and (layer != 30) and (layer != 32): 
                return(v-u,v,0)
            else :
                if u >= 0:
                    return (v,u,1)
                else :
                    return(-u,v-u,1)
        else:
            if  (layer <34):
                u,v = -v,v-u
            if (layer >= 34) and (layer%2 == 0):
                u,v = -v+1,u-v+1
            if (layer >= 34) and (layer%2 == 1):
                u,v = -v-1,u-v-1
            if Sector0(layer,u,v):
                if (layer != 28) and (layer != 30) and (layer != 32): 
                    return(v-u,v,1)
                else:
                    if u >= 0:
                        return (v,u,1)
                    else :
                        return(-u,v-u,1)
                    
            else : 
                if  (layer <34):
                    u,v = -v,v-u
                if (layer >= 34) and (layer%2 == 0):
                    u,v = -v+1,u-v+1
                if (layer >= 34) and (layer%2 == 1):
                    u,v = -v-1,u-v-1
                if Sector0(layer,u,v):
                    if (layer != 28) and (layer != 30) and (layer != 32): 
                        return(v-u,v,2)
                    else :
                        if u >= 0:
                            return (v,u,1)
                        else :
                            return(-u,v-u,1)
