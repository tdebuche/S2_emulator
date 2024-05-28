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


def read_xml_plot(Edges):
    if Edges == 'yes':
        tree = ET.parse('config_files/AllocationpTTsEdges.xml')
    if Edges == 'no':
        tree = ET.parse('config_files/AllocationpTTsNoEdges.xml')
        
    root = tree.getroot()
    data_pTT  = defaultdict(list)

    S1_index = 0
    for s1_element in root.findall('.//S1'):
        for channel_element in s1_element.findall('.//Channel'):
            channel = int(channel_element.get('aux-id'))
            for frame_element in channel_element.findall('.//Frame'):
                if all(attr in frame_element.attrib for attr in ['id','pTT']):
                    frame  = int(frame_element.get('id'))
                    pTT     = frame_element.get('pTT')
                    n_link = 14 + 14*math.floor(channel/2) + S1_index
                    S1Board,eta,phi,CEECEH = get_pTT_numbers(pTT)
                    data_pTT[(S1Board,eta,phi)].append((frame,n_link,channel))

        S1_index += 1
    return data_pTT


def create_energies(data_links,args):
    Edges = args.Edges
    etaphi_links = read_xml_plot(Edges)
    print(etaphi_links)
    if Edges == 'yes': nb_phi = 28
    else : nb_phi = 24
    energies = [[0 for phi in range(nb_phi)]for eta in range(20)]
    for S1Board in range(14):
        for eta in range(20):
            for phi in range(nb_phi):
                if data_links[etaphi_links[(S1Board,eta,phi)][0]] != []:
                    energies[eta][phi] += data_links[etaphi_links[(S1Board,eta,phi)][0]][0]
    return energies


def get_pTT_numbers(pTT):
    S1Board = int(pTT[4:6],16) & 0x3F
    phi = int(pTT,16) & 0x1F
    eta = int(pTT,16) & 0x3E0 //(16 * 2)
    CEECEH = int(pTT,16) & 0x400 //(16*16*4)
    return(S1Board,eta,phi,CEECEH)
    
