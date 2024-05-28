import numpy as np
import awkward as ak
import uproot
import math
import yaml
import os
import argparse
import numpy as np


def read_xml_plot(Edges):
    if Edges == 'yes':
        tree = ET.parse('config_files/AllocationpTTsEdges.xml')
    if Edges == 'no':
        tree = ET.parse('config_files/AllocationpTTsNoEdges.xml')
        
    root = tree.getroot()
    reversed_data_pTT  = defaultdict(list)

    S1_index = 0
    for s1_element in root.findall('.//S1'):
        for channel_element in s1_element.findall('.//Channel'):
            channel = int(channel_element.get('aux-id'))
            for frame_element in channel_element.findall('.//Frame'):
                if all(attr in frame_element.attrib for attr in ['id','pTT']):
                    frame  = int(frame_element.get('id'))
                    pTT     = hex(int(frame_element.get('pTT'),16))
                    n_link = 14 + 14*math.floor(channel/2) + S1_index
                    reversed_data_pTT[(frame,n_link,channel)].append({'S1board': S1board, 'eta': eta ,'phi': phi })

        S1_index += 1
    return reversed_data_pTT
