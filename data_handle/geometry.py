import xml.etree.ElementTree as ET
from collections import defaultdict
import math
import numpy as np
import awkward as ak

def read_xml():
    tree = ET.parse('config_files/S1.ChannelAllocation.xml')
    root = tree.getroot()
    reversed_data_si  = defaultdict(list)
    reversed_data_sci = defaultdict(list)

    S1_index = 0
    for s1_element in root.findall('.//S1'):
        for channel_element in s1_element.findall('.//Channel'):
            channel = int(channel_element.get('aux-id'))
            for frame_element in channel_element.findall('.//Frame'):
                if all(attr in frame_element.attrib for attr in ['id', 'column', 'Motherboard']):
                    frame  = int(frame_element.get('id'))
                    column = int(frame_element.get('column'))
                    MB     = frame_element.get('Motherboard')
                    n_link = 14 + 14*math.floor(channel/3) + S1_index
                    index  = int(frame_element.get('index'))
                    reversed_data_sci[MB].append({'column' : column,
                                                  'frame'  : frame, 
                                                  'channel': channel, 
                                                  'n_link' : n_link,
                                                  'index'  : index})
                if all(attr in frame_element.attrib for attr in ['id', 'column', 'Module']):
                    frame  = int(frame_element.get('id'))
                    column = int(frame_element.get('column'))
                    module = hex(int(frame_element.get('Module'),16))
                    n_link = 14 + 14*math.floor(channel/3) + S1_index
                    index  = int(frame_element.get('index'))
                    reversed_data_si[module].append({'column' : column,
                                                     'frame'  : frame, 
                                                     'channel': channel, 
                                                     'n_link' : n_link,
                                                     'index'  : index})
        S1_index += 1
    return [reversed_data_si, reversed_data_sci]



def MB_geometry():
    tree = ET.parse('config_files/Geometry.xml')
    root = tree.getroot()
    result_dict = {}
    
    for plane in root.findall('.//Plane'):
        plane_id = int(plane.attrib['id'])
        if plane_id < 34: continue
        
        for motherboard in plane.findall('Motherboard'):
            motherboard_id = motherboard.attrib['id']
            
            for module in motherboard.findall('Module')[1:]:
                vertices = module.attrib.get('Vertices')
                if len(vertices.split(";")) > 4: break
                result_dict.setdefault(plane_id, {}).setdefault(int(module.attrib['v']), motherboard_id)
                break

    return result_dict


def read_xml_pTTs(Edges):
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
                    pTT     = frame_element.get('pTT')
                    n_link = 14 + 14*math.floor(channel/3) + S1_index
                    reversed_data_pTT[pTT].append({'frame'  : frame, 
                                                  'channel': channel, 
                                                  'n_link' : n_link,})

        S1_index += 1
    return [reversed_data_pTT]

def get_pTT_id(Sector, S1Board, CEECEH, x):
    eta = x[x.find('eta')+3]
    if x[x.find('eta')+4] != '-':
        eta += x[x.find('eta')+4]
    phi = x[x.find('phi')+3]
    if x[x.find('phi')+4] != '*':
        phi += x[x.find('eta')+4]
    eta = int(eta)
    phi = int(phi)
    S1Board = (int(S1Board[2])*16 + int(S1Board[3])) & 0x3F
    return hex(0x00000000 | ((Sector & 0x3) << 29) | ((1 & 0x3) << 26)  | ((6 & 0xF) << 22) | ((S1Board & 0x3F) << 16) | ((CEECEH & 0x1) << 10) | ((eta & 0x1F) << 5) | ((phi & 0x1F) << 0))
    
def get_moduleCEE(x,Sector):
    start_cursor = 0
    end_cursor = x[start_cursor:].find(',') + start_cursor
    layer = int(x[start_cursor:end_cursor])
    start_cursor = x[end_cursor+1].find(',') + end_cursor + 1 +1
    end_cursor = x[start_cursor:].find(',') + start_cursor
    u = int(x[start_cursor:end_cursor])
    start_cursor = x[end_cursor+1].find(',') + end_cursor + 1 +1
    v = int(x[start_cursor:])
    module_id = hex(0x00000000 | ((Sector & 0x3) << 29) | ((0 & 0x3) << 26)  | ((0 & 0xF) << 22) | ((layer & 0x3F) << 16) |  ((u & 0xF) << 12) | ((v & 0xF) << 8))
    return(module_id,layer,u,v)
                                                                                                                                                 

def get_moduleCEH(x,Sector):
    start_cursor = 0
    end_cursor = x[start_cursor:].find(',') + start_cursor
    layer = int(x[start_cursor:end_cursor])
    start_cursor = x[end_cursor+1].find(',') + end_cursor + 1 +1
    end_cursor = x[start_cursor:].find(',') + start_cursor
    u = int(x[start_cursor:end_cursor])
    start_cursor = x[end_cursor+1].find(',') + end_cursor + 1 +1
    end_cursor = x[start_cursor:].find(',') + start_cursor
    v = int(x[start_cursor:end_cursor])
    start_cursor = x[end_cursor+1].find(',') + end_cursor + 1 +1
    stc = int(x[start_cursor:])
    module_id = hex(0x00000000 | ((Sector & 0x3) << 29) | ((0 & 0x3) << 26)  | ((0 & 0xF) << 22) | ((layer & 0x3F) << 16) |  ((u & 0xF) << 12) | ((v & 0xF) << 8))
    return(module_id,layer,u,v,stc)


def read_pTT(x,S1Board,CEECEH,Sector):
    pTT = {'pTT' :get_pTT_id(Sector,S1Board,CEECEH,x), 'Modules':[]}
    cursor = x.find('\t')+2
    nb_module = x[cursor]
    for k in range(nb_module): 
        start_module = x[cursor:].find('(') +1 +cursor
        end_module = x[cursor:].find(')')  +cursor
        energy = x[end_module+3,x[end_module+3:].find(',')]
        if CEECEH==0:
            module_id,layer,u,v = get_moduleCEE(x[start_module: end_module],Sector)
            pTT['Module'].append([{'module_id' : module_id, 'module_layer' : layer,'module_u' : u,'module_v' : v,'module_energy' : int(energy)}])
        if CEECEH==1:
            module_id,layer,u,v,stc_idx = get_moduleCEH(x[start_module: end_module],Sector)
            pTT['Module'].append([{'module_id' : module_id,'module_layer' : layer,'module_u' : u,'module_v' : v,'stc_idx': stc ,'module_energy' : int(energy)}])
        cursor = x[end_module+3:].find(',')
    return(pTT)
        
def read_txt_pTTs(Edges,Sector):
    if Edges == 'yes':
        fCEE = open('config_files/CE_E_allBoards_Edges.txt', 'r')
        data_CEE = fCEE.readlines()
        fCEE.close()
        fCEH = open('config_files/CE_H_allBoards_Edges.txt', 'r')
        data_CEH = fCEH.readlines()
        fCEH.close()
    if Edges == 'no':
        fCEE = open('config_files/CE_E_allBoards_NoEdges.txt', 'r')
        data_CEE = fCEE.readlines()
        fCEE.close()
        fCEH = open('config_files/CE_H_allBoards_NoEdges.txt', 'r')
        data_CEH = fCEH.readlines()
        fCEH.close()
    pTTs_CEE = []
    for x in data_CEE:
        if x[0:5] == 'Board':
            S1Board = x[6:16]
        if x[0] == '/':
            pTTs_CEE.append(read_pTT(x,S1Board,0,Sector))
    pTTs_CEH = []
    for x in data_CEH:
        if x[0:5] == 'Board':
            S1Board = x[6:16]
        if x[0] == '/':
            pTTs_CEH.append(read_pTT(x,S1Board,1,Sector))
    return(pTTs_CEE,pTTs_CEH)

print(read_txt_pTTs('no',0))
    
