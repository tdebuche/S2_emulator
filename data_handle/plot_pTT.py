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
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap



def read_xml_plot(Edges,Sector):
    if Edges == 'yes':
        tree = ET.parse('config_files/pTTs/Sector'+str(Sector)+'/AllocationpTTsEdges.xml')
    if Edges == 'no':
        tree = ET.parse('config_files/pTTs/Sector'+str(Sector)+'/AllocationpTTsNoEdges.xml')
        
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
                    Sector,S1Board,eta,phi,CEECEH = get_pTT_numbers(pTT)
                    data_pTT[(Sector,S1Board,eta,phi,CEECEH )].append((frame,n_link,channel%2))

    
    if Edges == 'yes':
        tree = ET.parse('config_files/pTTs/Sector'+str(Sector)+'/DuplicationpTTsEdges.xml')
    if Edges == 'no':
        tree = ET.parse('config_files/pTTs/Sector'+str(Sector)+'/DuplicationpTTsNoEdges.xml')
        
    root = tree.getroot()

    S1_index = 0
    for s1_element in root.findall('.//S1'):
        for channel_element in s1_element.findall('.//Channel'):
            channel = int(channel_element.get('aux-id'))
            for frame_element in channel_element.findall('.//Frame'):
                if all(attr in frame_element.attrib for attr in ['id','pTT']):
                    frame  = int(frame_element.get('id'))
                    pTT     = frame_element.get('pTT')
                    n_link = 14 + 14*math.floor(channel/2) + S1_index
                    Sector,S1Board,eta,phi,CEECEH = get_pTT_numbers(pTT)
                    data_pTT[(Sector,S1Board,eta,phi,CEECEH )].append((frame,n_link,channel%2))
                    

        S1_index += 1
    return data_pTT


def create_energies(data_links,etaphi_links,args):
    Edges = args.Edges
    if Edges == 'yes': 
        nb_phi = 28
        offset = 3
    else : 
        nb_phi = 24
        offset = 0
        
    Sector = args.Sector
    energies = [[0 for phi in range(36)]for eta in range(20)]
    for S1Board in range(14):
        for eta in range(20):
            for phi in range(36):
                if etaphi_links[(Sector,S1Board,eta,phi+offset,0)] != []:
                    if data_links[etaphi_links[(Sector,S1Board,eta,phi+offset,0)][0]] != []:
                        #energies[eta][phi] += data_links[etaphi_links[(Sector,S1Board,eta,phi+offset,0)][0]][0]
                        energies[eta][phi] += 1
    for S1Board in range(14):
        for eta in range(20):
            for phi in range(36):
                if etaphi_links[(Sector+1,S1Board,eta,phi-24+offset,0)] != []:
                    if data_links[etaphi_links[(Sector+1,S1Board,eta,phi-nb_phi+offset,0)][0]] != []:
                        #energies[eta][phi] += data_links[etaphi_links[(Sector+1,S1Board,eta,phi-nb_phi+offset,0)][0]][0]
                        energies[eta][phi] += 1
    return energies


def get_pTT_numbers(pTT):
    S1Board = int(pTT[4:6],16) & 0x3F
    phi = int(pTT,16) & 0x1F
    eta = (int(pTT,16) & 0x3E0) //(16 * 2)
    CEECEH = (int(pTT,16) & 0x400) //(16*16*4)
    Sector = (int(pTT[2],16) &  0x6)//2
    return(Sector,S1Board,eta,phi,CEECEH)

def create_bins(args):
    Edges = args.Edges
    if Edges == 'yes': 
        nb_phi = 28 + 9
        phimin =-15 * np.pi/180
    else : 
        nb_phi = 24 + 6
        phimin =  0 * np.pi/180
    nb_phi = 36
    phimin =  0 * np.pi/180
    etamin = 1.305
    L = [[[]for phi in range(nb_phi)]for eta in range(20)]
    BinsXY =[[[]for phi in range(nb_phi)]for eta in range(20)] 
    for eta in range(20):
        for phi in range(nb_phi):
            vertices = np.array([[eta * np.pi/36 + etamin,(eta+1) * np.pi/36 + etamin,
                    (eta+1) * np.pi/36 + etamin,eta * np.pi/36 + etamin],
                   [phi * np.pi/36  + phimin,phi * np.pi/36 + phimin,
                    (phi+1) * np.pi/36 + phimin,(phi+1) * np.pi/36 + phimin]])
            verticesXY = etaphitoXY(vertices[0],vertices[1],1)
            BinsXY[eta][phi].append(verticesXY[0].tolist()+[verticesXY[0][0]])
            BinsXY[eta][phi].append(verticesXY[1].tolist()+[verticesXY[1][0]])
    return BinsXY
            

    
def etaphitoXY(eta,phi,z):
    x = z * np.tan(2*np.arctan(np.exp(-eta))) * np.cos(phi)
    y = z * np.tan(2*np.arctan(np.exp(-eta))) * np.sin(phi)
    return(x,y)


def record_plot(data_links,etaphi_links,args,title):
    energies = create_energies(data_links,etaphi_links,args)
    BinXY = create_bins(args)
    plt.figure(figsize = (14,8))
    X =[]
    Y = []
    pointXY = [[],[]]
    weights = []
    weightmax = 0
    for eta in range(len(BinXY)):
        for phi in range(len(BinXY[0])):
            X.append(BinXY[eta][phi][0][0])
            Y.append(BinXY[eta][phi][1][0])
            if energies[eta][phi] != 100000:
                weights.append(energies[eta][phi])
                if  energies[eta][phi] > weightmax:
                    weightmax = energies[eta][phi]
            else : 
                weights.append(0)
            pointXY[0].append(np.sum(np.array(BinXY[eta][phi][0][0:4]))/4)
            pointXY[1].append(np.sum(np.array(BinXY[eta][phi][1][0:4]))/4)
    sc = plt.scatter(pointXY[0],pointXY[1],c=weights, vmin=0)
    plt.colorbar(sc)
    res = 0 
    colors = cm.get_cmap("viridis", 8)
    for eta in range(len(BinXY)):
        for phi in range(len(BinXY[0])):
            plt.plot(BinXY[eta][phi][0],BinXY[eta][phi][1],color = 'black')
            plt.fill(BinXY[eta][phi][0],BinXY[eta][phi][1],c = colors(weights[res]/weightmax))
            res +=1
            #if energies[eta][phi] != 100000:
                #plt.annotate(str(round(energies[eta][phi],2)),(np.sum(np.array(BinXY[eta][phi][0][0:4]))/4,np.sum(np.array(BinXY[eta][phi][1][0:4]))/4))
    plt.savefig('plot_pTTs/'+title +'.png')
    
