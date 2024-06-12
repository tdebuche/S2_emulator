import numpy as np
import awkward as ak
from data_handle.tools import getuvsector
from collections import defaultdict



def build_pTTsCEE(ts_energy, args, S1pTTCEE):
    if args.Edges == 'yes': nb_phi = 28
    else : nb_phi = 24
    Sector = args.Sector
    pTTsCEE = []
    for pTT_idx in range(len(S1pTTCEE)):
        energyCEE = 0
        pTT_id = S1pTTCEE[pTT_idx]['pTT'] 
        ModulesCEE = S1pTTCEE[pTT_idx]['Modules'] 
        for module_idx in range(len(ModulesCEE)):
            module_id = ModulesCEE[module_idx]['module_id']
            energy = ModulesCEE[module_idx]['module_energy']
            if ts_energy[module_id] != []:
                energyCEE += ts_energy[module_id][0] * energy/16
        pTTsCEE.append({'pTT_id' : pTT_id, 'energy': energyCEE})
            #if energyCEE != 0:
            #print({'pTT_id' : pTT_id, 'energy': energyCEE})
    return(pTTsCEE)


def build_pTTsCEH(stc_energy,args,S1pTTCEH):
    if args.Edges == 'yes': nb_phi = 28
    else : nb_phi = 24
    Sector = args.Sector
    pTTsCEH = []
    for pTT_idx in range(len(S1pTTCEH)):
        energyCEH = 0
        pTT_id = S1pTTCEH[pTT_idx]['pTT'] 
        STCsCEH = S1pTTCEH[pTT_idx]['Modules'] 
        for stc_idx in range(len(STCsCEH)):
            module_id = ModulesCEE[stc_idx]['stc_id']
            stc = ModulesCEE[stc_idx]['module_idx']
            energy = ModulesCEH[module_idx]['stc_energy']
            if stc_energy[module_id] != []:
                energyCEH += stc_energy[module_id][0][stc] * energy/16
        pTTsCEE.append({'pTT_id' : pTT_id, 'energy': energyCEH})
    return(pTTsCEH)

def add_TCs(pTTs,TCs,nb_selected_TCs, Sector,CEECEH):
    energytoadd= defaultdict(list)
    for module_idx in range(len(TCs.good_tc_layer)):
        if ((TCs.good_tc_layer[module_idx][0] < 27) and (CEECEH == 'CEE' )) or ((TCs.good_tc_layer[module_idx][0] >= 27) and (CEECEH == 'CEH' )):
            u,v,sector = getuvsector(TCs.good_tc_layer[module_idx][0],
                                        TCs.good_tc_waferu[module_idx][0],
                                        TCs.good_tc_waferv[module_idx][0])
            module = get_module_id(Sector, TCs.good_tc_layer[module_idx][0], u, v)
            if sector == Sector:
                for idx in range(len(nb_selected_TCs[module])):
                    eta,phi = getetaphi(TCs.good_tc_phi[module_idx][idx],TCs.r_over_z[module_idx][idx])
                    S1_Board = S1_Board(TCs.good_tc_layer[module_idx][idx])
                    if CEECEH == 'CEE': a = 0
                    if CEECEH == 'CEH': a = 1
                    pTT = get_pTT_id(Sector, S1Board, a, eta,phi)
                    energytoadd[pTT].append(TCs.good_tc_pt[module_idx][idx])
    return pTTs


def get_pTT_id(Sector, S1Board, CEECEH, eta,phi):
    S1Board = (int(S1Board[4],16)*16 + int(S1Board[5],16)) & 0x3F
    pTT_id = hex(0x00000000 | ((Sector & 0x3) << 29) | ((1 & 0x3) << 26)  | ((6 & 0xF) << 22) | ((S1Board & 0x3F) << 16) | ((CEECEH & 0x1) << 10) | ((eta & 0x1F) << 5) | ((phi & 0x1F) << 0))
    #while len(pTT_id) <10:
    #    pTT_id = '0x'+ str(0) +pTT_id[2:]
    return pTT_id

def get_module_id(Sector, plane, u, v):
    return hex(0x00000000 |  ((Sector & 0x3) << 29) | ((0 & 0x3) << 26)  | ((0 & 0xF) << 22) |  ((plane & 0x3F) << 16) | ((u & 0xF) << 12) | ((v & 0xF) << 8))

def getetaphi(phi,roverz):
    teta = np.arccos(1/roverz)
    eta = -np.log(np.tan(teta/2))
    eta = int((eta-1.305)/(np.pi/36)) #1.305 offset
    phi = int((phi+ (15*np.pi/180))/(np.pi/36) ) # -15° offset
    return(eta,phi)

