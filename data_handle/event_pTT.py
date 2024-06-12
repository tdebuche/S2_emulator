from collections import defaultdict
import numpy as np
import awkward as ak
import uproot
import math
import yaml
from data_handle.S1simulator import build_pTTsCEE
import cppyy
from cppyy.gbl import l1thgcfirmware, std

import data_handle.plot_tools as plot
from data_handle.tools import compress_value, printProgressBar, getuvsector


class EventData():
    def __init__(self, ds_si, ds_sci, gen):
        self.ds_si  = ds_si
        self.ds_sci  = ds_sci
        self.ds_ts = None
        self.ds_stc = None
        self.ds_pTTsCEE = None
        self.ds_pTTsdupCEE = None
        self.ds_pTTsCEH = None
        self.ds_pTTsdupCEH = None
        self.gen     = gen
        self.event   = gen.event
        self.eta_gen = gen.good_genpart_exeta[0]
        self.phi_gen = gen.good_genpart_exphi[0]
        self.pT_gen  = self._compute_pt(self.eta_gen,
                             gen.good_genpart_energy[0])

        self.data_packer = None
        self.pTT_packer = None          ######################add##################
        self.LSB = 1/10000 # 100 keV
        self.LSB_r_z = 0.7/4096
        self.LSB_phi = np.pi/1944
        self.offset_phi = -0.8

    def _compute_pt(self, eta, energy):
        return energy/np.cosh(eta)

    def ObjectType(self, object_type):
        return ((object_type & 0xF) << 22)
    
    def get_module_id(self,Sector, plane, u, v):
        return hex(0x00000000 |  ((Sector & 0x3) << 29) | ((0 & 0x3) << 26)  | ((0 & 0xF) << 22) |  ((plane & 0x3F) << 16) | ((u & 0xF) << 12) | ((v & 0xF) << 8))

    def get_pTT_id(self, Sector , S1Board, CEECEH, eta, phi):
        return hex(0x00000000 | ((Sector & 0x3) << 29) | ((1 & 0x3) << 26)  | ((6 & 0xF) << 22) | ((S1Board & 0x3F) << 16) | ((CEECEH & 0x1) << 10) | ((eta & 0x1F) << 5) | ((phi & 0x1F) << 0))
    
    def get_MB_id(self, plane, v, MB):
        return MB[int(plane)][int(v)]
                
            
    def get_pTT_allocation(self, xml_allocation, pTT):
        return xml_allocation[pTT]
        
    def get_pTT_duplication(self,xml_duplication,pTT):
        return xml_duplication[pTT]
        
    def get_TC_allocation(self, xml_data, module):
        return xml_data[module]

    def get_pTT_id(self, Sector, S1Board, CEECEH, eta,phi):
        S1Board = (int(S1Board[2])*16 + int(S1Board[3])) & 0x3F
        pTT_id = hex(0x00000000 | ((Sector & 0x3) << 29) | ((1 & 0x3) << 26)  | ((6 & 0xF) << 22) | ((S1Board & 0x3F) << 16) | ((CEECEH & 0x1) << 10) | ((eta & 0x1F) << 5) | ((phi & 0x1F) << 0))
        return pTT_id


    def provide_ts(self,args,xml):
        nb_selected_TCs = defaultdict(list)
        selected_TCs = self.ds_si
        for module_idx in range(len(self.ds_si.good_tc_layer)):
            u,v,sector = getuvsector(self.ds_si.good_tc_layer[module_idx][0],
                                        self.ds_si.good_tc_waferu[module_idx][0],
                                        self.ds_si.good_tc_waferv[module_idx][0])
            module_alloc = self.get_module_id(3,self.ds_si.good_tc_layer[module_idx][0],u,v)
            xml_alloc = self.get_TC_allocation(xml[0], module_alloc)
            module = self.get_module_id(sector,self.ds_si.good_tc_layer[module_idx][0],u,v)
            if xml_alloc: 
                n_TCs = xml_alloc[-1]['index']
                nb_selected_TCs[module].append[n_TCs]  

        TCs = self.ds_si
        ts = defaultdict(list)
        Sector = args.Sector
        for module_idx in range(len(self.ds_si.good_tc_layer)):
            u,v,sector = getuvsector(self.ds_si.good_tc_layer[module_idx][0],
                                        self.ds_si.good_tc_waferu[module_idx][0],
                                        self.ds_si.good_tc_waferv[module_idx][0])
            if u != -999:
                module = self.get_module_id(sector,self.ds_si.good_tc_layer[module_idx][0],u,v)
                if self.ds_si.good_tc_layer[module_idx][0] < 48:
                    if ts[module] == []:
                        ts[module].append(0)
                    #for idx in range(len(self.ds_si.good_tc_layer[module_idx])):
                    for idx in range(nb_selected_TCs[module][0],len(self.ds_si.good_tc_layer[module_idx])):
                        ts[module][0] += self.ds_si.good_tc_pt[module_idx][idx]
        self.ds_ts = ts
        return(nb_selected_TCs)
        
            
    def _process_eventpTT(self,args, xml_allocation,xml_duplication,S1pTTCEE,S1pTTCEH,S1pTTCEEdup,S1pTTCEHdup,xml):
        data_pTTs = defaultdict(list)
        Sector = args.Sector
        nb_selected_TCs = self.provide_ts(args,,xml)

        #CEE

        #Sector 0
        self.ds_pTTsCEE = build_pTTsCEE(self.ds_ts, args, S1pTTCEE)  #from module sums
        self.ds_pTTsCEE  = add_TCs(self.ds_pTTsCEE,self.ds_si,nb_selected_TCs,0,'CEE') #add selected TCs
        pTTs = self.ds_pTTsCEE
    

        #Sector 1
        self.ds_pTTsdupCEE = build_pTTsCEE(self.ds_ts, args, S1pTTCEEdup)
        self.ds_pTTsdupCEE  = add_TCs(self.ds_pTTsdupCEE,self.ds_si,nb_selected_TCs,1,'CEE')
        pTTsdup = self.ds_pTTsdupCEE

        #fill CEE links 
        
        for pTT_idx in range(len(pTTs)):
            pTT = pTTs[pTT_idx]['pTT_id']
            pTT_xml = self.get_pTT_allocation(xml_allocation, pTT)
            if pTT_xml != [] :    #if pTT is allocated in the 4 links
                data_pTTs[(pTT_xml[0]['frame'],pTT_xml[0]['n_link'],pTT_xml[0]['channel']%2)].append(pTTs[pTT_idx]['energy'])
        for pTT_idx in range(len(pTTsdup)):
            pTT = pTTsdup[pTT_idx]['pTT_id']
            pTT_xml = self.get_pTT_duplication(xml_duplication, pTT)
            if pTT_xml != [] :    #if pTT is allocated in the 2 links
                if data_pTTs[(pTT_xml[0]['frame'], pTT_xml[0]['n_link'],pTT_xml[0]['channel']%2)] == []:
                    data_pTTs[(pTT_xml[0]['frame'], pTT_xml[0]['n_link'],pTT_xml[0]['channel']%2)].append(pTTsdup[pTT_idx]['energy'])


        #CEH

        #Sector 0 
        self.ds_pTTsCEH = build_pTTsCEE(self.ds_ts, args, S1pTTCEH)
        self.ds_pTTsCEH  = add_TCs(self.ds_pTTsCEH,self.ds_si,nb_selected_TCs,0,'CEH')
        pTTs = self.ds_pTTsCEE

        #Sector 1
        self.ds_pTTsdupCEH = build_pTTsCEE(self.ds_ts, args, S1pTTCEHdup)
        self.ds_pTTsdupCEH  = add_TCs(self.ds_pTTsdupCEH,self.ds_si,nb_selected_TCs,1,'CEH')
        pTTsdup = self.ds_pTTsdupCEH

        #fill CEH links
        
        for pTT_idx in range(len(pTTs)):
            pTT = pTTs[pTT_idx]['pTT_id']
            pTT_xml = self.get_pTT_allocation(xml_allocation, pTT)
            if pTT_xml != [] :    #if pTT is allocated in the 4 links
                data_pTTs[(pTT_xml[0]['frame'],pTT_xml[0]['n_link'],pTT_xml[0]['channel']%2)].append(pTTs[pTT_idx]['energy'])
        for pTT_idx in range(len(pTTsdup)):
            pTT = pTTsdup[pTT_idx]['pTT_id']
            pTT_xml = self.get_pTT_duplication(xml_duplication, pTT)
            if pTT_xml != [] :    #if pTT is allocated in the 2 links
                if data_pTTs[(pTT_xml[0]['frame'], pTT_xml[0]['n_link'],pTT_xml[0]['channel']%2)] == []:
                    data_pTTs[(pTT_xml[0]['frame'], pTT_xml[0]['n_link'],pTT_xml[0]['channel']%2)].append(pTTsdup[pTT_idx]['energy'])

        return data_pTTs


    def _pTT_packer(self, args, xml_allocation,xml_duplication,S1pTTCEE,S1pTTCEH,S1pTTCEEdup,S1pTTCEHdup,xml):
        data_pTTs = self._process_eventpTT(args, xml_allocation,xml_duplication,S1pTTCEE,S1pTTCEH,S1pTTCEEdup,S1pTTCEHdup,xml)
        self.pTT_packer =  data_pTTs

################################TCs###################################################

    def _process_module(self, ds_TCs, idx, xml_alloc, data_TCs):
        n_TCs = xml_alloc[-1]['index']  # dangerous
        columns = [frame['column'] for frame in xml_alloc]
   
        # simulating the BC algorithm (ECON-T) and the phi sorting in the S1 FPGA
        mod_phi = ds_TCs.good_tc_phi[idx][:n_TCs+1]
        mod_energy = ds_TCs.good_tc_pt[idx][:n_TCs+1][ak.argsort(mod_phi)]
        mod_r_over_z = ds_TCs.r_over_z[idx][:n_TCs+1][ak.argsort(mod_phi)]
        mod_phi = ak.sort(mod_phi)

        # assigning each TCs to a columns
        xml_alloc = sorted(xml_alloc, key=lambda x: x['column'])
        
        for tc_idx, TC_xml in enumerate(xml_alloc):
            if tc_idx > len(mod_energy)-1: break
            n_link = TC_xml['n_link']
    
            value_energy, code_energy = compress_value(mod_energy[tc_idx]/self.LSB)
            value_r_z = int(mod_r_over_z[tc_idx]/self.LSB_r_z) & 0xFFF # 12 bits
            value_phi = int((mod_phi[tc_idx]-self.offset_phi)/self.LSB_phi) & 0xFFF # 12 bits

            data_TCs[(TC_xml['frame'],n_link,TC_xml['channel']%3)] = [ 
                code_energy, value_r_z, value_phi
                ]
 
    def _process_event(self, args, xml, MB_conv):
        data_TCs = defaultdict(list)
        for module_idx in range(len(self.ds_si.good_tc_layer)):
            u,v,sector = getuvsector(self.ds_si.good_tc_layer[module_idx][0],
                                            self.ds_si.good_tc_waferu[module_idx][0],
                                            self.ds_si.good_tc_waferv[module_idx][0])
            if sector ==0:
                module = self.get_module_id(3,self.ds_si.good_tc_layer[module_idx][0],u,v)
                xml_alloc = self.get_TC_allocation(xml[0], module)
                if xml_alloc: 
                    self._process_module(self.ds_si, module_idx, xml_alloc, data_TCs)
                    print(self.ds_si.good_tc_layer[module_idx][0],
                                                self.ds_si.good_tc_waferu[module_idx][0],
                                                self.ds_si.good_tc_waferv[module_idx][0])

        for MB_idx in range(len(self.ds_sci.good_tc_layer)):
            MB = self.get_MB_id(self.ds_sci.good_tc_layer[MB_idx][0],
                                self.ds_sci.MB_v[MB_idx][0], MB_conv)
            xml_alloc = self.get_TC_allocation(xml[1], MB)
            if xml_alloc: self._process_module(self.ds_sci, MB_idx, xml_alloc, data_TCs)
        return data_TCs


    
    def _data_packer(self, args, xml, xml_MB):
        data_TCs = self._process_event(args, xml, xml_MB)
        self.data_packer = data_TCs
         

#######################################################################################
############################### PROVIDE EVENTS ########################################
#######################################################################################
with open('config.yaml', "r") as afile:
    cfg_particles = yaml.safe_load(afile)["particles"]

def apply_sort(df, counts, axis):
    for field in df.fields:
        df[field] = ak.unflatten(df[field], counts, axis)
    return df

    
            

def provide_event(ev, gen):
    ev['r_over_z'] = np.sqrt(ev.good_tc_x**2 + ev.good_tc_y**2)/ev.good_tc_z
    ev['MB_v'] = np.floor((ev.good_tc_cellv-1)/4)
    #ev = ev[[x for x in ak.fields(ev) if not x in ["good_tc_x","good_tc_y","good_tc_z"]]]
    
    # dividing silicon and scintillators
    sci = ev[ev['good_tc_subdet'] == 10]
    si  = ev[ev['good_tc_subdet'] != 10]
    
    # selecting first 120 sector only
    sci = sci[sci['good_tc_cellv']<=48]
    #si = si[si['good_tc_layer'] < 27]

    # sorting by modules  
    sorted_waferu = si[ak.argsort(si['good_tc_waferu'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferu.good_tc_waferu), axis=None)
    sorted_si = apply_sort(sorted_waferu, counts, 1)

    sorted_waferv = sorted_si[ak.argsort(sorted_si['good_tc_waferv'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferv.good_tc_waferv), axis=None)
    sorted_si = apply_sort(sorted_waferv, counts, 2)

    sorted_layer = sorted_si[ak.argsort(sorted_si['good_tc_layer'])]
    counts = ak.flatten(ak.run_lengths(sorted_layer.good_tc_layer), axis=None)
    sorted_si = apply_sort(sorted_layer, counts, 3)
    sorted_si = ak.flatten(sorted_si, axis=3)
    sorted_si = ak.flatten(sorted_si, axis=2)

    # sorting by transverse energy, simulating the ECONT_T
    sorted_si = sorted_si[ak.argsort(sorted_si['good_tc_pt'], ascending=False)][0]
    
    # sorting sci by MB (cellv) and plane
    sorted_MB = sci[ak.argsort(sci['MB_v'])]
    counts = ak.flatten(ak.run_lengths(sorted_MB.MB_v), axis=None)
    sorted_sci = apply_sort(sorted_MB, counts, 1)

    sorted_layer = sorted_sci[ak.argsort(sorted_sci['good_tc_layer'])]
    counts = ak.flatten(ak.run_lengths(sorted_layer.good_tc_layer), axis=None)
    sorted_sci = apply_sort(sorted_layer, counts, 2)
    sorted_sci = ak.flatten(sorted_sci, axis=2)

    # sorting by transverse energy, simulating the ECONT_T
    sorted_sci = sorted_sci[ak.argsort(sorted_sci['good_tc_pt'], ascending=False)][0]
  
    return EventData(sorted_si, sorted_sci,gen)



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

    branches_gen = [
        'event', 'good_genpart_exeta', 'good_genpart_exphi', 'good_genpart_energy'
    ]

    tree = uproot.open(filepath)[name_tree]
    events_ds = []
    printProgressBar(0, n, prefix='Reading '+str(n)+' events from ROOT file:', suffix='Complete', length=50)
    for ev in range(n):
      data = tree.arrays(branches_tc, entry_start=ev, entry_stop=ev+1, library='ak')
      data_gen = tree.arrays(branches_gen, entry_start=ev, entry_stop=ev+1, library='ak')[0]
      events_ds.append(provide_event(data, data_gen))
      printProgressBar(ev+1, n, prefix='Reading '+str(n)+' events from ROOT file:', suffix='Complete', length=50)
    return events_ds

