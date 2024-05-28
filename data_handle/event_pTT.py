from collections import defaultdict
import numpy as np
import awkward as ak
import uproot
import math
import yaml
import data_handle.S1simulator

import data_handle.plot_tools as plot
from data_handle.tools import compress_value, printProgressBar
        

class EventData():
    def __init__(self, ds_si, ds_sci, gen):
        self.ds_si  = ds_si
        self.ds_sci  = ds_sci
        self.ds_ts = None
        self.ds_stc = None
        self.ds_pTTs = None
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
        # CMSSW to our u v convention u'=v-u, v'=v
        # print('Analysing module ', plane, v-u, v)
        if plane & ~0x3F : return 0 # raise Exception( "Invalid plane" )
        if v-u & ~0xF : return 0 # raise Exception( "Invalid u" )
        if v   & ~0xF : return 0 # raise Exception( "Invalid v" )
        return hex(0x00000000 |  ((Sector & 0x3) << 29) | ((0 & 0x3) << 26)  | ((0 & 0xF) << 22) | ((plane & 0x3F) << 16) | ((plane & 0x3F) << 16) | ((v-u & 0xF) << 12) | ((v & 0xF) << 8))

    def get_pTT_id(self, Sector , S1Board, CEECEH, eta, phi):
        return hex(0x00000000 | ((Sector & 0x3) << 29) | ((1 & 0x3) << 26)  | ((6 & 0xF) << 22) | ((S1Board & 0x3F) << 16) | ((CEECEH & 0x1) << 10) | ((eta & 0x1F) << 5) | ((phi & 0x1F) << 0))
    
    def get_MB_id(self, plane, v, MB):
        return MB[int(plane)][int(v)]

    def provide_ts(self,args):
        TCs = self.ds_si
        ts = defaultdict(list)
        Sector = args.Sector
        for module_idx in range(len(self.ds_si.good_tc_layer)):
            module = self.get_module_id(Sector,
                                        self.ds_si.good_tc_layer[module_idx][0],
                                        self.ds_si.good_tc_waferu[module_idx][0],
                                        self.ds_si.good_tc_waferv[module_idx][0])
            if ts[module] == []:
                ts[module].append(0)
            ts[module][0] += self.ds_si.good_tc_pt[module_idx][0]
        self.ds_ts = ts
            
                
            
    def get_pTT_allocation(self, xml_allocation, pTT):
        return xml_allocation[pTT]

    def get_pTT_id(self, Sector, S1Board, CEECEH, eta,phi):
        S1Board = (int(S1Board[2])*16 + int(S1Board[3])) & 0x3F
        pTT_id = hex(0x00000000 | ((Sector & 0x3) << 29) | ((1 & 0x3) << 26)  | ((6 & 0xF) << 22) | ((S1Board & 0x3F) << 16) | ((CEECEH & 0x1) << 10) | ((eta & 0x1F) << 5) | ((phi & 0x1F) << 0))
        #while len(pTT_id) <10:
        #    pTT_id = '0x'+ str(0) +pTT_id[2:]
        return pTT_id
        
    def _process_eventpTT(self,args, xml_allocation, xml_duplication,S1pTTCEE,S1pTTCEH):
        data_pTTs = std.map[int,std.map[int,std.map[int,'std::vector<long int>']]]()
        Sector = args.Sector
        self.ds_pTTs = S1simulator.build_pTTsCEE(self.ds_ts, args, S1pTTCEE)
        pTTs = self.ds_pTTs
        
        for pTT_idx in range(len(pTTs)):
            pTT = pTTs[pTT_index]['pTT']
            pTT_xml = self.get_pTT_allocation(xml_allocation, pTT)
            data_pTTs[pTT_xml['frame']][pTT_xml['n_link']][pTT_xml['channel']%2] = pTT = pTTs[pTT_index]['energy']

        return data_pTTs


    def _pTT_packer(self, args, xml_allocation,S1pTTCEE,S1pTTCEH):
        self.provide_ts(args)
        data_pTTs = self._process_eventpTT(args,args, xml_allocation,S1pTTCEE,S1pTTCEH)
        self.data_packerpTT =  data_pTTs


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
    ev = ev[[x for x in ak.fields(ev) if not x in ["good_tc_x","good_tc_y","good_tc_z"]]]
    
    # dividing silicon and scintillators
    sci = ev[ev['good_tc_subdet'] == 10]
    si  = ev[ev['good_tc_subdet'] != 10]
    
    # selecting first 120 sector only
    si  = si[(si['good_tc_waferv']-si['good_tc_waferu']>0) & (si['good_tc_waferv']>=0)]
    sci = sci[sci['good_tc_cellv']<=48]

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

