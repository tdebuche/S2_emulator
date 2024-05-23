import numpy as np
import awkward as ak
import uproot
import math
import yaml
import cppyy
from cppyy.gbl import l1thgcfirmware, std

import data_handle.plot_tools as plot
from data_handle.tools import compress_value, printProgressBar
        
cppyy.cppdef("""
#import "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h"
std::vector<std::unique_ptr<l1thgcfirmware::HGCalLinkTriggerCell>> fill_Links(std::map<int, std::vector<long int>>& data_links) {
std::vector<std::unique_ptr<l1thgcfirmware::HGCalLinkTriggerCell>> LinksInData;
for (int link_frame = 0; link_frame < 84 * 108; ++link_frame) {
    LinksInData.push_back(std::make_unique<l1thgcfirmware::HGCalLinkTriggerCell>());
    if (data_links.find(link_frame) != data_links.end()) {
        LinksInData.back()->data_     = data_links[link_frame][0];
        LinksInData.back()->r_over_z_ = data_links[link_frame][1];
        LinksInData.back()->phi_      = data_links[link_frame][2];
    }
};
return LinksInData;
}""")

cppyy.cppdef("""
std::vector<long int> create_link(const std::map<int, std::vector<long int>>& data_in) {
    int64_t energy = 0, r_over_z = 0, phi = 0;
    for (const auto& pair : data_in) {
        int channel = pair.first;
        energy   = energy | (pair.second[0] << (channel*15));
        r_over_z = r_over_z | (pair.second[1] << (channel*15));
        phi      = phi | (pair.second[2] << (channel*15));
    }
    return {energy, r_over_z, phi};
}

std::map<int, std::vector<long int>> pack_Links(std::map<int, std::map<int, std::map<int, std::vector<long int>>>>& data_TCs) {
std::map<int, std::vector<long int>> data_links;
// Iterate over data_TCs
for (const auto& frame_pair : data_TCs) {
    int frame = frame_pair.first;
    for (const auto& n_link_pair : frame_pair.second) {
        int n_link = n_link_pair.first;
        // Create link data
        std::vector<long int> link_data = create_link(n_link_pair.second);
        // Insert into data_links
        data_links[84 * frame + n_link] = link_data;
    }
}
return data_links;
}""")

class EventData():
    def __init__(self, ds_si, ds_sci, gen):
        self.ds_si  = ds_si
        self.ds_sci  = ds_sci
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
    
    def get_module_id(self, plane, u, v):
        # CMSSW to our u v convention u'=v-u, v'=v
        # print('Analysing module ', plane, v-u, v)
        if plane & ~0x3F : return 0 # raise Exception( "Invalid plane" )
        if v-u & ~0xF : return 0 # raise Exception( "Invalid u" )
        if v   & ~0xF : return 0 # raise Exception( "Invalid v" )
        return hex(0x60000000 | self.ObjectType(0) | ((plane & 0x3F) << 16) | ((v-u & 0xF) << 12) | ((v & 0xF) << 8))
    
    def get_MB_id(self, plane, v, MB):
        return MB[int(plane)][int(v)]
    
    def get_TC_allocation(self, xml_data, module):
        return xml_data[module]
    
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

            data_TCs[TC_xml['frame']][n_link][TC_xml['channel']%3] = [ 
                code_energy, value_r_z, value_phi
                ]

    def _process_pTT(self, ds_TCs, idx, xml_alloc, data_TCs):
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

            data_TCs[TC_xml['frame']][n_link][TC_xml['channel']%3] = [ 
                code_energy, value_r_z, value_phi
                ]


 
    def _process_event(self, args, xml, MB_conv):
        data_TCs = std.map[int,std.map[int,std.map[int,'std::vector<long int>']]]()

        for module_idx in range(len(self.ds_si.good_tc_layer)):
            module = self.get_module_id(self.ds_si.good_tc_layer[module_idx][0],
                                        self.ds_si.good_tc_waferu[module_idx][0],
                                        self.ds_si.good_tc_waferv[module_idx][0])
            xml_alloc = self.get_TC_allocation(xml[0], module)
            if xml_alloc: self._process_module(self.ds_si, module_idx, xml_alloc, data_TCs)

        for MB_idx in range(len(self.ds_sci.good_tc_layer)):
            MB = self.get_MB_id(self.ds_sci.good_tc_layer[MB_idx][0],
                                self.ds_sci.MB_v[MB_idx][0], MB_conv)
            xml_alloc = self.get_TC_allocation(xml[1], MB)
            if xml_alloc: self._process_module(self.ds_sci, MB_idx, xml_alloc, data_TCs)
        return data_TCs

    def _process_eventpTT(self,args, xml_allocation, xml_duplication,txt_S1pTTs)
        data_TCs = std.map[int,std.map[int,std.map[int,'std::vector<long int>']]]()
        Sector = args.sector
        pTTs = self.build_pTTs(args,txt_S1pTTs)

        for pTT_idx in range(len(pTTs)):
            pTT = self.get_pTT_id(pTTs['pTT_Board'][pTT_idx],pTTs['pTT_CEECEH'][pTT_idx],pTTs['pTT_eta'][pTT_idx],pTTs['pTT_phi'][pTT_idx])
            xml_alloc = self.get_pTT_allocation(xml_allocation, pTT)
            if xml_alloc: self._process_pTT(self.ds_si, module_idx, xml_alloc, data_TCs)

        for MB_idx in range(len(self.ds_sci.good_tc_layer)):
            MB = self.get_MB_id(self.ds_sci.good_tc_layer[MB_idx][0],
                                self.ds_sci.MB_v[MB_idx][0], MB_conv)
            xml_alloc = self.get_TC_allocation(xml[1], MB)
            if xml_alloc: self._process_module(self.ds_sci, MB_idx, xml_alloc, data_TCs)
        return data_TCs

    
    def _data_packer(self, args, xml, xml_MB):
        data_TCs = self._process_event(args, xml, xml_MB)
        data_links = cppyy.gbl.pack_Links(data_TCs)

        # filling data into emulator c++ variables
        self.data_packer = cppyy.gbl.fill_Links(data_links) 

        
    def _pTT_packer(self, args, xml_allocation, xml_duplication,txt_S1pTTs):
        data_pTTs = self._process_eventpTT(args,args, xml_allocation, xml_duplication,txt_S1pTTs)
        data_links = cppyy.gbl.pack_LinkspTT(data_pTTs)

        # filling data into emulator c++ variables
        self.data_packerpTT = cppyy.gbl.fill_LinkspTT(data_links) 
        




#######################################################################################
############################### PROVIDE EVENTS PTT ########################################
#######################################################################################

with open('config.yaml', "r") as afile:
    cfg_particles = yaml.safe_load(afile)["particles"]

def apply_sort(df, counts, axis):
    for field in df.fields:
        df[field] = ak.unflatten(df[field], counts, axis)
    return df

def provide_eventpTT(ev,ev_pTT, gen):
    ev['r_over_z'] = np.sqrt(ev.good_tc_x**2 + ev.good_tc_y**2)/ev.good_tc_z
    ev['MB_v'] = np.floor((ev.good_tc_cellv-1)/4)
    ev = ev[[x for x in ak.fields(ev) if not x in ["good_tc_x","good_tc_y","good_tc_z"]]]
    
    # dividing CE-E and CE-H
    ts = ev_pTT[ev_pTT['ts_layer']<27]
    stc = ev[ev['ts_layer']>26]

    # sorting ts
    sorted_waferu = ts[ak.argsort(ts['ts_waferu'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferu.ts_waferu), axis=None)
    sorted_ts = apply_sort(sorted_waferu, counts, 1)
        
    sorted_waferv = sorted_ts[ak.argsort(sorted_ts['ts_waferv'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferv.ts_waferv), axis=None)
    sorted_ts = apply_sort(sorted_waferv, counts, 2)   

    sorted_layer = sorted_ts[ak.argsort(sorted_ts['ts_layer'])]
    counts = ak.flatten(ak.run_lengths(sorted_layer.ts_layer), axis=None)
    sorted_ts = apply_sort(sorted_layer, counts, 3)
    sorted_ts = ak.flatten(sorted_ts, axis=3)
    sorted_ts = ak.flatten(sorted_ts, axis=2)


    # sorting stc
    sorted_waferu = stc[ak.argsort(stc['good_tc_waferu'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferu.good_tc_waferu), axis=None)
    sorted_stc = apply_sort(sorted_waferu, counts, 1)

    sorted_waferv = sorted_stc[ak.argsort(sorted_stc['good_tc_waferv'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferv.good_tc_waferv), axis=None)
    sorted_stc = apply_sort(sorted_waferv, counts, 2)   

    sorted_layer = sorted_stc[ak.argsort(sorted_stc['good_tc_layer'])]
    counts = ak.flatten(ak.run_lengths(sorted_layer.good_tc_layer), axis=None)
    sorted_stc = apply_sort(sorted_layer, counts, 3)
    sorted_stc = ak.flatten(sorted_stc, axis=3)
    sorted_stc = ak.flatten(sorted_stc, axis=2)

    return EventData(sorted_ts,sorted_stc, gen)

def provide_eventspTT(n, particles, PU):
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
    branches_pTT =  [
        'ts_energy', 'ts_waferu','ts_waferv','ts_layer',
    ]

    tree = uproot.open(filepath)[name_tree]
    events_ds = []
    printProgressBar(0, n, prefix='Reading '+str(n)+' events from ROOT file:', suffix='Complete', length=50)
    for ev in range(n):
      data = tree.arrays(branches_tc, entry_start=ev, entry_stop=ev+1, library='ak')
      data_pTT = tree.arrays(branches_pTT, entry_start=ev, entry_stop=ev+1, library='ak')
      data_gen = tree.arrays(branches_gen, entry_start=ev, entry_stop=ev+1, library='ak')[0]
      events_ds.append(provide_eventpTT(data, data_pTT, data_gen))
      printProgressBar(ev+1, n, prefix='Reading '+str(n)+' events from ROOT file:', suffix='Complete', length=50)
    return events_ds
        
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

    return EventData(sorted_si, sorted_sci, gen)

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
