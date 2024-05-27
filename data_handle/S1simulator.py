import numpy as np
import awkward as ak

def build_pTTsCEE(ts_erngy, args, S1pTTCEE):
    if args.Edges == 'yes': nb_phi = 28
    else : nb_phi = 24
    Sector = args.Sector
    pTTsCEE = []
    for S1Board in range(14):
        for pTT_idx in range(len(S1pTTCEE):
            energyCEE = 0
            pTT_id = S1pTTCEE[pTT_idx][0] 
            ModulesCEE = S1pTTCEE[pTT_idx][1] 
            for module_idx in range(len(ModulesCEE)):
                module_id = ModulesCEE[module_idx]['module_id']
                energy = ModulesCEE[module_idx]['module_energy']
                energyCEE += ts_energy[module_id] * energy/16
            pTTsCEE.append({'pTT' : pTT_id, 'energy': energy})
    return(pTTsCEE)


def build_pTTsCEH(stc_energy,args,S1pTTCEH):
    if args.Edges == 'yes': nb_phi = 28
    else : nb_phi = 24
    Sector = args.Sector
    pTTsCEH = []
    for pTT_idx in range(len(S1pTTCEH):
        energyCEH = 0
        pTT_id = S1pTTCEH[pTT_idx][0] 
        STCsCEH = S1pTTCEH[pTT_idx][1] 
        for stc_idx in range(len(STCsCEH)):
            module_id = ModulesCEE[stc_idx]['stc_id']
            stc = ModulesCEE[stc_idx]['module_idx']
            energy = ModulesCEH[module_idx]['stc_energy']
            energyCEH += stc_energy[module_id][stc] * energy/16
        pTTsCEE.append({'pTT' : pTT_id, 'energy': energy})
    return(pTTsCEH)
                    
