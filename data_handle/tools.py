import yaml
import cppyy
from cppyy.gbl import l1thgcfirmware, std

def compress_value(value, exponent_bits=4, mantissa_bits=3, truncation_bits=0):
    saturation_code = (1 << (exponent_bits + mantissa_bits)) - 1
    saturation_value = ((1 << (mantissa_bits + truncation_bits + 1)) - 1) << ((1 << exponent_bits) - 2)

    if value > saturation_value:
        return saturation_value, saturation_code

    bitlen = 0
    shifted_value = int(value) >> truncation_bits
    valcopy = shifted_value
    while valcopy != 0:
        valcopy >>= 1
        bitlen += 1

    if bitlen <= mantissa_bits:
        compressed_code = shifted_value
        compressed_value = shifted_value << truncation_bits
        return compressed_value, compressed_code

    # Build exponent and mantissa
    exponent = bitlen - mantissa_bits
    mantissa = (shifted_value >> (exponent - 1)) & ~(1 << mantissa_bits)
   
    # Assemble floating-point
    compressed_value = ((1 << mantissa_bits) | mantissa) << (exponent - 1)
    compressed_code = (mantissa << exponent_bits) | exponent #(exponent << mantissa_bits) | mantissa
    return compressed_value, compressed_code

def define_map():
    # reading the yaml file
    with open('config.yaml', "r") as afile:
      params = yaml.safe_load(afile)["s2emu_config"]

    enum = cppyy.gbl.l1thgcfirmware.Step
    map_custom = std.map[enum, "unsigned int"]()

    for d in params['stepLatency']:
      for latency in d.keys():
        map_custom[int(vars(enum)[latency])] = d[latency]

    params['stepLatency'] = map_custom
    return params

def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    if iteration == total: 
        print()

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
                u,v = -v,u-v
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
                    u,v = -v,u-v
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
