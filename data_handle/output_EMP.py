from collections import defaultdict
import numpy as np
import awkward as ak
import uproot
import math
import yaml
import os

def createEMPfile(event):
    pTTlinks = event.pTT_packer 
    
    with open('./EMPfiles/EMPtest.txt', 'a') as file:   
      #write the first comments for the EMP file
        num_columns=56
        file.write(f"ID: x1 \n")
        file.write(f"Metadata: (strobe,) start of orbit, start of packet, end of packet, valid \n \n")
        column_str = '          '.join(str(j).zfill(3).rjust(20) for j in range(num_columns))
        file.write(f"      Link {column_str}\n")
        for frame_idx in range(108):
            if frame_idx == 0:
                metadata = 1101 
            elif frame_idx != 107:
                metadata = 1
            elif frame_idx == 107:
                metadata = 11
            
            frame = 'Frame '+ str(frame_idx).zfill(4).rjust(4)
            #frame += ' '.join(str(metadata).zfill(4) + " " + str(f'{word(pTTlinks,frame_idx,j):016x}' ) + " " for j in range(84))
            frame += ' '.join(str(metadata).zfill(4) + " " + str(f'{0:016x}' ) + " " for j in range(84))
            file.write(f" {frame} \n")
        file.close()

def word(pTTlinks,frame_idx,nb_link):
    energypTT0 = pTTlinks((frame_idx,nb_link,0))[0]
    energypTT1 = pTTlinks((frame_idx,nb_link,1))[0]
    return(hex(0x0000000000000000|(energypTT0) << 53 | (energypTT0)<< 45))
    

    # with open('output.txt', 'a') as file:  
    #   #metadata = [1101 if i == 0 else 11 if i == 107 else 1 for i in range(num_rows)],below only for the last line of the EMP file
    #   metadata = 11 
    # # Iterate over each row and write the data for each link
    #   LinksInData_2d = [0] * 56       
    #   row_str = ' '.join(str(metadata).zfill(4) + " " + str(f'{LinksInData_2d[j]:016x}' ) + " " for j in range(56))
    #   # Write the last row to the file
    #   frame = 108*nevents
    #   file.write(f"Frame " + str(frame).zfill(4).rjust(4) + "    " + row_str + '\n')            
    #   file.close()
