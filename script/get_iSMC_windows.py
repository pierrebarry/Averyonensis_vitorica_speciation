import gzip
import numpy as np
import os
from tqdm import tqdm
import pandas as pd
from itertools import repeat
import sys
from os import listdir
from os.path import isfile, join

samples = os.listdir("/shared/projects/origamis/wgs/output/iSMC/")
samples = "sphegodes_hiseq"
chrom = os.listdir("/shared/projects/origamis/wgs/output/iSMC/"+samples)
chrom = [x for x in chrom if x[-1]=="1"]

for ss in tqdm([samples]):
  for cc in tqdm(chrom,leave=False):
    rho_files = os.popen("ls -tr /shared/projects/origamis/wgs/output/iSMC/"+str(ss)+"/"+str(cc)+"/my_dataset_diploid_1_block_*_estimated_rho.txt.gz").read().split("\n")[:-1]
    #rho_files = [x for x in rho_files if "block_2" in x]
    
    block_size = 2000000
    
    for window in [10000,50000,100000,500000,1000000]:
      
      if os.path.isfile("/shared/projects/origamis/wgs/output/iSMC/"+str(cc)+"_"+str(ss)+"_"+str(window)+"kb.csv")==False:
        pos = 1
        POS = []
        RHO = []
        for j in tqdm(rho_files,leave=False):
          #rho_files_block = [rho_files[x] for x in range(0,len(block)) if block[x] == j]
          file = gzip.open(j, "r")
          content = np.array(list(map(float,file.read().decode().split("\n")[:-1])))
          file.close()
          if len(content)!=block_size:
            if (len(content))>=window:
              content = content[0:int(len(content)/window)*window]
              avgResult = list(np.average(content.reshape(-1, window), axis=1)) 
          else:
            avgResult = list(np.average(content.reshape(-1, window), axis=1)) 
          pp = list(range(pos,len(avgResult)*window+pos-1,window))
          pos = pp[-1]+window
          POS += pp
          RHO += avgResult
                #BLOCK += list(repeat(str(j), len(pp)))
        df2 = pd.DataFrame(
        {
          "POS": POS,
          "RHO": RHO
        }
        )
        
        df2.to_csv("/shared/projects/origamis/wgs/output/iSMC/"+str(cc)+"_"+str(ss)+"_"+str(window)+"kb.csv", index=False)

for ss in tqdm([samples]):
  for cc in tqdm(chrom,leave=False):
    dd = pd.read_csv("/shared/projects/origamis/wgs/output/iSMC/"+str(cc)+"_"+str(ss)+"_"+str(1000000)+"kb.csv")
    for ww in [5000000]:
      window = range(0+5000000, dd["POS"].tolist()[-1],5000000)
      POS = []
      RHO = []
      for www in window:
        POS += [dd[dd["POS"]<www]["POS"].tolist()[-5:][0]]
        RHO += [np.mean(dd[dd["POS"]<www]["RHO"].tolist()[-5:])]
      
      df2 = pd.DataFrame(
      {
        "POS": POS,
        "RHO": RHO
      }
      )
      df2.to_csv("/shared/projects/origamis/wgs/output/iSMC/"+str(cc)+"_"+str(ss)+"_"+str(ww)+"kb.csv", index=False)

for ss in tqdm([samples]):
  for cc in tqdm(chrom,leave=False):
    dd = pd.read_csv("/shared/projects/origamis/wgs/output/iSMC/"+str(cc)+"_"+str(ss)+"_"+str(1000000)+"kb.csv")
    for ww in [10000000]:
      window = range(0+10000000, dd["POS"].tolist()[-1],10000000)
      POS = []
      RHO = []
      for www in window:
        POS += [dd[dd["POS"]<www]["POS"].tolist()[-5:][0]]
        RHO += [np.mean(dd[dd["POS"]<www]["RHO"].tolist()[-5:])]
      
      df2 = pd.DataFrame(
      {
        "POS": POS,
        "RHO": RHO
      }
      )
      df2.to_csv("/shared/projects/origamis/wgs/output/iSMC/"+str(cc)+"_"+str(ss)+"_"+str(ww)+"kb.csv", index=False)
