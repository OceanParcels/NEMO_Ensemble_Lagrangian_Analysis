#%%
import numpy as np
from glob import glob
import xarray as xr
from tqdm import tqdm
import os
import re
import pandas as pd

fname_general = '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/NATL025-CJMCYC3.001-S/1d/2010/NATL025-CJMCYC3.001_y2010m*.1d_gridU.nc'

#%% Check if the .npy is already created
savefile = 'check_corrupted_files.npy'

if os.path.isfile(savefile):
    checked_files = np.load(savefile).tolist()
else:
    checked_files = []

# Check if the .npy is already created. If not, create an empty list
# Corrupted files provides a list of the corrupted files
corrupted_files = []


# %%
variables = ['U', 'V', 'W']

for member in range(1, 51):
    print(f"Member {member:03d}")

    for year in range(2010, 2016):
        print(f"Year {year}")
        fname_general = f"/storage/shared/oceanparcels/input_data/NEMO_Ensemble/NATL025-CJMCYC3.{member:03d}-S/1d/{year}/NATL025-CJMCYC3.{member:03d}_y{year}m*.1d_gridU.nc"

        for fname in variables:
            print(fname)
            files = sorted(glob(fname_general.replace('_gridU', f'_grid{fname}')))
            I = {} # dictionary to store the NaNs
            for f in tqdm(files):
                if f not in checked_files: # check if the file has already been checked
                    # some files cannot be opened, so we need to use a try-except
                    try:
                        ds = xr.open_dataset(f, decode_cf=False) # open the file
                        
                        for vname, v in ds.items(): # loop over the variables
                                
                            try:
                                v = np.array(v) # convert to numpy array
                                if vname not in I: 
                                    I[vname] = np.argwhere(np.isnan(v))
                                else:  # check if NaNs at exactly the same positions
                                    try:
                                        I2 = np.argwhere(np.isnan(v)) # get the NaNs
                                        assert np.all(I[vname] == I2) # check if they are the same
                                        checked_files.append(f) # if they are the same, add to the list of checked files
                                        np.save(savefile, checked_files)
                                    except:
                                        print(f'file is corrupted (not same NaNs):  {f}') # if they are not the same, print the file
                                        corrupted_files.append(f) # add to the list of corrupted files
                        
                            except:
                                print(f'file is corrupted (unable to open Variable):  {f}') # if the file is corrupted, print the file
                                corrupted_files.append(f)
                                
                    except:
                        print(f'file is corrupted (unable to open Dataset):  {f}')
                        corrupted_files.append(f)
                        

np.save('corrupted_files.npy', corrupted_files) # save the list of corrupted files

# %% Analize the corrupted files
# read file line by line and extract path
# with open("NEMO_file_check.35524.o", "r") as f:
#     corrupted_files = f.readlines()

# # look for string with 'file is corrupted' and extract the path
# corrupted_files = [
#     re.findall(r"(/storage.*)", x)[0]
#     for x in corrupted_files
#     if "file is corrupted" in x
# ]

corrupted_files = np.load("corrupted_files.npy")
corrupted_files = np.unique(corrupted_files)

# remove the '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/' from the path
corrupted_files = [x.replace("/storage/shared/oceanparcels/input_data/NEMO_Ensemble/", "") for x in corrupted_files]

with open("corrupted_files_NEMO_Ensemble.txt", "w") as f:
    for item in corrupted_files:
        f.write("%s\n" % item)

# %% From corrupted files, extract the year and member and month and save it in a dataframe

corrupted_files_df = pd.DataFrame(corrupted_files, columns=["path"])

corrupted_files_df["member"] = corrupted_files_df["path"].apply(
    lambda x: int(re.findall(r"NATL025-CJMCYC3.(\d{3})-S", x)[0])
)
corrupted_files_df["year"] = corrupted_files_df["path"].apply(
    lambda x: int(re.findall(r"y(\d{4})m\d{2}", x)[0])
)
corrupted_files_df["month"] = corrupted_files_df["path"].apply(
    lambda x: int(re.findall(r"y\d{4}m(\d{2})", x)[0])
)
corrupted_files_df["variable"] = corrupted_files_df["path"].apply(
    lambda x: re.findall(r"grid(\w).nc", x)[0]
)

# save the dataframe as an .csv
corrupted_files_df.to_csv("corrupted_files_NEMO_Ensemble.csv", index=False)

# %% filter for corrupted files in 2015 
corrupted_files_2015 = corrupted_files_df[corrupted_files_df["year"] == 2015]

# %%
