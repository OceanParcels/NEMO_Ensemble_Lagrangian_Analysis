# %%
import numpy as np
from glob import glob
import xarray as xr
from tqdm import tqdm
import os
import re
import pandas as pd

fname_general = "/storage/shared/oceanparcels/input_data/NEMO_Ensemble/NATL025-CJMCYC3.001-S/1d/2010/NATL025-CJMCYC3.001_y2010m01.1d_gridU.nc"

# look for character 'm' foloowed by two numbers in the string
# re.findall(r'm\d{2}', fname_general)

# %% Check if the .npy is already created
savefile = "missing_NEMO_ensemble_files.npy"

missing_files = []

# %%
variables = ["U", "V", "W"]

for member in range(1, 51):
    # print(f"Member {member:03d}")
    for year in range(2010, 2016):
        for var in variables:
            for mm in range(1, 13):
                file = f"/storage/shared/oceanparcels/input_data/NEMO_Ensemble/NATL025-CJMCYC3.{member:03d}-S/1d/{year}/NATL025-CJMCYC3.{member:03d}_y{year}m{mm:02d}.1d_grid{var}.nc"

                if not os.path.isfile(file):
                    missing_files.append(file)

missing_files = [x.replace("/storage/shared/oceanparcels/input_data/NEMO_Ensemble/", "") for x in missing_files]
# %% Test missing files actually don't exist
existing_files = 0
not_existing_files = 0

for f in tqdm(missing_files):
    if not os.path.isfile(f):
        print(f"File {f} does not exist")
        not_existing_files += 1
    else:
        existing_files += 1
# %% export missing files to .txt file with each element in new row
with open("missing_files_NEMO_Ensemble.txt", "w") as f:
    for item in missing_files:
        f.write("%s\n" % item)

# %%
