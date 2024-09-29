# %% Load the packages
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde
import pickle
import os

def marginal_entropy(P):
    # Shannon entropy
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    P_safe = np.where(P > 0, P, np.finfo(float).eps)
    return -np.nansum(P_safe * np.log2(P_safe))

def cross_entropy(P, Q):
    # Cross entropy H_p(q) = - sum_x p(x) log q(x)
    # P is the true distribution, Q is the estimated distribution.
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    P_safe = np.where(P > 0, P, np.finfo(float).eps)
    return -np.nansum(Q * np.log2(P_safe))


# %% Load all dataframes

all_mix_temp = {}
all_mix_space = {}
all_temp = {}
all_space = {}
location = "Cape_Hatteras"
base_path = "/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/"


Latitude_limit = 53 # 44 or 53
Longitude_limit = None # -40 

if Latitude_limit is not None:
    criterium_string = f"_{Latitude_limit}N"
elif Longitude_limit is not None:
    criterium_string = f"_{abs(Longitude_limit)}W"

for week in [4, 12, 20]:
    
    df_path = base_path + f"analysis/connectivity/Stats/Stats_W{week:02d}" + criterium_string + ".csv"
    _df = pd.read_csv(df_path)
    all_temp[week] = _df

    df_mix = base_path + f"analysis/connectivity/Stats/Stats_mix_W{week:02d}" + criterium_string + ".csv"
    _df = pd.read_csv(df_mix)
    all_mix_temp[week] = _df
    
for delta_r in [0.1, 1., 2.]:
    df_path = base_path + f"analysis/connectivity/Stats/Stats_dr{delta_r*100:03.0f}" + criterium_string + ".csv"
    _df = pd.read_csv(df_path)
    all_space[delta_r] = _df
    
    df_path = base_path + f"analysis/connectivity/Stats/Stats_mix_dr{delta_r*100:03.0f}" + criterium_string + ".csv"
    _df = pd.read_csv(df_path)
    all_mix_space[delta_r] = _df

# %% Build KDEs to compute entropy etc.. 
def evaluate_kde(kde, x):
    pdf = kde(x)
    pdf = np.where(pdf > 0, pdf, np.nan)
    return pdf

x = np.linspace(0, 7500, 7500)

kde_temp = {}
kde_space = {}
kde_mix_temp = {}
kde_mix_space = {}

for week in [4, 12, 20]:
    kde_single = gaussian_kde(all_temp[week]["counts"])
    kde_mix = gaussian_kde(all_mix_temp[week]["counts"])
    
    kde_temp[week] = evaluate_kde(kde_single, x)
    kde_mix_temp[week] = evaluate_kde(kde_mix, x)
    
for delta_r in [0.1, 1., 2.]:
    kde_single = gaussian_kde(all_space[delta_r]["counts"])
    kde_mix = gaussian_kde(all_mix_space[delta_r]["counts"])
    kde_space[delta_r] = evaluate_kde(kde_single, x)
    kde_mix_space[delta_r] = evaluate_kde(kde_mix, x)
    

# %% Compute the entropy of the distributions
ensemble_entropy = np.zeros(6)
mixture_entropy = np.zeros(6)

i = 0
for week in [4, 12, 20]:
    ensemble_entropy[i] = marginal_entropy(kde_temp[week])
    mixture_entropy[i] = marginal_entropy(kde_mix_temp[week])
    i += 1

for delta_r in [0.1, 1., 2.]:
    ensemble_entropy[i] = marginal_entropy(kde_space[delta_r])
    mixture_entropy[i] = marginal_entropy(kde_mix_space[delta_r])
    i += 1

# np.where(pdf > 0, pdf, np.finfo(float).eps)

#%%
Latitude_limit = 53 # 44 or 53
Longitude_limit = None
path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/"


N_members = 50

all_messages = {}
all_entropy = {}
time_range = np.arange(0, 6*365-1, 1)

for delta_r in [0.1, 1., 2.]:
    
    messages = np.zeros((50, len(time_range)))
    entropy = np.zeros(50)
    
    for member in range(1, N_members+1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"analysis/connectivity/dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"analysis/connectivity/dr_{delta_r*100:03.0f}_{abs(Longitude_limit)}W/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
            
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                distributions = pickle.load(f)
            
            drift_time = distributions["drift_time"]
            depths = distributions["depths"]
            trajectory = distributions["trajectory"]
            
        else:
            print(f"--EMPTY--")

        drift_time = np.sort(drift_time)

        unique, counts = np.unique(drift_time, return_counts=True)

        sms = np.zeros_like(time_range)

        for k in time_range:
            if k in unique:
                sms[k] = 1 #counts[np.where(unique == k)[0][0]]

        sms = sms / np.sum(sms)
        messages[member - 1, :] = sms 
        entropy[member - 1] = marginal_entropy(sms)
    
    all_messages[delta_r] = messages
    all_entropy[delta_r] = entropy

# %%
for delta_r in [0.1, 1., 2.]:
    average_message = np.mean(all_messages[delta_r], axis=0)
    print(sum(average_message))
    plt.plot(average_message, label=f"dr = {delta_r}")
    
    print(f"Entropy dr = {delta_r}: {marginal_entropy(average_message)}")

plt.legend()    
# %%

# %% Load the packages
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde
import pickle
import os

def marginal_entropy(P):
    # Shannon entropy
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    P_safe = np.where(P > 0, P, np.finfo(float).eps)
    return -np.nansum(P_safe * np.log2(P_safe))

def cross_entropy(P, Q):
    # Cross entropy H_p(q) = - sum_x p(x) log q(x)
    # P is the true distribution, Q is the estimated distribution.
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    P_safe = np.where(P > 0, P, np.finfo(float).eps)
    return -np.nansum(Q * np.log2(P_safe))


# %% Load all dataframes

all_mix_temp = {}
all_mix_space = {}
all_temp = {}
all_space = {}
location = "Cape_Hatteras"
base_path = "/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/"


Latitude_limit = 53 # 44 or 53
Longitude_limit = None # -40 

if Latitude_limit is not None:
    criterium_string = f"_{Latitude_limit}N"
elif Longitude_limit is not None:
    criterium_string = f"_{abs(Longitude_limit)}W"

for week in [4, 12, 20]:
    
    df_path = base_path + f"analysis/connectivity/Stats/Stats_W{week:02d}" + criterium_string + ".csv"
    _df = pd.read_csv(df_path)
    all_temp[week] = _df

    df_mix = base_path + f"analysis/connectivity/Stats/Stats_mix_W{week:02d}" + criterium_string + ".csv"
    _df = pd.read_csv(df_mix)
    all_mix_temp[week] = _df
    
for delta_r in [0.1, 1., 2.]:
    df_path = base_path + f"analysis/connectivity/Stats/Stats_dr{delta_r*100:03.0f}" + criterium_string + ".csv"
    _df = pd.read_csv(df_path)
    all_space[delta_r] = _df
    
    df_path = base_path + f"analysis/connectivity/Stats/Stats_mix_dr{delta_r*100:03.0f}" + criterium_string + ".csv"
    _df = pd.read_csv(df_path)
    all_mix_space[delta_r] = _df

# %% Build KDEs to compute entropy etc.. 
def evaluate_kde(kde, x):
    pdf = kde(x)
    pdf = np.where(pdf > 0, pdf, np.nan)
    return pdf

x = np.linspace(0, 7500, 7500)

kde_temp = {}
kde_space = {}
kde_mix_temp = {}
kde_mix_space = {}

for week in [4, 12, 20]:
    kde_single = gaussian_kde(all_temp[week]["counts"])
    kde_mix = gaussian_kde(all_mix_temp[week]["counts"])
    
    kde_temp[week] = evaluate_kde(kde_single, x)
    kde_mix_temp[week] = evaluate_kde(kde_mix, x)
    
for delta_r in [0.1, 1., 2.]:
    kde_single = gaussian_kde(all_space[delta_r]["counts"])
    kde_mix = gaussian_kde(all_mix_space[delta_r]["counts"])
    kde_space[delta_r] = evaluate_kde(kde_single, x)
    kde_mix_space[delta_r] = evaluate_kde(kde_mix, x)
    

# %% Compute the entropy of the distributions
ensemble_entropy = np.zeros(6)
mixture_entropy = np.zeros(6)

i = 0
for week in [4, 12, 20]:
    ensemble_entropy[i] = marginal_entropy(kde_temp[week])
    mixture_entropy[i] = marginal_entropy(kde_mix_temp[week])
    i += 1

for delta_r in [0.1, 1., 2.]:
    ensemble_entropy[i] = marginal_entropy(kde_space[delta_r])
    mixture_entropy[i] = marginal_entropy(kde_mix_space[delta_r])
    i += 1

# np.where(pdf > 0, pdf, np.finfo(float).eps)

#%%
Latitude_limit = None # 44 or 53
Longitude_limit = -40
path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/"


N_members = 50

all_messages = {}
all_entropy = {}
time_range = np.arange(0, 6*365-1, 1)

for delta_r in [0.1, 1., 2.]:
    
    messages = np.zeros((50, len(time_range)))
    entropy = np.zeros(50)
    
    for member in range(1, N_members+1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"analysis/connectivity/dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"analysis/connectivity/dr_{delta_r*100:03.0f}_{abs(Longitude_limit)}W/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
            
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                distributions = pickle.load(f)
            
            drift_time = distributions["drift_time"]
            depths = distributions["depths"]
            trajectory = distributions["trajectory"]
            
        else:
            print(f"--EMPTY--")

        drift_time = np.sort(drift_time)

        unique, counts = np.unique(drift_time, return_counts=True)

        sms = np.zeros_like(time_range)

        for k in time_range:
            if k in unique:
                sms[k] = 1 #counts[np.where(unique == k)[0][0]]

        # sms = sms / np.sum(sms)
        messages[member - 1, :] = sms 
        entropy[member - 1] = marginal_entropy(sms)
    
    all_messages[delta_r] = messages
    all_entropy[delta_r] = entropy

# %%
P_frequency = np.zeros((6, len(time_range)))

i = 0
for delta_r in [0.1, 1., 2.]:
    average_message = np.mean(all_messages[delta_r], axis=0)
    average_message = average_message / np.sum(average_message)
    P_frequency[i] = average_message
    i += 1
    print(sum(average_message))
    plt.plot(average_message, label=f"dr = {delta_r}")
    
    print(f"Entropy dr = {delta_r}: {marginal_entropy(average_message)}")

plt.legend()    
# %% ########## MIXTURE DISTRIBUTIONS ############
###############################################

Latitude_limit = 53 # 44 or 53
Longitude_limit = None
path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/"


N_members = 50

all_messages = {}
all_entropy = {}
time_range = np.arange(0, 6*365-1, 1)

for delta_r in [0.1, 1., 2.]:
    
    messages = np.zeros((50, len(time_range)))
    entropy = np.zeros(50)
    
    for member in range(1, N_members+1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"analysis/connectivity/mix_dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_mix_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"analysis/connectivity/mix_dr_{delta_r*100:03.0f}_{abs(Longitude_limit)}W/Distributions_mix_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
            
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                distributions = pickle.load(f)
            
            drift_time = distributions["drift_time"]
            depths = distributions["depths"]
            trajectory = distributions["trajectory"]
            
        else:
            print(f"--EMPTY--")

        drift_time = np.sort(drift_time)

        unique, counts = np.unique(drift_time, return_counts=True)

        sms = np.zeros_like(time_range)

        for k in time_range:
            if k in unique:
                sms[k] = 1 #counts[np.where(unique == k)[0][0]]

        sms = sms / np.sum(sms)
        messages[member - 1, :] = sms 
        entropy[member - 1] = marginal_entropy(sms)
    
    all_messages[delta_r] = messages
    all_entropy[delta_r] = entropy

