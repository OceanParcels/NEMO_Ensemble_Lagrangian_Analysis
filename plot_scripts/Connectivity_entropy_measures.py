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

# # %% Build KDEs to compute entropy etc.. 
# def evaluate_kde(kde, x):
#     pdf = kde(x)
#     pdf = np.where(pdf > 0, pdf, np.nan)
#     return pdf

# x = np.linspace(0, 7500, 7500)

# kde_temp = {}
# kde_space = {}
# kde_mix_temp = {}
# kde_mix_space = {}

# for week in [4, 12, 20]:
#     kde_single = gaussian_kde(all_temp[week]["counts"])
#     kde_mix = gaussian_kde(all_mix_temp[week]["counts"])
    
#     kde_temp[week] = evaluate_kde(kde_single, x)
#     kde_mix_temp[week] = evaluate_kde(kde_mix, x)
    
# for delta_r in [0.1, 1., 2.]:
#     kde_single = gaussian_kde(all_space[delta_r]["counts"])
#     kde_mix = gaussian_kde(all_mix_space[delta_r]["counts"])
#     kde_space[delta_r] = evaluate_kde(kde_single, x)
#     kde_mix_space[delta_r] = evaluate_kde(kde_mix, x)
    

# %% Compute the entropy of the distributions
# ensemble_entropy = np.zeros(6)
# mixture_entropy = np.zeros(6)

# i = 0
# for week in [4, 12, 20]:
#     ensemble_entropy[i] = marginal_entropy(kde_temp[week])
#     mixture_entropy[i] = marginal_entropy(kde_mix_temp[week])
#     i += 1

# for delta_r in [0.1, 1., 2.]:
#     ensemble_entropy[i] = marginal_entropy(kde_space[delta_r])
#     mixture_entropy[i] = marginal_entropy(kde_mix_space[delta_r])
#     i += 1

# np.where(pdf > 0, pdf, np.finfo(float).eps)

#%%
Latitude_limit = None # 44 or 53
Longitude_limit = -40
path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/"

N_members = 50

all_messages = {}
all_prob = {}

time_range = np.arange(0, 6*365-1, 1)
all_prob['time'] = time_range


for delta_r in [0.1, 1., 2.]:
    
    messages = np.zeros((50, len(time_range)))
    
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
                sms[k] = 1 # counts[np.where(unique == k)[0][0]]

        sms = sms / np.sum(sms)
        messages[member - 1, :] = sms
    
    all_messages[delta_r] = messages
    all_prob[delta_r] = np.mean(messages, axis=0)

# %%
for delta_r in [0.1, 1., 2.]:
    
    plt.plot(all_prob[delta_r], label=f"$\delta_r = {delta_r}^o$",alpha=1)
    print(sum(all_prob[delta_r]))
    print(f"Entropy dr = {delta_r}: {all_prob[delta_r]}")

plt.xlabel("Time (days)")
plt.ylim(0, 0.0013)
plt.ylabel("Probability")
# plt.title("Probability of particles crossing the 40W meridian")

plt.legend()    
# %%

Latitude_limit = None # 44 or 53
Longitude_limit = -40
path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/"


N_members = 50

all_mix_messages = {}
all_prob_mix = {}
all_prob_mix['time'] = time_range
# all_mix_entropy = {}

for delta_r in [0.1, 1., 2.]:
    
    messages = np.zeros((50, len(time_range)))
    # entropy = np.zeros(50)
    
    for member in range(1, N_members+1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"analysis/connectivity/mix_dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_mix_dr{delta_r*100:03.0f}_s{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"analysis/connectivity/mix_dr_{delta_r*100:03.0f}_{abs(Longitude_limit)}W/Distributions_mix_dr{delta_r*100:03.0f}_s{member:03d}.pkl"

            
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                distributions = pickle.load(f)
            
            drift_time = distributions["drift_time"]
            depths = distributions["depths"]
            trajectory = distributions["trajectory"]
            
        else:
            print(f"--EMPTY--")
            # print("member", member, "dr", delta_r)

        drift_time = np.sort(drift_time)
        # print(drift_time)

        unique, counts = np.unique(drift_time, return_counts=True)

        sms = np.zeros_like(time_range)

        for k in time_range:
            if k in unique:
                sms[k] = 1 #counts[np.where(unique == k)[0][0]]

        # sms = sms / np.sum(sms)
        messages[member - 1, :] = sms 
        # entropy[member - 1] = marginal_entropy(sms)
    
    all_mix_messages[delta_r] = messages
    all_prob_mix[delta_r] = np.mean(messages, axis=0)/np.sum(np.mean(messages, axis=0))
    # all_mix_entropy[delta_r] = entropy

# %%
for delta_r in [0.1, 1., 2.]:
   
    print(sum(all_prob_mix[delta_r]))
    plt.plot(all_prob_mix[delta_r], label=f"Mix. $\delta_r = {delta_r}^o$", alpha=1)
    
    print(f"Entropy dr = {delta_r}: {marginal_entropy(all_prob_mix[delta_r])}")

plt.xlabel("Time (days)")
plt.ylabel("Probability")
plt.ylim(0, 0.0013)
plt.legend()

# %% Dataframe with the statistics
ensemble = pd.DataFrame(all_prob)
ensemble.set_index('time', inplace=True)

mixture = pd.DataFrame(all_prob_mix)
mixture.set_index('time', inplace=True)

# %% Entropy of the ensemble and mixture
ensemble_cross_entropy = {}
ensemble_KLD = {}
ensemble_entropy = np.zeros(3)

for j, delta in enumerate([0.1, 1., 2.]):
    cross = np.zeros(3)
    KLD = np.zeros(3)
    for i, delta_ref in enumerate([0.1, 1., 2.]):
        cross[i] = cross_entropy(all_prob_mix[delta_ref], all_prob[delta]) #- marginal_entropy(all_prob[delta])
        KLD[i] = cross_entropy(all_prob[delta_ref], all_prob_mix[delta]) - marginal_entropy(all_prob[delta])

    ensemble_cross_entropy[delta] = cross
    ensemble_KLD[delta] = KLD
    
    ensemble_entropy[j] = marginal_entropy(all_prob[delta])
    
H_cross = pd.DataFrame(ensemble_cross_entropy, index=['Mix 0.1', 'Mix 1.', 'Mix 2.'])
KLD_ensemble = pd.DataFrame(ensemble_KLD, index=['Mix 0.1', 'Mix 1.', 'Mix 2.'])

H_ensemble = pd.DataFrame(ensemble_entropy, index=['0.1', '1.', '2.'])


# %%
fig, ax = plt.subplots(figsize=(4, 4))
sns.heatmap(H_cross.T, annot=True, fmt=".3f", cmap='magma', ax=ax, vmin=10.64, vmax=11.24, cbar_kws={'label': '(bits)'})
ax.set_title('Cross Entropy')
plt.show()
  
# %%
fig, ax = plt.subplots(figsize=(2, 4))
sns.heatmap(H_ensemble, annot=True, fmt=".3f", cmap='magma', ax=ax, vmin=10.64, vmax=11.24, cbar_kws={'label': '(bits)'})
ax.set_title('Marginal Entropy')
ax.set_xticklabels([])  # Remove x tick labels
plt.show()

# %%
fig, ax = plt.subplots(figsize=(4, 4))
sns.heatmap(KLD_ensemble.T, annot=True, fmt=".3f", cmap='magma', ax=ax, cbar_kws={'label': '(bits)'})
ax.set_title('Kullback Leibler Divergence')
plt.show()

# %%
