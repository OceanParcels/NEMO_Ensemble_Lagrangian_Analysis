# %% Load the packages
import numpy as np
import matplotlib.pyplot as plt
# from tqdm import tqdm
import pandas as pd
import seaborn as sns
# from scipy.stats import gaussian_kde
import pickle
import os
import cmocean.cm as cmo

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

def KLDivergence(P, Q):
    # Kullback-Leibler divergence D_KL(P||Q)
    # P is the true distribution, Q is the estimated distribution.
    return cross_entropy(P, Q) - marginal_entropy(Q)
    
#%% ################# Single Member #################
Latitude_limit = None # 53 or 44
Longitude_limit = -40 # -40

if Latitude_limit is not None:
    criterium_string = f"_{Latitude_limit}N"
elif Longitude_limit is not None:
    criterium_string = f"_{abs(Longitude_limit)}W"

path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/"

N_members = 50

all_messages = {}
all_prob = {}

time_range = np.arange(0, 6*365-20*7, 1)
all_prob['time'] = time_range


for delta_r in [0.1, 1., 2.]:
    
    messages = np.zeros((50, len(time_range)))
    
    for member in range(1, N_members+1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"dr_{delta_r*100:03.0f}_{abs(Longitude_limit)}W/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
            
        ping = np.zeros_like(time_range)
        
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                distributions = pickle.load(f)
            
            drift_time = distributions["drift_time"]
            drift_time = np.sort(drift_time)
            depths = distributions["depths"]
            trajectory = distributions["trajectory"]
            
            unique, counts = np.unique(drift_time, return_counts=True)
            
            for k in time_range:
                if k in unique:
                    ping[k] = 1 # counts[np.where(unique == k)[0][0]]
            
            if np.sum(ping) > 0:
                ping = ping / np.sum(ping)
            # ping = ping / np.sum(ping)
            
        else:
            print(f"Delta_r {delta_r}, member {member} --EMPTY--")
        
        messages[member - 1, :] = ping
    
    all_messages[delta_r] = messages
    all_prob[delta_r] = np.mean(messages, axis=0)/np.sum(np.mean(messages, axis=0))


for week in [4, 12, 20]:
    
    messages = np.zeros((50, len(time_range)))
    
    for member in range(1, N_members+1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"W_{week:02d}_{Latitude_limit}N/Distributions_W{week:02d}_m{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"W_{week:02d}_{abs(Longitude_limit)}W/Distributions_W{week:02d}_m{member:03d}.pkl"
        
        ping = np.zeros_like(time_range)
        
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                distributions = pickle.load(f)
            
            drift_time = distributions["drift_time"]
            drift_time = np.sort(drift_time)
            depths = distributions["depths"]
            trajectory = distributions["trajectory"]
            
            unique, counts = np.unique(drift_time, return_counts=True)
            
            for k in time_range:
                if k in unique:
                    ping[k] =  1# counts[np.where(unique == k)[0][0]]
            
            if np.sum(ping) > 0:
                ping = ping / np.sum(ping)
            # ping = ping / np.sum(ping)
            
        else:
            print(f"Week {week}, member {member} --EMPTY--")


        messages[member - 1, :] = ping
    
    all_messages[week] = messages
    all_prob[week] = np.mean(messages, axis=0)/np.sum(np.mean(messages, axis=0))

    
# %% ################ MIXTURES ################
N_members = 50

all_mix_messages = {}
all_prob_mix = {}
all_prob_mix['time'] = time_range

for delta_r in [0.1, 1., 2.]:
    
    messages = np.zeros((50, len(time_range)))
    
    for member in range(1, N_members+1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"mix_dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_mix_dr{delta_r*100:03.0f}_s{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"mix_dr_{delta_r*100:03.0f}_{abs(Longitude_limit)}W/Distributions_mix_dr{delta_r*100:03.0f}_s{member:03d}.pkl"

        ping = np.zeros_like(time_range)
        
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                distributions = pickle.load(f)
            
            drift_time = distributions["drift_time"]
            drift_time = np.sort(drift_time)
            depths = distributions["depths"]
            trajectory = distributions["trajectory"]
            
            unique, counts = np.unique(drift_time, return_counts=True)
            
            for k in time_range:
                if k in unique:
                    ping[k] = 1 # counts[np.where(unique == k)[0][0]]
            
            if np.sum(ping) > 0:
                ping = ping / np.sum(ping)
            # ping = ping / np.sum(ping)
            
        else:
            print(f"Delta_r {delta_r}, member {member} --EMPTY--")


        messages[member - 1, :] = ping 
    
    all_mix_messages[delta_r] = messages
    all_prob_mix[delta_r] = np.mean(messages, axis=0)/np.sum(np.mean(messages, axis=0))
    
    
for week in [4, 12, 20]:
    
    messages = np.zeros((50, len(time_range)))
    
    for member in range(1, N_members+1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"mix_W_{week:02d}_{Latitude_limit}N/Distributions_mix_W{week:02d}_s{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"mix_W_{week:02d}_{abs(Longitude_limit)}W/Distributions_mix_W{week:02d}_s{member:03d}.pkl"
            
        ping = np.zeros_like(time_range)
        
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                distributions = pickle.load(f)
            
            drift_time = distributions["drift_time"]
            drift_time = np.sort(drift_time)
            depths = distributions["depths"]
            trajectory = distributions["trajectory"]
            
            unique, counts = np.unique(drift_time, return_counts=True)
            
            for k in time_range:
                if k in unique:
                    ping[k] = 1 # counts[np.where(unique == k)[0][0]]
            
            if np.sum(ping) > 0:
                ping = ping / np.sum(ping)
            # ping = ping / np.sum(ping)
            
        else:
            print(f"Week {week}, member {member} --EMPTY--")

        messages[member - 1, :] = ping
    
    all_mix_messages[week] = messages
    all_prob_mix[week] = np.mean(messages, axis=0)
    
#%% Check that the sum of the probabilities is 1
for key in all_prob.keys():
    print(f"Sum of probabilities for {key}: {np.sum(all_prob[key])}")
    
for key in all_prob_mix.keys():
    print(f"Sum of probabilities for {key}: {np.sum(all_prob_mix[key])}")
    

# %% ############### PLOT Single members and Mixtures ####################
fig, axs = plt.subplots(2, 2, figsize=(8, 6), sharex=True, sharey=True)
colors_space = ['mediumblue', 'blueviolet', 'teal']
colors_time = ['darkred', 'orangered', 'orange']

# Top left: Single members (delta_r)
i = 0
for delta_r in [0.1, 1., 2.]:
    axs[0, 0].plot(all_prob[delta_r], label=f"$\delta_r = {delta_r}^o$", alpha=0.5, color=colors_space[i], lw=0.3)
    i += 1
    
axs[0, 0].set_ylabel("Probability")
axs[0, 0].legend(fontsize="small")
axs[0, 0].text(0.05, 0.05, 'A', transform=axs[0, 0].transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='left', fontweight='bold')

# Bottom left: Single members (week)
i = 0
for week in [4, 12, 20]:
    axs[1, 0].plot(all_prob[week], label=f"Week {week}", alpha=0.5, color=colors_time[i], lw=0.3)
    i += 1
axs[1, 0].set_xlabel("Particle Age (days)")
axs[1, 0].set_ylabel("Probability")
axs[1, 0].legend(fontsize='small')
axs[1, 0].text(0.05, 0.05, 'C', transform=axs[1, 0].transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='left', fontweight='bold')

# Top right: Mixtures (delta_r)
i = 0
for delta_r in [0.1, 1., 2.]:
    axs[0, 1].plot(all_prob_mix[delta_r], label=f"Mix. $\delta_r = {delta_r}^o$", alpha=0.5, color=colors_space[i], lw=0.3)
    i += 1
axs[0, 1].legend(fontsize='small')
axs[0, 1].text(0.05, 0.05, 'B', transform=axs[0, 1].transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='left', fontweight='bold')

# Bottom right: Mixtures (week)
i = 0
for week in [4, 12, 20]:
    axs[1, 1].plot(all_prob_mix[week], label=f"Mix. Week {week}", alpha=0.5, color=colors_time[i], lw=0.3)
    i += 1
    
axs[1, 1].set_xlabel("Particle Age (days)")
axs[1, 1].legend(fontsize='small')
axs[1, 1].text(0.05, 0.05, 'D', transform=axs[1, 1].transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='left', fontweight='bold')

# Adjust layout
plt.tight_layout()
plt.show()
fig.savefig("../figs/FigX_Frequency_probability" + criterium_string + ".png", dpi=300)

# %% Dataframe with the statistics
ensemble = pd.DataFrame(all_prob)
ensemble.set_index('time', inplace=True)
# ensemble.columns = [r'$\delta_r 0.1^o$', 'dr 1.', 'dr 2.', '4 weeks', '12 weeks', '20 weeks']

mixture = pd.DataFrame(all_prob_mix)
mixture.set_index('time', inplace=True)

#%% Comparisson mesage by message
ensemble_rel_H = {}
df_labels = ['$\\delta_r = 0.1^o$', '$\\delta_r = 1^o$', '$\\delta_r = 2^o$', '4 weeks', '12 weeks', '20 weeks']


for j, delta in enumerate([0.1, 1., 2., 4, 12, 20]):
    
    KLD = np.zeros(6)
    
    for i, delta_ref in enumerate([0.1, 1., 2., 4, 12, 20]):
        
        relative_entropy = []
        
        for member in range(0, N_members):
            for set in range(0, N_members):
                
                full = all_mix_messages[delta][set, :]
                aprox = all_messages[delta_ref][member, :]
        
                relative_entropy.append(KLDivergence(aprox, full))

        KLD[i] = np.mean(relative_entropy)
    ensemble_rel_H[df_labels[j]] = KLD


# %% DATAFRAMES Entropy of the ensemble and mixture
ensemble_cross_entropy = {}
ensemble_KLD = {}
ensemble_entropy = np.zeros(6)
ensemble_mix_entropy = np.zeros(6)


for j, delta in enumerate([0.1, 1., 2., 4, 12, 20]):
    
    cross = np.zeros(6)
    KLD = np.zeros(6)
    
    for i, delta_ref in enumerate([0.1, 1., 2., 4, 12, 20]):
        
        full = all_prob_mix[delta]
        aprox = all_prob[delta_ref]
        
        # cross[i] = cross_entropy(all_prob_mix[delta_ref], all_prob[delta]) #- marginal_entropy(all_prob[delta])
        # KLD[i] = cross_entropy(all_prob_mix[delta_ref], all_prob[delta]) - marginal_entropy(all_prob[delta])
        
        cross[i] = cross_entropy(aprox, full) #- marginal_entropy(all_prob[delta])
        KLD[i] = KLDivergence(aprox, full)

    ensemble_cross_entropy[delta] = cross
    ensemble_KLD[delta] = KLD
    
    ensemble_entropy[j] = marginal_entropy(all_prob[delta])
    ensemble_mix_entropy[j] = marginal_entropy(all_prob_mix[delta])
    

#redefine keys in ensemble_cross_KLD
ensemble_KLD = {r'$\delta_r = 0.1^o$': ensemble_KLD[0.1],
                            r'$\delta_r = 1^o$': ensemble_KLD[1.],
                            r'$\delta_r = 2^o$': ensemble_KLD[2.],
                            '4 weeks': ensemble_KLD[4],
                            '12 weeks': ensemble_KLD[12],
                            '20 weeks': ensemble_KLD[20]}


H_cross = pd.DataFrame(ensemble_cross_entropy, index=[r'Mix $\delta_r = 0.1^o$', r'Mix $\delta_r = 1^o$', r'Mix $\delta_r = 2^o$', 
                                                           'Mix 4 weeks', 'Mix 12 weeks', 'Mix 20 weeks'])
KLD_ensemble = pd.DataFrame(ensemble_KLD, index=[r'Mix $\delta_r = 0.1^o$', r'Mix $\delta_r = 1^o$', r'Mix $\delta_r = 2^o$', 
                                                           'Mix 4 weeks', 'Mix 12 weeks', 'Mix 20 weeks'])

H_ensemble = pd.DataFrame(ensemble_entropy, index=[r'$\delta_r = 0.1^o$', r'$\delta_r = 1^o$', r'$\delta_r = 2^o$',
                                                   '4 weeks', '12 weeks', '20 weeks'])
H_ensemble_mix = pd.DataFrame(ensemble_mix_entropy, index=[r'Mix $\delta_r = 0.1^o$', r'Mix $\delta_r = 1^o$', r'Mix $\delta_r = 2^o$', 
                                                           'Mix 4 weeks', 'Mix 12 weeks', 'Mix 20 weeks'])


# %% PLOT CROSS ENTROPY 
# Create a 2x2 grid of subplots with space for colorbars
fig, axs = plt.subplots(2, 3, figsize=(8, 7), gridspec_kw={'width_ratios': [1/7, 6/7, 1/7], 'height_ratios': [6/7, 1/7]})

cmmap = 'Greys_r'
central_plot_data = H_cross.T
Max_ent = np.max(central_plot_data.max())
Min_ent = np.min(central_plot_data.min())

# Marginal entropy single members (upper left)
sns.heatmap(H_ensemble, annot=True, fmt=".3f", cmap=cmmap, ax=axs[0, 0], cbar=False, vmin=Min_ent, vmax=Max_ent) 
# axs[0, 0].set_title('Marginal Entropy')
axs[0, 0].set_xticklabels([r'$H(P_i)$'])  # Remove x tick labels
# axs[0, 0].set_xticks([])  # Remove x ticks

# Marginal entropy mixtures (lower right)
sns.heatmap(H_ensemble_mix.T, annot=True, fmt=".3f", cmap=cmmap, ax=axs[1, 1], cbar=False, vmin=Min_ent, vmax=Max_ent)
# axs[1, 1].set_title('Marginal Entropy')
axs[1, 1].set_yticklabels([r'$H(P_{Mix})$'])  # Remove y tick labels
# axs[1, 1].set_yticks([])  # Remove y ticks
axs[1, 1].set_xticklabels(axs[1, 1].get_xticklabels(), rotation=15, ha='center')  # Rotate x tick labels

# Hide the bottom left subplot (lower left)
axs[1, 0].text(0.5, 0.5, "Marginal\nEntropy", ha='center', va='center', rotation=-45, fontsize=12)
axs[1, 0].axis('off')
axs[1, 2].axis('off')
axs[0, 2].axis('off')
axs[0, 1].axis('off')

# Add colorbars
fig.subplots_adjust(right=0.85)
cbar_ax2 = fig.add_axes([0.87, 0.28, 0.03, 0.66])

sns.heatmap(central_plot_data, annot=True, fmt=".3f", cmap=cmmap, ax=axs[0, 1], cbar_ax=cbar_ax2)
cbar_ax2.set_ylabel('(bits)')
axs[0, 1].set_title(r'Cross Entropy, $H_{P_{i}}(P_{Mix})$')

# Adjust layout
plt.tight_layout()
plt.show()

fig.savefig("../figs/FigX_Cross_entropy" + criterium_string + ".png", dpi=300)

# %% Kullback-Leibler divergence plot 
fig, ax = plt.subplots(figsize=(6, 5))

cmmap = 'Blues_r' #cmo.algae_r

sns.heatmap(KLD_ensemble.T, annot=True, fmt=".3f", cmap=cmmap, ax=ax, cbar=True, vmin=0)

# Rotate y tick labels 90 degrees
ax.set_yticklabels(ax.get_yticklabels(), rotation=90, ha='center', fontsize=9)
# Rotate x tick labels 15 degrees
ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha='center', fontsize=9)

# Add colorbar label
cbar = ax.collections[0].colorbar
cbar.set_label('Relative Entropy, $D(P_{Mix}||P_i)$ (bits)')

# ax.set_xlabel(r'$P_{Mix}$')

plt.tight_layout()
plt.show()

fig.savefig("../figs/FigX_KLDivergence" + criterium_string + ".png", dpi=300)

# %%
