# %% Load the packages
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import seaborn as sns

# %% Load all dataframes

all_mix_temp = {}
all_mix_space = {}
all_temp = {}
all_space = {}
location = "Cape_Hatteras"
base_path = "/Volumes/Claudio SSD/Ensemble_article_data/"


Latitude_limit = None #44 # 44 or 53
Longitude_limit =  -40

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

K_h = 10
df_path = base_path + f"analysis/connectivity/Stats/Stats_Kh_{K_h:01d}" + criterium_string + ".csv"
all_diff = pd.read_csv(df_path)


#%% Plot the percentage of subpolar trajectories for temporal and spatial members
fig, ax = plt.subplots(2, 3, figsize=(10, 6))

ax = ax.flatten()
# Plot distributions for percentage of subpolar trajectories

colors_space = ["midnightblue", "blueviolet", "teal"]
ls_space = [(0, (1, 1)), '--', '-.']
j = 0

#spatial
for delta_r in [0.1, 1., 2.]:
    sns.kdeplot(all_space[delta_r]["counts"], ax=ax[0], label=f"$\delta_r = {delta_r}^o$", clip=(0, 7500),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    # kde_space[delta_r] = _counts_dr.get_lines()[j].get_data() # save the KDE values for entropy computation
    
    sns.kdeplot(all_space[delta_r]["median_time"]/365, ax=ax[1], label=f"Spatial dr{delta_r}", clip=(0, 6),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    sns.kdeplot(all_space[delta_r]["median_depth"], ax=ax[2], label=f"Spatial dr{delta_r}", clip=(0, 1000),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    j += 1

#diffusion
sns.kdeplot(all_diff["counts"], ax=ax[0], label=f"$K_h = 10 \ m^2 s^{-1}$",
            clip=(0, 7500), fill=False, color='black')
sns.kdeplot(all_diff["median_time"]/365, ax=ax[1], 
            clip=(0, 6), fill=False, color='black')
sns.kdeplot(all_diff["median_depth"], ax=ax[2], 
            clip=(0, 1000), fill=False, color='black')

#temporal
colors_temp = ["darkred", "orangered", "orange"]
ls_time = [(0, (1, 1)), '--', '-.', (0, (3, 1, 1, 1, 1, 1))]
j = 0
for week in [4, 12, 20]:
    sns.kdeplot(all_temp[week]["counts"], ax=ax[3], label=f"{week} weeks", clip=(0, 7500),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    # kde_temp[week] = _counts.get_lines()[j].get_data() # save the KDE values for entropy computation
    
    sns.kdeplot(all_temp[week]["median_time"]/365, ax=ax[4], label=f"Temporal W{week}", clip=(0, 6),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    sns.kdeplot(all_temp[week]["median_depth"], ax=ax[5], label=f"Temporal W{week}", clip=(0, 1000),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    j += 1

# Add labels 'A', 'B', ... for each subplot in the top left corner
labels = ['A', 'B', 'C', 'D', 'E', 'F']
for i, label in enumerate(labels):
    ax[i].text(0.05, 0.95, label, transform=ax[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

ax[0].set_xlim(0, 6500)
ax[3].set_xlim(0, 6500)
ax[1].set_xlim(0, 6)
ax[4].set_xlim(0, 6)
ax[2].set_xlim(0, 200)
ax[5].set_xlim(0, 200)

ax[0].set_xlabel("Counts")
ax[0].set_ylabel("Density")
ax[0].legend(fontsize=7, shadow=False)

ax[1].set_xlabel("Median Particle Age (years)")
ax[1].set_ylabel("Density")

ax[2].set_xlabel("Median Depth (m)")
ax[2].set_ylabel("Density")

ax[3].set_xlabel("Counts")
ax[3].set_ylabel("Density")
ax[3].legend(fontsize=7, shadow=False)

ax[4].set_xlabel("Median Particle Age (years)")
ax[4].set_ylabel("Density")

ax[5].set_xlabel("Median Depth (m)")
ax[5].set_ylabel("Density")

# set a title
if Latitude_limit is not None:
    plt.suptitle(f"Particles Crossing {Latitude_limit}°N", fontsize=14)
elif Longitude_limit is not None:
    plt.suptitle(f"Particles Crossing {abs(Longitude_limit)}°W", fontsize=14)
plt.tight_layout()
# save the figure
plt.savefig("../figs/FigS5_Connect_tempNspace" + criterium_string + ".png", dpi=300)

# %% Plot the percentage of subpolar trajectories for Mixture temporal and spatial members
fig, ax = plt.subplots(2, 3, figsize=(10, 6))

ax = ax.flatten()
# Plot distributions for percentage of subpolar trajectories

colors_space = ["midnightblue", "blueviolet", "teal"]
ls_space = [(0, (1, 1)), '--', '-.']
j=0

for delta_r in [0.1, 1., 2.]:
    sns.kdeplot(all_mix_space[delta_r]["counts"], ax=ax[0], label=f"Mix. $\delta_r = {delta_r}^o$", clip=(0, 5000),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    # kde_mix_space[delta_r] = _counts.get_lines()[0].get_data() # save the KDE values for entropy computation
    
    sns.kdeplot(all_mix_space[delta_r]["median_time"]/365, ax=ax[1], label=f"Mix. $\delta_r = {delta_r}^o$", clip=(0, 6),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    # sns.kdeplot(all_mix_space[delta_r]["std_time"]/365, ax=ax[5], label=f"Mix. $\delta_r = {delta_r}^o$", clip=(0, 6),
    #             fill=False, color=colors_space[j], linestyle=ls_space[j])
    sns.kdeplot(all_mix_space[delta_r]["median_depth"], ax=ax[2], label=f"Mix. $\delta_r = {delta_r}^o$", clip=(0, 1000),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    j += 1

colors_temp = ["darkred", "orangered", "orange"]
ls_time = [(0, (1, 1)), '--', '-.', (0, (3, 1, 1, 1, 1, 1))]
j = 0
for week in [4, 12, 20]:
    sns.kdeplot(all_mix_temp[week]["counts"], ax=ax[3], label=f"Mix. {week} weeks", clip=(0, 5000),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    # kde_mix_temp[week] = _counts.get_lines()[j].get_data() # save the KDE values for entropy computation
    
    sns.kdeplot(all_mix_temp[week]["median_time"]/365, ax=ax[4], label=f"Mix. {week} weeks", clip=(0, 6),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    # sns.kdeplot(all_mix_temp[week]["std_time"]/365, ax=ax[2], label=f"Mix. {week} weeks", clip=(0, 6),
    #             fill=False, color=colors_temp[j], linestyle=ls_time[j])
    sns.kdeplot(all_mix_temp[week]["median_depth"], ax=ax[5], label=f"Mix. {week} weeks", clip=(0, 1000),
               fill=False, color=colors_temp[j], linestyle=ls_time[j])
    j += 1

# Add labels 'A', 'B', ... for each subplot in the top left corner
labels = ['A', 'B', 'C', 'D', 'E', 'F']
for i, label in enumerate(labels):
    ax[i].text(0.05, 0.95, label, transform=ax[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
    

ax[0].set_xlim(0, 6500)
ax[3].set_xlim(0, 6500)
ax[1].set_xlim(0, 6)
ax[4].set_xlim(0, 6)
ax[2].set_xlim(0, 200)
ax[5].set_xlim(0, 200)

ax[0].set_xlabel("Counts")
ax[0].set_ylabel("Density")
ax[0].legend(fontsize=7, shadow=False)
# ax[0].grid()

ax[1].set_xlabel("Median Particle Age (years)")
ax[1].set_ylabel("Density")
# ax[1].grid()

# ax[2].set_xlabel("STD Drift Time (years)")
ax[2].set_xlabel("Median Depth (m)")
ax[2].set_ylabel("Density")
# ax[2].grid()

ax[3].set_xlabel("Counts")
ax[3].set_ylabel("Density")
ax[3].legend(fontsize=7, shadow=False)
# ax[3].grid()

ax[4].set_xlabel("Median Particle Age (years)")
ax[4].set_ylabel("Density")
# ax[4].grid()

# ax[5].set_xlabel("STD Drift Time (years)")
ax[5].set_xlabel("Median Depth (m)")
ax[5].set_ylabel("Density")
# ax[5].grid()

# set a title
if Latitude_limit is not None:
    plt.suptitle(f"Mixture Particles Crossing {Latitude_limit}°N", fontsize=14)
elif Longitude_limit is not None:
    plt.suptitle(f"Mixture Particles Crossing {abs(Longitude_limit)}°W", fontsize=14)

plt.tight_layout()
plt.savefig("../figs/FigS6_Connect_MIX_tempNspace" + criterium_string + ".png", dpi=300)

###################################################################




#%% ECDFs of the counts and mean and 99% confidence interval for temporal and spatial members
# 
fig, ax = plt.subplots(2, 2, figsize=(8, 6))
ax = ax.flatten()
# Plot CDFs for percentage of subpolar trajectories

colors_space = ["midnightblue", "blueviolet", "teal"]
ls_space = [(0, (1, 1)), '--', '-.']

for j, delta_r in enumerate([0.1, 1., 2.]):
    sns.ecdfplot(all_space[delta_r]["counts"], ax=ax[0], label=f"$\delta_r = {delta_r}^o$", 
                 color=colors_space[j], linestyle=ls_space[j])
    sns.ecdfplot(all_space[delta_r]["median_time"]/365, ax=ax[1], label=f"Spatial dr{delta_r}", 
                 color=colors_space[j], linestyle=ls_space[j])

#diffusion
sns.ecdfplot(all_diff["counts"], ax=ax[0], label=f"$K_h = 10 $",
             color='black', linestyle='-')
sns.ecdfplot(all_diff["median_time"]/365, ax=ax[1], label=f"Spatial dr{delta_r}", 
                 color='black', linestyle='-')


for j, dr_ref in enumerate([0.1, 1., 2.]):
    mean_counts = all_mix_space[dr_ref]["counts"].mean()
    std_counts = 2.576 * all_mix_space[dr_ref]["counts"].std() / np.sqrt(len(all_mix_space[dr_ref]["counts"]))

    ax[0].fill_betweenx([0, 1], mean_counts-std_counts, mean_counts+std_counts,
                        color=colors_space[j], alpha=0.5, label=f"Mix. $\delta_r = {dr_ref}^o$",
                        edgecolor='none')
    ax[0].set_xlim(1000, 5100)
    
    mean_median_time = all_mix_space[dr_ref]["median_time"].mean()/365
    std_median_time = 2.576 * all_mix_space[dr_ref]["median_time"].std()/365 / np.sqrt(len(all_mix_space[dr_ref]["median_time"]))

    ax[1].fill_betweenx([0, 1], mean_median_time-std_median_time, mean_median_time+std_median_time,
                        color=colors_space[j], alpha=0.5, label=f"Mix. $\delta_r = {dr_ref}^o$",
                        edgecolor='none')
    ax[1].set_xlim(2, 5)

colors_temp = ["darkred", "orangered", "orange"]
ls_time = [(0, (1, 1)), '--', '-.', (0, (3, 1, 1, 1, 1, 1))]

for j, week in enumerate([4, 12, 20]):
    sns.ecdfplot(all_temp[week]["counts"], ax=ax[2], label=f"{week} weeks", color=colors_temp[j], linestyle=ls_time[j])
    sns.ecdfplot(all_temp[week]["median_time"]/365, ax=ax[3], label=f"Temporal W{week}", color=colors_temp[j], linestyle=ls_time[j])

for j, dr_ref in enumerate([4, 12, 20]):
    mean_counts = all_mix_temp[dr_ref]["counts"].mean()
    std_counts = 2.576 * all_mix_temp[dr_ref]["counts"].std() / np.sqrt(len(all_mix_temp[dr_ref]["counts"]))

    ax[2].fill_betweenx([0, 1], mean_counts-std_counts, mean_counts+std_counts,
                        color=colors_temp[j], alpha=0.5, label=f"Mix. {dr_ref} weeks", edgecolor='none')
    ax[2].set_xlim(1000, 5100)
    
    mean_median_time = all_mix_temp[dr_ref]["median_time"].mean()/365
    std_median_time = 2.576 * all_mix_temp[dr_ref]["median_time"].std()/365 / np.sqrt(len(all_mix_temp[dr_ref]["median_time"]))

    ax[3].fill_betweenx([0, 1], mean_median_time-std_median_time, mean_median_time+std_median_time,
                        color=colors_temp[j], alpha=0.5, label=f"Mix. {dr_ref} weeks", edgecolor='none')
    ax[3].set_xlim(2, 5)

# Add labels 'A', 'B', 'C', 'D' for each subplot in the top left corner
labels = ['A', 'B', 'C', 'D']
for i, label in enumerate(labels):
    ax[i].text(0.05, 0.95, label, transform=ax[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

ax[0].set_xlabel("Counts")
ax[0].set_ylabel("ECDF")
ax[0].legend(fontsize=7)
ax[0].grid(True, linestyle=':', alpha=0.8)

ax[1].set_xlabel("Median Particle Age (years)")
ax[1].set_ylabel("ECDF")
ax[1].grid(True, linestyle=':', alpha=0.8)

ax[2].set_xlabel("Counts")
ax[2].set_ylabel("ECDF")
ax[2].legend(fontsize=7, loc='lower right')
ax[2].grid(True, linestyle=':', alpha=0.8)

ax[3].set_xlabel("Median Particle Age (years)")
ax[3].set_ylabel("ECDF")
ax[3].grid(True, linestyle=':', alpha=0.8)

plt.tight_layout()
# save the figure

plt.savefig("../figs/Fig3_ECDF_Comparisson" + criterium_string + ".png", dpi=300)


# then compute the entropy of the KDE. No plotting here
    
# Compute the entropy of the KDE
from scipy.stats import entropy

types = ["Spatial 0.1", "Spatial 1.0", "Spatial 2.0", "Temporal 4", "Temporal 12", "Temporal 20"]

entropy_single = np.zeros(6)
entropy_mixture = np.zeros(6)

j = 0
for delta_r in [0.1, 1., 2.]:
    entropy_single[j] = entropy(kde_space[delta_r][1])
    # entropy_mixture[j] = entropy(kde_mix_space[delta_r][1])
    j += 1
    
for j, week in enumerate([4, 12, 20]):
    entropy_single[j+3] = entropy(kde_temp[week][1])
    # entropy_mixture[j+3] = entropy(kde_mix_temp[week][1])
    j += 1

# # Create the dataframe
# entropy_df = pd.DataFrame({"Type": types, "Entropy": [entropy_space[0.1], entropy_space[1.], entropy_space[2.],
#                                                             entropy_temp[4], entropy_temp[12], entropy_temp[20]]})

# entropy_mix_df = pd.DataFrame({"Type": types, "Entropy": [entropy_mix_space[0.1], entropy_mix_space[1.], entropy_mix_space[2.],
#                                                             entropy_mix_temp[4], entropy_mix_temp[12], entropy_mix_temp[20]]})

