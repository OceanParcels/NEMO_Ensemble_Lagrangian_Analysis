# %%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

#%%
K_h = 1000
location = 'Cape_Hatteras'
member = 48  # member
path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}_15000.nc"

P_m = xr.open_dataset(path)

K_h = 10
location = 'Cape_Hatteras'
member = 48  # member
path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}.nc"

P_m_10 = xr.open_dataset(path)

#%%
fig, ax1 = plt.subplots()

color1 = 'tab:blue'
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Entropy (bits)', color=color1)
l1, = ax1.plot(P_m['time'], P_m['entropy'], label=r'$K_h=1000 \ m^2s^{-1}$', color=color1)
l2, = ax1.plot(P_m_10['time'], P_m_10['entropy'], label=r'$K_h=10 \ m^2s^{-1}$', color=color1, ls='--')
ax1.tick_params(axis='y', labelcolor=color1)
ax1.set_xlim(0, 2189)
# ax1.semilogx()

ax2 = ax1.twinx()
color2 = 'tab:orange'
ax2.set_ylabel('Number of particles binned', color=color2)
ax2.plot(P_m['time'], P_m['number_particles_binned'], label='Number of particles binned', color=color2)
ax2.plot(P_m_10['time'], P_m_10['number_particles_binned'], label='Number of particles binned K_h=10', color=color2, ls='--')
ax2.tick_params(axis='y', labelcolor=color2)

ax3 = ax1.twinx()
color3 = 'tab:green'
ax3.spines['right'].set_position(('outward', 60))
ax3.set_ylabel('Entropy (bits)/ Number of particles', color=color3)
ax3.plot(P_m['time'], P_m['entropy']/P_m['number_particles_binned'], label='Entropy/Number of particles', color=color3)
ax3.plot(P_m_10['time'], P_m_10['entropy']/P_m_10['number_particles_binned'], label='Entropy/Number of particles K_h=10', color=color3, ls='--')
ax3.tick_params(axis='y', labelcolor=color3)


# # Combine legends
lines = [l1, l2]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc='lower right')

plt.tight_layout()

plt.savefig('../figs/FigS10-diff1000vsdiff10Entropy_curves_15000.png', dpi=300)


# plt.ylabel('Value')

# %%
K_h = 1000
location = 'Cape_Hatteras'
member = 48  # member
path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}_subsample_7500.nc"

P_m_sub = xr.open_dataset(path)

# K_h = 10
# location = 'Cape_Hatteras'
# member = 48  # member
# path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}.nc"

# P_m_10 = xr.open_dataset(path)

#%%
fig, ax1 = plt.subplots()

color1 = 'tab:blue'
color2 = 'tab:orange'
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Entropy (bits)')
l1, = ax1.plot(P_m['time'], P_m['entropy'], label=r'$K_h=1000 \ m^2s^{-1}$')
l2, = ax1.plot(P_m_10['time'], P_m_10['entropy'], label=r'$K_h=10 \ m^2s^{-1}$', ls='--')
ax1.plot(P_m_sub['time'], P_m_sub['entropy'], label=r'$K_h=1000 \ m^2s^{-1}$, subsampled 7500 particles')
# ax1.plot(P_m_10['time'], P_m_10['entropy'], label=r'$K_h=10 \ m^2s^{-1}$, subsampled 4000 particles', ls='--')


# ax1.tick_params(axis='y', labelcolor=color1)
ax1.set_xlim(0, 2189)
# ax1.semilogx()

# ax2 = ax1.twinx()
# color2 = 'tab:orange'
# ax2.set_ylabel('Number of particles binned', color=color2)
# ax2.plot(P_m['time'], P_m['number_particles_binned'], label='Number of particles binned', color=color2)
# ax2.plot(P_m_10['time'], P_m_10['number_particles_binned'], label='Number of particles binned K_h=10', color=color2, ls='--')
# ax2.tick_params(axis='y', labelcolor=color2)

# ax3 = ax1.twinx()
# color3 = 'tab:green'
# ax3.spines['right'].set_position(('outward', 60))
# ax3.set_ylabel('Entropy (bits)/ Number of particles', color=color3)
# ax3.plot(P_m['time'], P_m['entropy']/P_m['number_particles_binned'], label='Entropy/Number of particles', color=color3)
# ax3.plot(P_m_10['time'], P_m_10['entropy']/P_m_10['number_particles_binned'], label='Entropy/Number of particles K_h=10', color=color3, ls='--')
# ax3.tick_params(axis='y', labelcolor=color3)


# # Combine legends
lines = [l1, l2]
labels = [line.get_label() for line in lines]
ax1.legend(loc='lower right')

plt.savefig('../figs/FigS10-diff1000vsdiff10Entropy_curves_subsample_7500.png', dpi=300)
# %%
