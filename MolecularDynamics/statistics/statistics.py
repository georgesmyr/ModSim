import pandas as pd
import matplotlib.pyplot as plt

def plot_summary(filename, tit):

    statistics = pd.read_csv(filename).set_index('time')
    means = statistics[statistics.index > 20].mean()

    fig, axes= plt.subplots(nrows=5, figsize=(10, 15))
    plt.suptitle(tit)

    axes[0].plot(statistics.index, statistics.potential_energy_pp, label = 'Potential energy per particle')
    ax.axhline(y=means.potential_energy_pp , color='r')
    axes[0].set_xlabel('Time')
    axes[0].set_ylabel('Potential energy per particle')

    axes[1].plot(statistics.index, statistics.total_energy_pp, label = 'Total energy per particle')
    axes[1].set_xlabel('Time')
    axes[1].set_ylabel('Total energy per particle')

    axes[2].plot(statistics.index, statistics.temperature, label = 'Temperature')
    axes[2].set_xlabel('Time')
    axes[2].set_ylabel('Temperature')

    axes[3].plot(statistics.index, statistics.velo_autocor, label = 'Velocity autocorrelation')
    axes[3].set_xlabel('Time')
    axes[3].set_ylabel('Velocity autocorrelation')

    axes[4].plot(statistics.index, statistics.diffusion_coef, label = 'Diffusion coefficient')
    axes[4].set_xlabel('Time')
    axes[4].set_ylabel('Diffusion coefficient')

    plt.savefig(filename[:-4] + "_summary.png", dpi = 300)


    titles = {'statistics_01_1_01_NVE.csv' : 'NVE | Density = 0.1 | Temperature = 1',
          'statistics_01_1_01_NVT.csv' : 'NVT | Density = 0.1 | Temperature = 1',
          'statistics_01_025_01_NVE.csv' : 'NVE | Density = 0.1 | Temperature = 0.25',
          'statistics_01_025_01_NVT.csv' : 'NVT | Density = 0.1 | Temperature = 0.25',
          'statistics_12_1_01_NVE.csv' : 'NVE | Density = 1.2 | Temperature = 1',
          'statistics_12_1_01_NVT.csv' : 'NVT | Density = 1.2 | Temperature = 1',
          'statistics_0025_125_01_NVE.csv' : 'NVE | Density = 0.025 | Temperature = 1.25',
          'statistics_0025_125_01_NVT.csv' : 'NVT | Density = 0.025 | Temperature = 1.25',
          'statistics_085_1_01_NVE.csv' : 'NVE | Density = 0.85 | Temperature = 1',
          'statistics_085_1_01_NVT.csv' : 'NVT | Density = 0.85 | Temperature = 1'}

for file, title in titles.items():
    plot_summary(file, title)
