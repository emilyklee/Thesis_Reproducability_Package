import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import brewer2mpl

if __name__ == '__main__':
    
    species_names_units = ['Oxygen (mmol O$_2$/m$^3$)', 
                            'Phosphate (mmol P/m$^3$)',
                            'Nitrate (mmol N/m$^3$)', 
                            'Ammonium (mmol N/m$^3$)',
                            'Nitrogen sink (mmol N/m$^3$)', 
                            'Silicate (mmol Si/m$^3$)',
                            'Reduction Equivalents (mmol S/m$^3$)',
                            'Pelagic Bacteria (mg C/m$^3$)',
                            'Pelagic Bacteria (mmol N/m$^3$)',
                            'Pelagic Bacteria (mmol P/m$^3$)',
                            'Diatoms (mg C/m$^3$)',
                            'Diatoms (mmol N/m$^3$)',
                            'Diatoms (mmol P/m$^3$)',
                            'Diatoms (mg Chl-a/m$^3$)',
                            'Diatoms (mmol Si/m$^3$)',
                            'Flagellates (mg C/m$^3$)',
                            'Flagellates (mmol N/m$^3$)',
                            'Flagellates (mmol P/m$^3$)',
                            'Flagellates (mg Chl-a/m$^3$)',
                            'PicoPhytoplankton (mg C/m$^3$)',
                            'PicoPhytoplankton (mmol N/m$^3$)',
                            'PicoPhytoplankton (mmol P/m$^3$)',
                            'PicoPhytoplankton (mg Chl-a/m$^3$)',
                            'Large Phytoplankton (mg C/m$^3$)',
                            'Large Phytoplankton (mmol N/m$^3$)',
                            'Large Phytoplankton (mmol P/m$^3$)',
                            'Large Phytoplankton (mg Chl-a/m$^3$)',
                            'Carnivorous Mesozooplankton  (mg C/m$^3$)',
                            'Carnivorous Mesozooplankton  (mmol N/m$^3$)',
                            'Carnivorous Mesozooplankton  (mmol P/m$^3$)',
                            'Omnivorous Mesozooplankton (mg C/m$^3$)',
                            'Omnivorous Mesozooplankton (mmol N/m$^3$)',
                            'Omnivorous Mesozooplankton (mmol P/m$^3$)',
                            'Microzooplankton (mg C/m$^3$)',
                            'Microzooplankton (mmol N/m$^3$)',
                            'Microzooplankton (mmol P/m$^3$)',
                            'Heterotrophic Nanoflagellates (mg C/m$^3$)',
                            'Heterotrophic Nanoflagellates (mmol N/m$^3$)',
                            'Heterotrophic Nanoflagellates (mmol P/m$^3$)',
                            'Labile Dissolved Organic Carbon (mg C/m$^3$)',
                            'Labile Dissolved Organic Nitrogen (mmol N/m$^3$)',
                            'Labile Dissolved Organic Phosphate (mmol P/m$^3$)',
                            'Semi-labile Dissolved Organic Carbon (mg C/m$^3$)',
                            'Semi-refractory Dissolved Organic Carbon (mg C/m$^3$)',
                            'Particulate Organic Carbon (mg C/m$^3$)',
                            'Particulate Organic Nitrate (mmol N/m$^3$)',
                            'Particulate Organic Phosphate (mmol P/m$^3$)',
                            'Particulate Organic Silicate (mmol Si/m$^3$)',
                            'Dissolved Inorganic Carbon (mg C/m$^3$)',
                            'Total Alkalinity (mmol Eq/m$^3$)'
                            ]
    
    # Set style for plots
    plt.rc('font', family='serif', size=22)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rc('axes', labelsize=20, linewidth=2)
    # legend defaults
    plt.rc('legend', framealpha=1.0, facecolor='white', frameon=True, edgecolor='white')
    # Plotting colors
    bmap = brewer2mpl.get_map('Paired', 'qualitative', 10)
    colors = bmap.mpl_colors
    
    # Read in data for plots
    full_50_t = numpy.genfromtxt('BFM_Reduced_Models_Data/full_50_t_data.csv')
    full_50_c = numpy.genfromtxt('BFM_Reduced_Models_Data/full_50_c_data.csv', delimiter=' ')
    reduced_44_t = numpy.genfromtxt('BFM_Reduced_Models_Data/reduced_44_t_data.csv')
    reduced_44_c = numpy.genfromtxt('BFM_Reduced_Models_Data/reduced_44_c_data.csv', delimiter=' ')
    reduced_21_t = numpy.genfromtxt('BFM_Reduced_Models_Data/reduced_21_t_data.csv')
    reduced_21_c = numpy.genfromtxt('BFM_Reduced_Models_Data/reduced_21_c_data.csv', delimiter=' ')
    reduced_41_t = numpy.genfromtxt('BFM_Reduced_Models_Data/reduced_41_t_data.csv')
    reduced_41_c = numpy.genfromtxt('BFM_Reduced_Models_Data/reduced_41_c_data.csv', delimiter=' ')
    reduced_4_t = numpy.genfromtxt('BFM_Reduced_Models_Data/reduced_4_t_data.csv')
    reduced_4_c = numpy.genfromtxt('BFM_Reduced_Models_Data/reduced_4_c_data.csv', delimiter=' ')
    
    # ------------------------------------------------------------------------
    # Reduced 21 species model plot
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.grid(linestyle='--')
    ax.plot(full_50_t/86400/365, full_50_c[48], 'o', markersize=2, color = colors[7], label='Full model')
    ax.plot(reduced_21_t/86400/365, reduced_21_c[48], 's', markersize=2, color = colors[1], label='Reduced model')
    ax.set_ylabel(species_names_units[48])
    ax.set_xlabel('Time (year)')
    ax.set_xlim(5,10)
    ax.set_ylim(bottom = 27000)
    ax.yaxis.set_ticks(numpy.arange(27000, 29500, 500))
    ax.legend(loc='center', bbox_to_anchor=(0.5,1.05), frameon=False, ncol=2, markerscale=5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
    plt.savefig('Plots/reduced_21_DIC.pdf', bbox_inches = "tight")
    plt.close('all')
    
    # ------------------------------------------------------------------------
    # Reduced 41 species model plot
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.grid(linestyle='--')
    ax.plot(full_50_t/86400/365, full_50_c[45], 'o', markersize=2, color = colors[7], label='Full model')
    ax.plot(reduced_41_t/86400/365, reduced_41_c[45], 's', markersize=2, color = colors[1], label='Reduced model')
    ax.set_ylabel(species_names_units[45])
    ax.set_xlabel('Time (year)')
    ax.set_xlim(5,10)
    ax.set_ylim(bottom=4.6, top = 5.6)
    ax.legend(loc='center', bbox_to_anchor=(0.5,1.05), frameon=False, ncol=2, markerscale=5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
    plt.savefig('Plots/reduced_41_PON.pdf', bbox_inches = "tight")
    plt.close('all')
    
    # ------------------------------------------------------------------------
    # Reduced 4 species model plot
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.grid(linestyle='--')
    ax.plot(full_50_t/86400/365, full_50_c[0], 'o', markersize=2, color = colors[7], label='Full model')
    ax.plot(reduced_4_t/86400/365, reduced_4_c[0], 's', markersize=2, color = colors[1], label='Reduced model')
    ax.set_ylabel(species_names_units[0])
    ax.set_xlabel('Time (year)')
    ax.set_xlim(5,10)
    ax.set_ylim(top=300)
    ax.legend(loc='center', bbox_to_anchor=(0.5,1.05), frameon=False, ncol=2, markerscale=5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
    plt.savefig('Plots/reduced_4_oxygen.pdf', bbox_inches = "tight")
    plt.close('all')
    
    # ------------------------------------------------------------------------
    # Reduced 44 species model plots
    # increase font size
    plt.rc('font', family='serif', size=30)
    plt.rc('xtick', labelsize=28)
    plt.rc('ytick', labelsize=28)
    plt.rc('axes', labelsize=28, linewidth=2)
    
    # Concentration data for plots
    phyto_c_sum_full = full_50_c[10] + full_50_c[15] + full_50_c[19] + full_50_c[23]
    phyto_c_sum_reduced = reduced_44_c[10] + reduced_44_c[15] + reduced_44_c[19] + reduced_44_c[23]
    phyto_chl_sum_full = full_50_c[13] + full_50_c[18] + full_50_c[22] + full_50_c[26]
    phyto_chl_sum_reduced = reduced_44_c[13] + reduced_44_c[18] + reduced_44_c[22] + reduced_44_c[26]
    bact_zoo_c_sum_full = full_50_c[7] + full_50_c[27] + full_50_c[30] + full_50_c[33] + full_50_c[36]
    bact_zoo_c_sum_reduced = reduced_44_c[7] + reduced_44_c[27] + reduced_44_c[30] + reduced_44_c[33] + reduced_44_c[36]
    
    fig, axes = plt.subplots(3, 1, figsize=(20,15), sharex=True, gridspec_kw={'hspace': 0.1}
                             )
    # Phyto chlorophyll
    ax = axes[0]
    ax.plot(full_50_t/86400/365, phyto_chl_sum_full, 'o', markersize=1, color = colors[7], label='Full model')
    ax.plot(reduced_44_t/86400/365, phyto_chl_sum_reduced, 's', markersize=1, color = colors[1], label='Reduced model')
    ax.set_ylabel('(mg Chl-a/m$^3$)')
    ax.legend(loc='upper right', bbox_to_anchor=(1,1.1), markerscale=12)#, ncol=2, bbox_to_anchor=(0.5,1.05)
    ax.set_xlim(5,10)
    ax.set_ylim(top=1)
    
    # Phyto carbon
    ax = axes[1]
    ax.plot(full_50_t/86400/365, phyto_c_sum_full, 'o', markersize=1, color = colors[7], label='Full model')
    ax.plot(reduced_44_t/86400/365, phyto_c_sum_reduced, 's', markersize=1, color = colors[1], label='Reduced model')
    ax.set_ylabel('(mg C/m$^3$)')
    ax.set_xlim(5,10)
    ax.set_ylim(top=40)
    
    # Non-photosynthesizers (bact + zoo) carbon
    ax = axes[2]
    ax.plot(full_50_t/86400/365, bact_zoo_c_sum_full, 'o', markersize=1, color = colors[7], label='Full model')
    ax.plot(reduced_44_t/86400/365, bact_zoo_c_sum_reduced, 's', markersize=1, color = colors[1], label='Reduced model')
    ax.set_ylabel('(mg C/m$^3$)')
    ax.set_xlabel('Time (year)')
    ax.set_xlim(5,10)
    ax.set_ylim(top=30)
    
    # Put lable on subplots
    axes[0].text(5.1, 0.9, 'Phytoplankton chlorophyll', horizontalalignment='left', verticalalignment='center')
    axes[1].text(5.1, 35, 'Phytoplankton carbon', horizontalalignment='left', verticalalignment='center')
    axes[2].text(5.1, 27, 'Non-photosynthesizers carbon', horizontalalignment='left', verticalalignment='center')
    
    for ax in axes.reshape(-1):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
        ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
        ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
        ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
        
        ax.grid(axis='both', color="0.9", linestyle='-', linewidth=1)
    
    plt.savefig('Plots/reduced_44.pdf', bbox_inches = "tight")
    