''' 
plotters for DUNE-GAMPix 
author: Henry Purcell
'''

import os
import glob
import h5py 
from particle import PDGID
from particle import Particle
from GAMPixTools import edepsim_tools, dune_tools, consts
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import scipy
from scipy import stats

import argparse


################
# plot options #
################
mpl.rcParams['axes.facecolor'] = 'white'
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'serif'
#mpl.rc('text', usetex=True)

mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 2
#mpl.rcParams['xtick.minor.size'] = 4
#mpl.rcParams['xtick.minor.width'] = 2
#mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.top'] = True

mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 2
#mpl.rcParams['ytick.minor.size'] = 4
#mpl.rcParams['ytick.minor.width'] = 2
#mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.right'] = True

mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['axes.grid'] = False

np.set_printoptions(precision=16)


class Plotters:
    def __init__(self,
                 input_edepsim_file, 
                 importance_method = 'primary',
                 single_event = 0,
                 eventIDs = [0,1,2],):
    
        self.input_edepsim_file = input_edepsim_file
        self.eventIDs           = eventIDs
        self.importance_method  = importance_method
        self.single_event       = single_event
        
        self.h5file = h5py.File(input_edepsim_file, 'r')
        self.segments = self.h5file['segments']
        self.trajectories = self.h5file['trajectories']
        
    # make bar plot of energy deposited from each particle type to electrons over whole event
    #method = 'initial' (finds fraction by initial energies of particles) 
    # or 'deposited' (finds fraction by integrating deposited charge from particles)
    #or 'count' (finds fraction by counting number of particles)
    #fraction = True means answer will be given as a fraction of total (unless method = 'count')
    def energyPerEventPlot(self,
                           removePrimaries = True,
                           trueEnergies = True,
                           energy_bounds = None,
                           savefig = False, 
                           title_append = '',
                           show = True):
        
        input_edepsim_file = self.input_edepsim_file
        eventIDs = self.eventIDs
        method = self.importance_method
        
        #initialize ER, extract segments
        ER = dune_tools.EnergyResolution(input_edepsim_file, eventIDs, depth=0)
        eventTrajectories = ER.eventTrajectories
        
        #if removePrimaries, mask out the primary so its energy doesnt dwarf everything else
        if removePrimaries:
            noPrimaryMask = eventTrajectories['trackID'] != 0
            eventTrajectories = eventTrajectories[noPrimaryMask]
            
        #intialize dictionary to fill with energies for different particle types
        pIDs = np.unique(eventTrajectories['pdgId'])
        
        #loop thru pIDs and add importances to dictionary
        pID_dict = {}
        for pID in pIDs:
            pID = int(pID)
            eventDict = ER.particleImportance(method, pID,
                                              energy_bounds = energy_bounds,
                                              removePrimaries = removePrimaries)
            pID_dict[pID] = eventDict

        #plot
        f, ax = plt.subplots()
        bottom = np.zeros_like(eventIDs, dtype = float)
        
        for pID, eventDict in pID_dict.items():
            #check if nucleon, bundle all nucleons together into one color
            if pID >= 1E9:
                kwargs = {'color': 'black'}
            else:
                kwargs = {'label': r'$' + str(Particle.from_pdgid(pID).latex_name) + r'$'}
            events = list(eventDict.keys())
            vals = list(eventDict.values())
            ax.bar(events,
                   vals, 
                   bottom = bottom,
                   **kwargs)
            bottom += vals

        ax.bar(events, np.zeros_like(events), bottom = bottom, color = 'black', label = 'Nuclei') # to add a legend entry for nuclei
        ax.legend(bbox_to_anchor=(1., 1.), fontsize = 9, reverse = True)
        ax.set_xlabel('Event Energy (MeV)')
        ax.set_ylabel(method.capitalize() + ' Energy (MeV)')
        ax.set_title(method.capitalize() + ' Energy in each Particle Type '+title_append)

        #label ticks by energies
        if trueEnergies:
                energyDict = ER.trueEnergy()
                energies = [round(energyDict[event],1) for event in eventIDs]
                ax.set_xticks(eventIDs, energies,rotation = 'vertical') 
        else:
            ax.set_xticks(eventIDs)      
        if method == 'count':
            ax.set_ylabel('Number of Particles')
            ax.set_title( 'Number of Particles of each Type '+ title_append)
        if savefig:
            f.savefig(savefig + '.png', dpi = 300, bbox_inches = 'tight')
        
        if show:
            plt.show()




    # make histogram of particles in energy bins for an event
    def energyHistogramPlot(self,
                            pIDs, 
                            bins = 100, energy_bounds = [0,1000], #MeV
                            Normalize = True, 
                            savefig = False, 
                            title_append = '',
                            show = True):
        
        event = self.single_event
        input_edepsim_file = self.input_edepsim_file
        ER = dune_tools.EnergyResolution(input_edepsim_file, event, depth=0)      
        binBoundaries = np.linspace(energy_bounds[0], energy_bounds[1], bins+1)
            
        #make dictionary for each particle type
        if type(pIDs) == (int):
            pIDs = [pIDs]
        if type(pIDs[0]) == (str):
            pIDs = [int(Particle.from_name(pID).pdgid) for pID in pIDs]
        
        pID_dict = {}
        for pID in pIDs:
            energies = ER.energyList(pID, energy_bounds = None, eventID = event)
            pID_dict[pID] = energies
        
        
        #find particle names 
        pID_names = [r'$'+ Particle.from_pdgid(pID).latex_name +r'$' for pID in list(pID_dict.keys())]
        
        f, ax = plt.subplots()
        
        ax.hist(list(pID_dict.values()),
                label = pID_names,
                bins = binBoundaries,
                histtype='step')
        
        ax.semilogy()
        
        ax.legend(bbox_to_anchor=(1., 1.), fontsize = 9, reverse = True)
        
        ax.set_xlabel('Energy Bins (MeV)')
        if Normalize:
            ax.set_ylabel('Fraction of Particles of Energy')
        else:
            ax.set_ylabel('Number of Particles of Energy')

        ax.set_title('Energy distribution for each particle '+ title_append)
        
    
        if savefig:
            f.savefig(savefig + '.png', dpi = 300, bbox_inches = 'tight')
        
        if show:
            plt.show()


    # make bar plot of energy deposited from each particle type to electrons at each generation for one event
    def integratedEnergyVsGenerationPlot(self,
                                         trueEnergy = True,
                                         savefig = False, 
                                         show = True, 
                                         title_append = ''
                                         ):
        segments = self.segments
        trajectories = self.trajectories
        event = self.single_event
        
        eventMask = segments['eventID'] == event      
        eventSegs  = segments[eventMask]
        
        eventMask = trajectories['eventID'] == event
        eventTrajectories = trajectories[eventMask]
        #intialize dictionary to fill with energies for different particle types
        pIDs = np.unique(eventSegs['pdgId'])
        pIDs = {int(pID) if pID < 1E9 else 'Nuclei'
                for pID in pIDs }  
        
        pID_dict = {pID:{0:0.} for pID in pIDs}
        
        #loop through secondaries, find initial energies and append to corresponding particle type
        for segment in eventSegs:
            
            #find trackID
            trackID = segment['trackID']
            
            ### calculate generation
            generation = 0 
            trajectory = eventTrajectories[trackID]
            parent_id  = trajectory['parentID']
            
            while parent_id != -1:
                trajectory = eventTrajectories[parent_id]
                parent_id  = trajectory['parentID']
                generation += 1
                
            #find particle type and energy in this segment
            pID = segment['pdgId']
            if pID >= 1E9:
                pID = 'Nuclei'
            dE = segment['dE']
            
            #add generation to dictionaries if not in keys
            if generation not in pID_dict[pID].keys():
                for pIDs in pID_dict.keys():
                    pID_dict[pIDs][generation] = 0.
            pID_dict[pID][generation] += dE
        
        #plot 
        f, ax = plt.subplots()
        #bottom like number of generations
        bottom = np.zeros(( len( 
                                list( pID_dict.values() )[0].keys() 
                                )  
                        ))
        
        
        for pID, gen_dict in pID_dict.items():
            #extract values to plot
            generations = list(gen_dict.keys())
            energies = list(gen_dict.values())
            #check if nucleus and generate printable string for legend
            if pID != 'Nuclei':
                pID = r'$'+str(Particle.from_pdgid(pID).latex_name)+r'$'            
           
            #plot bar
            ax.bar(generations,
                energies, 
                label = pID,
                bottom = bottom,
                )
            #increase bottom to add next particle type ontop
            bottom += energies

        ax.legend(bbox_to_anchor=(1., 1.), fontsize = 9, reverse = True)
        ax.set_xlabel('Generation')
        ax.set_ylabel('Energy (MeV)')

        ax.set_xticks(generations)
        
        #print true initial particle energy   
        if trueEnergy:
            eventMask = trajectories['eventID'] == event
            eventTrajectories  = trajectories[eventMask]
            primary = eventTrajectories[eventTrajectories['trackID'] == 0]
            #find energy in this segment
            pID = primary['pdgId']
            m = Particle.from_pdgid(pID).mass #MeV
            if type(m) == type(None): #its a neutrino
                m = 0 #approximate, for GeV neutrinos. This is just to get the energy of the primary 
            p_abs = np.linalg.norm(primary['pxyz_start']) #MeV
            
            #only KE if a baryon or nucleus
            Energy = np.sqrt(m**2 + p_abs**2) - m 
            ax.set_title(f'Deposited Energy at each Generation. Initial Energy: {Energy:.2f} MeV '+ title_append)
        else:
            ax.set_title(f'Deposited Energy at each Generation '+ title_append)
            
        if savefig:
            f.savefig(savefig + '.png', dpi = 300, bbox_inches = 'tight')
        
        if show:
            plt.show()
            

    # make bar plot of summed initial energy of each particle type at each generation for one event
    def initialEnergyVsGenerationPlot(self,
                                      savefig = False,
                                      trueEnergy = True,
                                      removePrimaries = True,
                                      title_append = '',
                                      show = True):
        
        trajectories = self.trajectories
        event = self.single_event
        eventMask = trajectories['eventID'] == event      
        eventTrajectories  = trajectories[eventMask]
        
        #intialize dictionary to fill with energies for different particle types
        pIDs = np.unique(eventTrajectories['pdgId'])
        pIDs = {int(pID) if pID < 1E9 else 'Nuclei'
                for pID in pIDs }  
        pID_dict = {pID:{1:0.} for pID in pIDs}
        
        #loop through secondaries, find initial energies and append to corresponding particle type
        for traj in eventTrajectories:
                    
            ### calculate generation
            generation = 0 
            parent_id  = traj['parentID']
            
            while parent_id != -1:
                trajectory = eventTrajectories[parent_id]
                parent_id  = trajectory['parentID']
                generation += 1
            
            #append to dict (skip primary)
            if removePrimaries:
                if generation == 0: 
                    continue   
                
            #find particle type and energy in this traj
            pID = traj['pdgId']
            
            m = Particle.from_pdgid(pID).mass #MeV
            if type(m) == type(None): #its a neutrino
                m = 0 #approximate, for GeV neutrinos. This is just to get the energy of the primary 
            p_abs = np.linalg.norm(traj['pxyz_start']) #MeV
            
            #only KE if a baryon or nucleus
            if PDGID(pID).is_baryon or PDGID(pID).is_nucleus or pID == 11:
                Energy = np.sqrt(m**2 + p_abs**2) - m 
            else:
                Energy = np.sqrt(m**2 + p_abs**2)
            
            #check if nucleus
            if pID >= 1E9:
                pID = 'Nuclei'
            
            #add generation if not in the dicitonaries yet
            if generation not in pID_dict[pID].keys():
                for pIDs in pID_dict.keys():
                    pID_dict[pIDs][generation] = 0.
            pID_dict[pID][generation] += Energy
    
        #plot
        f, ax = plt.subplots()
        
        #bottom like number of generations
        bottom = np.zeros(( len( 
                                list( pID_dict.values() )[0].keys() 
                                )  
                        ))
        
        
        for pID, gen_dict in pID_dict.items():
            #extract values
            generations = list(gen_dict.keys())
            energies = list(gen_dict.values())
            #check if nucleus, generate printable string
            if pID != 'Nuclei':
                pID = r'$'+str(Particle.from_pdgid(pID).latex_name)+r'$'
            ax.bar(generations,
                energies, 
                label = pID,
                bottom = bottom,
                    )
            bottom += energies

        ax.legend(bbox_to_anchor=(1., 1.), fontsize = 9, reverse = True)
        ax.set_xlabel('Generation')
        ax.set_ylabel('Energy (MeV)')
        ax.set_xticks(generations)
        
        #print true initial particle energy   
        if trueEnergy:
            eventMask = trajectories['eventID'] == event
            eventTrajectories  = trajectories[eventMask]
            primary = eventTrajectories[eventTrajectories['trackID'] == 0]
            #find energy in this segment
            pID = primary['pdgId']
            m = Particle.from_pdgid(pID).mass #MeV
            if type(m) == type(None): #its a neutrino
                m = 0 #approximate, for GeV neutrinos. This is just to get the energy of the primary 
            p_abs = np.linalg.norm(primary['pxyz_start']) #MeV
            
            #only KE if a baryon or nucleus
            Energy = np.sqrt(m**2 + p_abs**2) - m 
            ax.set_title( f'Initial Energy at each Generation. Initial Energy: {Energy:.2f} MeV '+ title_append)
        else:
            ax.set_title( f'Initial Energy at each Generation '+ title_append)
            
        if savefig:
            f.savefig(savefig + '.png', dpi = 300, bbox_inches = 'tight')
        
        if show:
            plt.show()
            
 
 
#plots the reconstructed energy/charge, normalized by the true energy, as a function of importance of some particle

#importance_method = 'initial' (finds fraction by initial energies of particles) 
# or 'deposited' (finds fraction by integrating deposited charge from particles)
# or 'count' (finds fraction by counting number of particles)
# or 'eventID' (finds fraction by event ID)
# or 'primary' (finds fraction by energy of incoming particle)
    
def energy_particle_plot( 
                        input_dict, #dict from energy_vs_particle output or file from dune_tools.main()
                        charge_readout, # = 'Gampix-Pixels' or 'Gampix-Tiles' or 'LArPix' or 'Anode-Grid' or a list of these
                        importance_method, #see above                      
                        recombination = False, #if True, estimated constant recombination (monte carlo truth) is used to correct for recombination. if a number, then this number is used. if an array of the same length as the readouts, then this array is used.
                        divide_by = None, #None, 'primary' or 'deposited'. divides the readout (after factoring out recombination if specified) by either the true primary energy or the deposited charge or energy (depending on if recombination was factored out: youll get either Q_observed/Q_deposited or E_observed/E_deposited)
                        ylims = None,
                        xlims = None,
                        binned_stat = False, hist2d = False,
                        bins = 30, #for binned_stat or hist2d
                        title_append = '', #prepended to title
                        event_type = '', #if False, all events are used. if b"NC" or b"CC", only neutral / charged current events are used
                        savefig = False
                        ):
    
        
    #extract dictionary if input_dict is an h5 file
    if type(input_dict) == str:
        assert(input_dict[-3:] == '.h5'), "input_dict must be a dictionary or an h5 file"
        input_dictionary = {**h5py.File(input_dict, 'r')}
        input_dictionary = {key: input_dictionary[key][()] for key in input_dictionary.keys()} #unpack
        h5py.File(input_dict, 'r').close()
    else: #the dictionary has been passed directly
        input_dictionary = input_dict.copy()
        
    #apply event_type mask to dictionary data to filter out CC or NC events
    if event_type and 'type' in input_dictionary.keys():
        eventMask = np.array(input_dictionary['type']) == event_type
        mask_shape = np.array(input_dictionary['type']).shape
        
        #loop through and find data corresponding to event arrays
        for key in input_dictionary.keys(): 
            if np.array(input_dictionary[key]).shape==mask_shape:
                input_dictionary[key] = input_dictionary[key][eventMask]
        
    ###extract data from dictionary
    #true energy of primaries
    energy_true_list = input_dictionary['energy_true']
    
    #particle importances (x axis data)
    if importance_method == 'primary': #importance catalogued by primary energy
        importance_list = energy_true_list
    else:
        importance_list  = input_dictionary[importance_method]
    
    #readout values (y axis data)
    readout_dict     = {}
    for readout_type in charge_readout:
        readout_dict[readout_type] = input_dictionary[readout_type]    
    
    #recombination
    if recombination == True: #use estimated factors included in the data
        recomb_factors = input_dictionary['recombination']
    elif type(recombination) == float or type(recombination) == int: #use a constant factor
        recomb_factors = np.ones_like(energy_true_list)*recombination
    elif type(recombination) == np.ndarray: #use a factor for each readout
        recomb_factors = recombination
    
    #factor to divide by
    if divide_by == 'primary':
        denominator = energy_true_list #normalize by primary energy
    elif divide_by == 'deposited':
        if recombination:
            denominator = input_dictionary['deposited_energy'][()] #normalize by deposited energy
        else:
            denominator = input_dictionary['deposited_charge'][()] #normalize by deposited charge
    
    #other parameters
    importance_energy_bounds = input_dictionary['importance_energy_bounds'] #for x label
    pID = input_dictionary['pID'] #for x label
    
    if type(pID) == h5py._hl.dataset.Dataset: #extract parameters from h5 dataset
        pID = int(pID[()])
        importance_energy_bounds = importance_energy_bounds[()]
    #find latex string corresponding to pdgID
    latex_particle_str = Particle.from_pdgid(pID).latex_name
    
    #make labels
    #x axis
    xlabel_dict1 = {'initial': 'Summed Initial Energy - Gamma Rays',   
                    'deposited': 'Deposited Energy - Gamma Rays',
                    'count': 'Particle Count - Gamma Rays',
                    'eventID': 'event ID',
                    'primary': 'Primary Energy (MeV)'} 
    if importance_method == 'eventID' or 'primary':
        xlabel = xlabel_dict1[importance_method]
    else:
        xlabel = r'$'+ latex_particle_str + r'$' + ' - ' + xlabel_dict1[importance_method]
    
    #y axis
    ylabel_dict1 = {True:  r'$E_{reco}$', False: r'$Q_{obs}$' }
    ylabel_dict2 = {'primary': {True:r'$ / E_{true}$', False:r'$ / E_{true}$'},
                    'deposited':{True:r'$ / E_{dep}$', False:r'$ / Q_{dep}$'},
                    None: {True:'', False:''}}
    #ylabel_dict3 = {True: '(range = ' + str(importance_energy_bounds)+')', False: ''}     
    ylabel = ylabel_dict1[recombination] \
            + ' ' + ylabel_dict2[divide_by][recombination] #\
     #       + ' ' + ylabel_dict3[importance_energy_bounds != None]
    
    #color dict
    color_dict = {'Gampix-Pixels':"#1A5276", 'Gampix-Tiles': "#02C25F", 'LArPix': "#ED343D", 'Anode-Grid': "#DC12F3"}
    colormap_dict = {'Gampix-Pixels':"Blues", 'Gampix-Tiles': "Greens", 'LArPix': "Reds", 'Anode-Grid': "BuPu"}
    #loop through readout types, factor out recombination and division factor, then plot
    for n,readout_type in enumerate(readout_dict.keys()):
        energy_reco_list = readout_dict[readout_type]
        
        if recombination: #factor out recombination
            recomb_factors[np.where(recomb_factors==0)]=1            
            energy_reco_list = np.array(energy_reco_list)/np.array(recomb_factors)
        
        if divide_by: #divide by factor
            denominator[np.where(denominator==0)]=1
            energy_reco_list = np.array(energy_reco_list)/np.array(denominator)
            
        if binned_stat: #ise scipy binned_statistic to find medians and quantiles for events within a bin over the x axis
            if not 'f' in locals():
                f, ax = plt.subplots(1, 1, squeeze = False, figsize = (6,5))
            
            #define function that defines function for quantiles, taking in a list of values x
            def makeQuantile(q):
                def quantile(x):
                    if len(x) <= 1:
                        return np.nan 
                    else:
                        return np.quantile(x, q) 
                return quantile  
            
            #define median, low quantile and upper quantile using this
            median = makeQuantile(0.5)
            low_quantile = makeQuantile(0.16)
            up_quantile  = makeQuantile(0.84)
            
            #find quantiles
            avgs, edges = stats.binned_statistic(importance_list, energy_reco_list, statistic= median, bins=bins, range=None)[0:2]
            upper_quantile                = stats.binned_statistic(importance_list, energy_reco_list, statistic=up_quantile,  bins=bins, range=None)[0]
            lower_quantile                = stats.binned_statistic(importance_list, energy_reco_list, statistic=low_quantile, bins=bins, range=None)[0]
            
            #plot
            bin_centers = edges[:-1] + 1/2*(edges[1]-edges[0])
            ax[0,0].errorbar(bin_centers, avgs, yerr = (avgs - lower_quantile, upper_quantile - avgs),fmt = 'o', label = readout_type, color = color_dict[readout_type])
            ax[0,0].set_xlabel(xlabel)
            ax[0,0].legend(loc = 'best', reverse = True)

        elif hist2d: #plot a 2d histogram (needs a lot of points)
            if not 'f' in locals():
                f, ax = plt.subplots(1,  len(readout_dict.keys()), squeeze = False,sharey = True, sharex = True, figsize = (len(readout_dict.keys())*5+2,5))
            histo = ax[0,n].hist2d(importance_list, energy_reco_list, bins = bins, label = readout_type,norm='log', cmap = colormap_dict[readout_type])
            ax[0,n].set_title(readout_type)
            ax[0,1].set_xlabel(xlabel)
            #f.colorbar(histo[3])    
        
        else: #scatter plot of all the datapoints
            if not 'f' in locals():
                f,ax = plt.subplots(1,1, squeeze=False, figsize = (6,5))
            ax[0,0].scatter(importance_list, energy_reco_list,s=20, label = readout_type, color = color_dict[readout_type])
            #make labels
            ax[0,0].set_xlabel(xlabel)
            ax[0,0].legend()
    
    if hasattr(input_dictionary, 'pdgMask'): #pdgMask was applied to the data, show this in the y label
        if input_dictionary['pdgMask']:
            pdgMask = input_dictionary['pdgMask']
            latex_str = r'$' + Particle.from_pdgid(pdgMask).latex_name + r'$'
            ax[0,0].set_ylabel(ylabel + ', ' + latex_str + "'s Only")
    else:
        ax[0,0].set_ylabel(ylabel)
        
    #title
    event_type_dict = {b'CC':'Charged-Current ', b'NC':'Neutral-Current ', '':''}
    f.suptitle('' + event_type_dict[event_type] + title_append + ' Events', fontsize = 16)
    
    #ylims
    if ylims and divide_by:
        ax[0,0].set_ylim(ylims)
    #xlims
    if xlims:
        ax[0,0].set_xlim(xlims)
        
    f.tight_layout()
    #save
    if savefig:
        f.savefig(savefig + '.png', dpi = 300, bbox_inches = 'tight')



      
# main 
def main(args):
    savefig_directory = args.output_dir
    if not os.path.exists(savefig_directory):
        os.makedirs(savefig_directory)
            
    f = h5py.File(args.h5file, 'r')
    events = args.multi_events
    single_event = args.single_event
    show = args.show
    which = args.which
    energy_bounds = args.energy_bounds
    importance_method = args.importance_method
    pIDs = args.pids
    charge_readout = args.charge_readout
    remove_primaries = not args.include_primaries
    bins = args.bins
    title_append = args.title
    
    if args.event_type: #plot CC and NC events separately
        event_types = {'Charged-Current':b'CC', 'Neutral-Current':b'NC', 'Both-Types':''}
    else:
        event_types = {'Both_types':''}
        
    if which == 'reco':
        print('starting reco plot')
        input_dict = {**h5py.File(args.h5file, 'r')}
        h5py.File(args.h5file, 'r').close()
        input_dict = {key: input_dict[key][()] for key in input_dict.keys()} #unpack
        binned_stat = args.binned_stat
        divide_by = args.divide_by
        hist2d = args.hist2d
        xlims  = args.xlims
        ylims  = args.ylims
        reconstruction_label_dict = {True: 'E-reco', False: 'Q-deposited'}
    
        for method in np.array(importance_method): 
            for recombination in [True, False]:
                for event_type in event_types.keys():    
                    print('Starting ' + method + ' ' + event_type + ' ' +  reconstruction_label_dict[recombination] + ' plot')
                    energy_particle_plot(input_dict, #dict from energy_vs_particle output or file from dune_tools.main()
                                charge_readout,
                                method,
                                recombination = recombination,
                                divide_by = divide_by,
                                ylims = ylims,
                                xlims = xlims,
                                binned_stat = binned_stat, hist2d = hist2d,
                                bins = bins,
                                event_type = event_types[event_type],
                                title_append = title_append,
                                savefig = savefig_directory + '/' +  method + '_reco_plot_' + reconstruction_label_dict[recombination]+ '_' + event_type)
        print('finished')
        return
    
    for method in np.array(importance_method):  
        if not method in ['initial', 'deposited', 'count']:
            print('importance_methods other than initial, deposited, and count only generate plots for which = reco')
            continue   
        
        PlotterClass = Plotters(args.h5file,
                                importance_method = method,
                                single_event = single_event,
                                eventIDs = events) 

        if (which == 'unique' or which == 'all'):
            print('starting initial energy plot')
            PlotterClass.initialEnergyVsGenerationPlot(
                                removePrimaries = remove_primaries,
                                savefig = savefig_directory + '/initialEnergyVsGenerationPlot',
                                title_append = title_append,
                                show = show)
            print('starting deposited energy plot')
            PlotterClass.integratedEnergyVsGenerationPlot(
                                savefig = savefig_directory + '/integratedEnergyVsGenerationPlot',
                                title_append = title_append,
                                show = show)
            print('starting histogram energy plot')
            PlotterClass.energyHistogramPlot(pIDs,
                                bins = 50, energy_bounds = energy_bounds, #MeV
                                Normalize = True, 
                                title_append = title_append,
                                savefig = savefig_directory + '/energyHistogramPlot', 
                                show = show)
                
        if which == 'multi' or which == 'all':  
            
            PlotterClass.energyPerEventPlot(
                                removePrimaries = remove_primaries,
                                energy_bounds = None,
                                title_append = title_append,
                                savefig = savefig_directory + '/' + method + 'EnergyPlot',
                                show = show)

    print('finished')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    #input
    parser.add_argument('h5file',                type = str,
                        help = "input edep-sim dumped h5 to convert")
    #output
    parser.add_argument('-o', '--output_dir',    type = str,   default = './plots',
                        help = "unique event for generation plots")
    #which plots to generate
    parser.add_argument('-w', '--which',         type = str,   default = 'all',
                        help="which plots to generate. 'all' for all of the particle distribution plots (not reco), 'unique' for unique-event  particle distribution plots, 'multi' for multi event plots")
    #single-event params
    parser.add_argument('-se', '--single_event', type = int,   default = 0,
                        help = " event to use for single-event plots")
    #multi event params
    parser.add_argument('-me', '--multi_events', type = int,   default = [0,1,2,3,4,5], nargs='+',
                        help = "events for multi-event plots, default is all")
    #show
    parser.add_argument('-s', '--show', action = 'store_true', default = False,
                        help = "show plots")
    #energy_bounds
    parser.add_argument('-eb', '--energy_bounds', type = float, default = [0,100],    nargs='+',
                        help = 'energy bounds for energy histogram plot and particle importance (MeV)')
    #pID list
    parser.add_argument('-p', '--pids',           type = int, default = [22, 11, 2112, 2212,],  nargs='+',
                        help = 'particle types for energy histogram plot')
    #remove primaries
    parser.add_argument('-ip', '--include_primaries', action = 'store_true',
                        help = 'include primaries in energy plots. need to call this for neutrino plots where the data doesnt have a primary neutrino in it!')
    ### reco arguments
    #importance method
    parser.add_argument('-im', '--importance_method', type = str, default = ['deposited'], nargs='+',
                        help = "methods for energyPerEventPlot and importance calculations. only deposited, count, and initial work for energyPerEventPlot")
    #charge readout
    parser.add_argument('-cr', '--charge_readout', type = str, default = ['Gampix-Pixels', 'Gampix-Tiles', 'LArPix', 'Anode-Grid'], nargs='+',
                        help = "charge readout to plot in energy_particle_plot, one of Gampix-Pixels, Gampix-Tiles, 'Anode-Grid', or LArPix")
    #divide_by
    parser.add_argument('-db', '--divide_by', type = str, default = None,
                        help = "divide readout by this factor in energy_particle_plot. one of 'primary', 'deposited', or None")
    #binned_stat
    parser.add_argument('-bs', '--binned_stat', action = 'store_true',
                        help = "plot binned statistics in energy_particle_plot")
    #hist2d
    parser.add_argument('-h2d', '--hist2d', action = 'store_true',
                        help = "plot 2d histogram in energy_particle_plot")
    #xlims
    parser.add_argument('-xl', '--xlims', type = float, default = None, nargs='+',
                        help = 'x limits for energy_particle_plot')
    #ylims
    parser.add_argument('-yl', '--ylims', type = float, default = None, nargs='+',
                        help = 'y limits for energy_particle_plot')
    #bins
    parser.add_argument('-b', '--bins', type = int, default = 25, 
                        help = 'number of bins for binned_stat and hist2d')
    #event type
    parser.add_argument('-et', '--event_type', action = 'store_true',
                        help = 'plot event types in energy_particle_plot')
    #title
    parser.add_argument('-t', '--title', type = str, default = '',
                        help = 'append this string to titles (usually the type of primary, eg "Neutrino")')
    args = parser.parse_args()

    main(args)