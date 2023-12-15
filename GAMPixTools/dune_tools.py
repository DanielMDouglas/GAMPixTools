import numpy as np
import h5py
import pickle
import os
from tqdm import tqdm
import math
import time
from math import log, isnan
from particle import Particle, PDGID, ParticleNotFound
from GAMPixTools import edepsim_tools, consts, electron_track_tools

#calculate total integrated charge readout 
class EnergyResolution:
    def __init__(self,input_edepsim_file, eventIDs, depth=0, allEvents = False):
        
        self.input_edepsim_file = input_edepsim_file
        self.eventIDs = eventIDs
        if type(self.eventIDs) == int or type(self.eventIDs) == float:
            self.eventIDs = [self.eventIDs]
        self.depth = depth
    
        self.h5file = h5py.File(self.input_edepsim_file, 'r')
        
        self.segments     = self.h5file['segments']
        self.trajectories = self.h5file['trajectories']
        
        if not allEvents:
            self.eventTrajMask = np.sum([self.trajectories['eventID'] == event for event in self.eventIDs], axis = 0).astype(bool)             
            self.eventSegMask  = np.sum([self.segments['eventID']     == event for event in self.eventIDs], axis = 0).astype(bool)            
            
            self.eventTrajectories = self.trajectories[self.eventTrajMask]
            self.eventSegments     = self.segments[self.eventSegMask]
        else:
            self.eventTrajectories = self.trajectories
            self.eventSegments     = self.segments
            
        if 'primaries' in self.h5file.keys():
            primaries = self.h5file['primaries']
            if not allEvents:
                self.eventPrimaryMask = np.sum([primaries['eventID'] == event for event in self.eventIDs], axis = 0).astype(bool)
                self.eventPrimaries = primaries[self.eventPrimaryMask]
            else:
                self.eventPrimaries = primaries
                
#returns dictionary primary energies for each event
    def trueEnergy(self, dictionary = True):
        
        if dictionary:
            eventsEnergy = {}
        else: 
            eventsEnergy = []
        
        if hasattr(self, 'eventPrimaries'):
            primaries = self.eventPrimaries
        else:    
            primaryMask = self.eventTrajectories['parentID'] == -1
            primaries =  self.eventTrajectories[primaryMask]
            
        for eventID in self.eventIDs:
            
            #find primary for this event
            eventMask = primaries['eventID'] == eventID
            eventPrimaries = primaries[eventMask]
            summedPrimaryEnergy = 0 #usually there will be one primary unless (e.g. Genie data) the primary isnt shown in which case this sums the energy of the secondaries to get the primary energy
            for primary in eventPrimaries:
                if 'E' in primary.dtype.names:
                    summedPrimaryEnergy += primary['E']*1000 #MeV
                    continue
                #find energy in this segment
                pID = primary['pdgId']
                m = Particle.from_pdgid(pID).mass #MeV
                if type(m) == type(None): #its a neutrino
                    m = 0 #approximate, for GeV neutrinos. This is just to get the energy of the primary
                p_abs = np.linalg.norm(primary['pxyz_start']) #MeV
                
                #only KE if a baryon or nucleus
                Energy = np.sqrt(m**2 + p_abs**2) - m 
                summedPrimaryEnergy += Energy
                
            if dictionary:    
                eventsEnergy[eventID] = summedPrimaryEnergy
            else:
                eventsEnergy.append(summedPrimaryEnergy)
        
        return eventsEnergy
    
    
#returns dictionary of total deposited charge or energy for each event
    def trueDeposited(self, 
                      dictionary = True):
        
        #initialize dictionary/list to fill with values per event
        if dictionary:
            eventsCharge = {}
            eventsEnergy = {}
        else: 
            eventsCharge = []
            eventsEnergy = []
            
        segs = self.eventSegments
        
        #loop through events and find the deposited energy/charge for each  
        for eventID in self.eventIDs:
            
            #find primary for this event
            eventMask = segs['eventID'] == eventID
            eventIDsegs = segs[eventMask]
                        
            #quench to generate n_electrons counts at each segment and sum over these to find total number of deposited electrons
            edepsim_tools.quench(eventIDsegs, consts.BIRKS)
            totalCharge = np.sum(eventIDsegs['n_electrons'])
            
            # sum over deposited energies of each segment
            totalEnergy = 0
            for seg in eventIDsegs:
                if 'dE_ionizing' in seg.dtype.names:
                    dE = seg["dE_ionizing"]
                else:
                    dE = seg['dE']
                totalEnergy += dE
        
           #append to dict or list
            if dictionary:    
                eventsCharge[eventID] = totalCharge
                eventsEnergy[eventID] = totalEnergy
            else:
                eventsCharge.append(totalCharge)
                eventsEnergy.append(totalEnergy)
        
        #output
        return {'deposited_energy':eventsEnergy, 'deposited_charge':eventsCharge}
    
    
    #finds average recombination factor, for reconstruction
    #recombination_method = 'Birks' or 'Box'
    def averageRecombination(self, eventID, recombination_method = consts.BIRKS):
        #find relevant trajectories
        eventMask = self.segments['eventID'] == eventID
        segments  = self.segments[eventMask]
        
        #quench for n_electrons (applies recombination)
        edepsim_tools.quench(segments, recombination_method) 

        #find total deposited dE
        if 'dE_ionizing' in segments.dtype.names:
            dE = segments["dE_ionizing"]
        else:
            dE = segments['dE']
        
        #find deposited charge
        dQ = segments['n_electrons']
        
        if np.sum(dE) == 0 or np.sum(dQ) == 0: # no deposited charge. return None which will lead to this event being skipped
            return 1/ consts.W_ION
        
        E = np.sum(dE)
        Q = np.sum(dQ)
        
        QoverE = Q / E
        
        #average recombination factor
        return QoverE
    
    
    #returns dictionary of average recombination for each event
    def eventRecombination(self, recombination_method = consts.BIRKS, dictionary = True):
        
        if dictionary:
            Recombs = {}
        else: 
            Recombs = []
        
        #loop thru events, find recombination factor using above function
        for eventID in self.eventIDs:
            
            recomb = self.averageRecombination(eventID, recombination_method = recombination_method)
                
            if dictionary:    
                Recombs[eventID] = recomb
            else:
                Recombs.append(recomb)
        
        return Recombs
        
        
    #finds total charge or energy at the pixels and Gampix-Tiles
    #returns libraries with events as keys
    def reconstruct(self, 
                    charge_readout = 'Gampix-Pixels',
                    sum = True,
                    remove_attenuation = True,
                    pdgMask = None,
                    simulation_bounds = [[-29000,29000], [-7250,7250], [-12000,0]],
                    randomize = False,
                    recombination_method = consts.BIRKS,
                    dictionary = True,
                    starttime = 0.):
        
        #make charge_readout a list
        if type(charge_readout) == str:
            charge_readout = [charge_readout]
       
        #initialize track objects per event
        print('initializing track objects')
        tracks = {eventID:electron_track_tools.Track(self.input_edepsim_file,
                                    eventID,
                                    input_format = 'dumpTree',
                                    pdgMask = pdgMask,
                                    simulation_bounds = simulation_bounds,
                                    randomize = randomize)
                  for eventID in self.eventIDs}
        
        #initialize dictionary of readout types that will be outputed
        readout_dict = {}
          
        #loop through readout types
        for readout_type in charge_readout:
            print('reconstructing ' + readout_type + ' readout, t = ' + str(time.time() - starttime) + ' s')
            #initialize charge dictionary to fill with values per event
            if dictionary:
                charge = {}
            else:
                charge = []
            #loop thru events and events and find triggered pixel charges
            for eventID in self.eventIDs:
                print('event ' + str(eventID) + ', t = ' + str(time.time() - starttime) + ' s')
                if readout_type == 'LArPix':
                    tracks[eventID].reset_params(charge_readout='LArPix')
                elif readout_type == 'Anode-Grid':
                    tracks[eventID].reset_params(charge_readout='AnodeGridD')
                else:
                    tracks[eventID].reset_params(charge_readout='GAMPixD')
                
                tracks[eventID].readout_charge(self.depth)

                if readout_type   == 'Gampix-Pixels' or readout_type == 'LArPix':
                    charge_instance = tracks[eventID].pixel_samples["samples_triggered"]
                elif readout_type == 'Gampix-Tiles':
                    charge_instance  = tracks[eventID].coarse_tiles_samples["samples_triggered"]
                elif readout_type == 'Anode-Grid':
                    charge_instance  = tracks[eventID].anode_grid_samples["samples_triggered"]
                else:
                    raise ValueError('charge_readout must be "Gampix-Pixels", "Gampix-Tiles", "LAriPix", or "Anode-Grid"')
                
                #remove attenuation
                 #TODO: implement for Anode-Grid. currently no r values with which to calculate attenuation from anode grid.
                        #locations  = tracks[eventID].anode_grid_samples["r_triggered"]
                if remove_attenuation and readout_type != 'Anode-Grid':
                    #find locations of hits
                    if readout_type  == 'Gampix-Pixels' or readout_type == 'LArPix':
                        locations = tracks[eventID].pixel_samples["r_triggered"]
                    elif readout_type == 'Gampix-Tiles':
                        locations  = tracks[eventID].coarse_tiles_samples["r_triggered"]
                    
                    #find depth of hits
                    depth = - (locations[2] - self.depth)
                    
                    #coefficient for probability
                    drift_length_constant = tracks[eventID].params.charge_drift['drift_length']
                    
                    #calculate survival probs
                    survival_prob = np.exp(-depth / drift_length_constant)
                    
                    #divide measured charge by survival probs to reconstruct true values
                    charge_instance /= survival_prob

                #sum over pixels and tiles
                if sum:
                    charge_instance = charge_instance.sum()
                                
                #append to dict or list
                if dictionary:
                    charge[eventID] = charge_instance
                else:
                    charge.append(charge_instance)
            
            #append to readout_dict
            readout_dict[readout_type] = charge
            
        return readout_dict
    
        
   
    
    #finds fraction of energy of event corresponding to particle with pdgID
    #method = 'initial' (finds fraction by initial energies of particles) 
    # or 'deposited' (finds fraction by integrating deposited charge from particles)
    #or 'count' (finds fraction by counting number of particles)
    #fraction = True means answer will be given as a fraction of total (unless method = 'count')
    #energy_bounds = [min,max] in MeV
    #returns dictionary catalogued by event
    def particleImportance(self, method, pdgID,
                           energy_bounds = None,
                           removePrimaries = True,): 
        
        #check input type
        assert(type(pdgID) == str or type(pdgID) == int), 'pdgID must be a string or int'
        try:
            pdgID = int(Particle.from_name(pdgID).pdgid) 
        except ParticleNotFound:
            pdgID = int(pdgID)
        
        #particle mass
        m = Particle.from_pdgid(pdgID).mass #MeV
        if type(m) == type(None): #its a neutrino
            m = 0 #approximate, for GeV neutrinos. This is just to get the energy of the primary

        if method == 'count':
     
            #loop through trajectories and find the total count of chosen particle
            particleCount      = {eventID:0 for eventID in self.eventIDs}
            
            pdgTrajectories = self.eventTrajectories[self.eventTrajectories['pdgId'] == pdgID]
            
            for traj in pdgTrajectories:
                
                #eventID
                eventID = traj['eventID']
                
                if energy_bounds:
                    #find energy in this segment  
                    p_abs = np.linalg.norm(traj['pxyz_start']) #MeV
                    
                    #include mass energy if the particle came from pair production (currently a bad fix, since this could include compton scattering)
                    if PDGID(pdgID).is_lepton:
                        parentID = traj['parentID']
                        parentTraj = self.eventTrajectores
                        if parentTraj['pdgId'] == 22:
                            Energy = np.sqrt(m**2 + p_abs**2) 
                        else:
                            Energy = np.sqrt(m**2 + p_abs**2) - m
                    else:
                        Energy = np.sqrt(m**2 + p_abs**2) - m
                    
                    if Energy < energy_bounds[0] or Energy > energy_bounds[1]:
                        continue
                particleCount[eventID] += 1
                   
            #output
            if len(self.eventIDs) == 1:
                return particleCount[self.eventIDs[0]]
            else:
                return particleCount
         
        #add up initial energy method   
        elif method == 'initial':       

            #loop through trajectories and find the total sum of initial energies and the sum of initial energies of chosen particle
            pdgInitialEnergy = {eventID:0. for eventID in self.eventIDs} 
            pdgTrajectories = self.eventTrajectories[self.eventTrajectories['pdgId'] == pdgID]
            
            #remove primaries
            if removePrimaries:
                noPrimaryMask = pdgTrajectories['parentID'] != -1
                pdgTrajectories = pdgTrajectories[noPrimaryMask]
        
            for traj in pdgTrajectories:
                
                #eventID
                eventID = traj['eventID']
                
                #find energy in this segment 
                p_abs = np.linalg.norm(traj['pxyz_start']) #MeV
                
                #only KE if a baryon or nucleus
                if PDGID(pdgID).is_baryon or PDGID(pdgID).is_nucleus:
                    Energy = np.sqrt(m**2 + p_abs**2) - m 
                else:
                    Energy = np.sqrt(m**2 + p_abs**2) - m
                
                #skip if not within bounds
                if energy_bounds:
                    if Energy < energy_bounds[0] or Energy > energy_bounds[1]:
                        continue
                    
                #add to pdgInitialEnergy
                pdgInitialEnergy[eventID] += Energy
            
            #output
            if len(self.eventIDs) == 1:
                return pdgInitialEnergy[self.eventIDs[0]]
            else:
                return pdgInitialEnergy
    
    
        
        elif method == 'deposited':
            
            #loop through trajectories and find the total sum of deposited dEs and the sum of initial energies of chosen particle
            pdgDepositedEnergy   = {eventID:0. for eventID in self.eventIDs}

            pdgSegments = self.eventSegments[self.eventSegments['pdgId'] == pdgID]
            if removePrimaries:
                noPrimaryMask = pdgSegments['trackID'] != 0
                pdgSegments = pdgSegments[noPrimaryMask]

            for seg in pdgSegments:
                
                #eventID
                eventID = seg['eventID']
                #find energy deposited in this segment
                if 'dE_ionizing' in seg.dtype.names:
                    dE = seg["dE_ionizing"]
                else:
                    dE = seg['dE']
            
                #append to pdgDepositedEnergy 
                pdgDepositedEnergy[eventID] += dE
                
            #return
            if len(self.eventIDs) == 1:
                return pdgDepositedEnergy[self.eventIDs[0]]
            else:
                return pdgDepositedEnergy
        
        
    #makes a list of the energies of the chosen particle for an event (for histogram plot)
    def energyList(self, pdgID, eventID, energy_bounds = None):
        #loop through trajectories and make a list of initial energies of chosen particle for given event
        energyList = []
        
        eventMask = self.trajectories['eventID'] == eventID
        eventTrajectories = self.eventTrajectories
        pdgMask = eventTrajectories['pdgId'] == pdgID
        pdgTrajectories = eventTrajectories[pdgMask]
        
        for traj in pdgTrajectories:
            
            #find energy in this segment
            m = Particle.from_pdgid(pdgID).mass
            if type(m) == type(None): #its a neutrino
                m = 0 #approximate, for GeV neutrinos. This is just to get the energy of the primary 
            p_abs = np.linalg.norm(traj['pxyz_start']) #MeV

            #include mass energy if the particle is a lepton whose parent was a photon (would be from pair production)
            if PDGID(pdgID).is_lepton:
                parentID = traj['parentID']
                parentTraj = eventTrajectories[parentID]
                if parentTraj['pdgId'] == 22:
                    Energy = np.sqrt(m**2 + p_abs**2) 
                else:
                    Energy = np.sqrt(m**2 + p_abs**2) - m
            else:
                Energy = np.sqrt(m**2 + p_abs**2) - m
            
            #skip if not within bounds
            if energy_bounds:
                if Energy < energy_bounds[0] or Energy > energy_bounds[1]:
                    continue
            #append energy
            energyList.append(Energy)
        return energyList            
        
    
    
        
    
                
#computes lists of true energies, reconstructed energies and particle importances for some pID -- data to be fed into energy_particle_plot

#charge_readout = 'Gampix-Pixels' or 'Gampix-Tiles' or 'LArPix' or 'Anode-Grid' or a list of these
#importance_method = 'initial' (finds fraction by initial energies of particles) 
# or 'deposited' (finds fraction by integrating deposited charge from particles)
# or 'count' (finds fraction by counting number of particles)
# or 'eventID' (just catalogues by eventID)
# fraction = True means answer will be given as a fraction of total (unless method = 'count')  

def energy_vs_particle(input_edepsim_files, #input files to process
                         pID, #particle ID for which to compute importance parameters for each event
                         events, #2d array of events with the first axis giving the file number and the second axis giving the event number, or a 1d list which will be applied to all files
                         importance_method = ['eventID', 'initial', 'count', 'deposited'], #methods for finding importance of particle given by pID
                         charge_readout = ['LArPix', 'Gampix-Pixels', 'Gampix-Tiles'], #charge readouts to process
                         remove_attenuation = True, #factor out attenuation when reconstructing
                         recombination_method = consts.BIRKS, #factor out recombination when reconstructing
                         event_type = True, #write event type (charged vs neutral) to output file (only does something for neutrino files)
                         depth = 0, #depth of charge readout
                         pdgMask = None, #mask all segments not belonging to this particle ID to specifically look at the contribution coming from this particle
                         allEvents = False, #if true, will use all events in the file
                         randomize = False, #randomize the edepsim data positions and momenta (should be true for GENIE data)
                         simulation_bounds = [[-29000,29000], [-7250,7250], [-12000,0]], #bounds of simulation in mm
                         importance_energy_bounds = None,  #energy bounds for importance_method = "initial" or "deposited" or "count"
                         ):
    #starttime
    starttime = time.time()
    
    #get data in the right format
    if type(charge_readout) == str:
        charge_readout = [charge_readout]  
    if type(importance_method) == str:
        importance_method = [importance_method]  
    if type(input_edepsim_files) != list:
        if os.path.isdir(input_edepsim_files):
            contained_files = os.listdir(input_edepsim_files)
            input_edepsim_files = [os.path.join(input_edepsim_files, file) for file in contained_files]
        else:
            input_edepsim_files = [input_edepsim_files]       

    #check pID has a valid form
    assert(type(pID) == str or type(pID) == int), 'pID must be a string or int'
    if type(pID) == str:
        pID = int(Particle.from_name(pID).pdgid) 
    
    #apply pdgMask if given
    if pdgMask:
        assert(type(pdgMask) == str or type(pdgMask) == int), 'pdgMask must be a string or int'
        if type(pdgMask) == str:
            pdgMask = int(Particle.from_name(pdgMask).pdgid) 
    
    #check events has a valid form
    assert type(events) == list or type(events) == np.ndarray, 'events must be a list or numpy array'
    events = np.array(events) 
    if events.shape[0] == 1 or events.ndim == 1:
        events = np.tile(events, (len(input_edepsim_files),1))
    elif events.shape[0] != len(input_edepsim_files):
        raise ValueError('events must be a list of length 1 or the same length as input_edepsim_files')

    #Initialize EnergyResolution objects
    ERs     = [EnergyResolution(input_edepsim_files[n], events[n], depth = depth, allEvents = allEvents) for n in range(len(input_edepsim_files))]
    
    # initialize lists to fill with energies and importances
    energy_true_list = []
    deposited_energy_list = []
    deposited_charge_list = []
    recombination_list = []
    type_list = []
    importance_dict  = {importance_type:[] for importance_type in importance_method} #one list per importance method
    readout_dict = {readout_type:[] for readout_type in charge_readout} #one list per readout
    
    #loop through files
    print('Beginning reconstruction: t={0:.2f}s'.format(time.time()-starttime))
    for n, ER in enumerate(ERs):
        print('new file processing')
        
        # append true energy for each event
        if pdgMask:
            energy_true_dict  = ER.particleImportance('deposited', pdgMask,
                                               energy_bounds = importance_energy_bounds,
                                               fraction = False,
                                               removePrimaries = True)
            energy_true_instance = [energy_true_dict[k] for k in events[n]]
            energy_true_list += energy_true_instance
        else:
            energy_true_instance = ER.trueEnergy(dictionary=False)
            energy_true_list += energy_true_instance
            
        # append deposited energy and charge for each event
        print('Finding true deposited energy: t={0:.2f}s'.format(time.time()-starttime))
        deposited = ER.trueDeposited(dictionary=False)
        deposited_energy_list += deposited['deposited_energy']
        deposited_charge_list += deposited['deposited_charge']
        
        # append recombination factor for each event
        print('Finding recombination factors: t={0:.2f}s'.format(time.time()-starttime))
        recombination_instance = ER.eventRecombination(recombination_method = recombination_method, dictionary=False)
        recombination_list += recombination_instance
        
        #append particle importances
        for importance_type in importance_method:
            if importance_type == 'eventID': #catalogue by eventID
                importance_dict[importance_type] += list(events[n])
            else: # catalogue importance of particle for each event
                print('Finding particle importances: ' + importance_type + ', t={0:.2f}s'.format(time.time()-starttime))
                importance_event_dict = ER.particleImportance(
                                                importance_type, pID,
                                                energy_bounds = importance_energy_bounds, 
                                                )  
                importance_dict[importance_type] += [importance_event_dict[k] for k in events[n]]
            
        # append reconstructed energy for each readout   
        print('reconstructing. t={0:.2f}s'.format(time.time()-starttime))
        reco_event_dict       = ER.reconstruct(      
                                            remove_attenuation    = remove_attenuation,
                                            recombination_method  = recombination_method,
                                            charge_readout        = charge_readout,
                                            dictionary            = False,
                                            randomize             = randomize,
                                            simulation_bounds     = simulation_bounds,
                                            pdgMask               = pdgMask,
                                            starttime             = starttime, )
        
        for readout_type in charge_readout:     
            readout_dict[readout_type] += reco_event_dict[readout_type]          
        
        # append event characteristic (charged vs neutral) for each readout  
        if event_type and hasattr(ER, 'eventPrimaries'):
            eventPrimaries = ER.eventPrimaries
            type_dict = {eventPrimaries[n]['eventID']:eventPrimaries[n]['code'] for n in range(len(eventPrimaries))}
            type_list += [type_dict[k] for k in events[n]]
            
    #create output dictionary
    output_dict = {'pID':pID, 
                   'energy_true':energy_true_list, 
                   'deposited_energy':deposited_energy_list, 
                   'deposited_charge':deposited_charge_list,
                   'recombination':recombination_list,
                   **importance_dict, 
                   **readout_dict}
    
    output_dict['importance_energy_bounds'] = importance_energy_bounds
    output_dict['remove_attenuation']       = int(remove_attenuation)
    
    if pdgMask:
        output_dict['pdgMask'] = pdgMask
    else:
        output_dict['pdgMask'] = 0
        
    if event_type:
        output_dict['type'] = type_list

    return output_dict
            
            

def main(args):
    if os.path.isdir(args.input_directory):
        input_edepsim_files = os.listdir(args.input_directory)
        input_edepsim_files = [os.path.join(args.input_directory, file) for file in input_edepsim_files]
    else:
        input_edepsim_files = [args.input_directory]
    pID = args.pID
    charge_readout = args.charge_readout
    importance_method = args.importance_methods
    
    events = args.events
    if events ==-1: #use all events
        h5_Files = [h5py.File(input_edepsim_files[n], 'r') for n in range(len(input_edepsim_files))]
        events = [np.unique(h5_Files[n]['segments']['eventID']) for n in range(len(h5_Files))]
        allEvents = True
    elif len(events) == 2:
        events = [np.arange(events[0], events[1])]
        allEvents = False
    else:
        allEvents = False

    remove_attenuation = not args.keep_attenuation
    pdgMask = args.pdgMask
    energy_bounds = args.energy_bounds
    randomize = args.randomize
    write_event_type = args.write_event_type

    
        
    output_dict = energy_vs_particle(input_edepsim_files, 
                         pID, 
                         events,
                         importance_method, 
                         charge_readout,
                         allEvents = allEvents,
                         event_type = write_event_type,
                         remove_attenuation = remove_attenuation,
                         recombination_method = consts.BIRKS,
                         depth = 0,
                         pdgMask = pdgMask,
                         randomize = randomize,
                         importance_energy_bounds = energy_bounds, 
                         )
    
    output_file = h5py.File(args.output_filename + '.h5', 'w')
    output_file.update(output_dict)
    
    output_file.close()
    
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input_directory', type = str,
                        help = "input edep-sim dumped h5, or directory of these, to process")
    parser.add_argument('-o', '--output_filename', type = str,
                        help = "output h5 file path")
    parser.add_argument('-pid', '--pID', type = int,
                        help = 'particle ID to process')
    parser.add_argument('-e', '--events', nargs='+', type = int,
                        default = -1,
                        help = "events to process. if just 2 numbers, will use all events between those 2 indices. default is all")
    parser.add_argument('-cr', '--charge_readout', default=['Gampix-Pixels', 'Gampix-Tiles', 'LArPix', 'Anode-Grid'], nargs='+',
                        help="charge readouts to process. One or multiple of 'Gampix-Pixels', 'Gampix-Tiles', 'LArPix', or 'Anode-Grid'. default is all 3")
    parser.add_argument('-im', '--importance_methods', default = ['count', 'initial', 'deposited', 'eventID'], nargs='+',
                        help="methos for finding importance of particle given by pID. default is all, can be any of 'count','initial', 'deposited', 'eventID' (last just catalogues by eventID and ignores importance)")
    parser.add_argument('-ka', '--keep_attenuation', action='store_true',
                        default = False,
                        help = "don't factor out attenuation when reconstructing")
    parser.add_argument('-pm', '--pdgMask', type = int, default= None,
                        help = 'str, particle ID mask to ignore everything else, like "gamma" or "e-"')
    parser.add_argument('-eb', '--energy_bounds', nargs='+', type = float, default = None,
                        help = 'energy bounds for importance_method = "initial" or "deposited"')   
    parser.add_argument('-r', '--randomize', action='store_true',default = False,
                        help = "randomize the edepsim data positions and momenta (should be true for GENIE data)") 
    parser.add_argument('-et', '--write_event_type', action='store_true', default = True,
                        help = "write event type (charged vs neutral) to output file")
    args = parser.parse_args()

    main(args)
