#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 22:25:33 2022

Colletion of tools to read .sim file and crate events

Parsing has long list of problems / needed improvements - see comments

@author: tshutt
"""

class Sim_File:
    """"
    Class for sim file handling

    At initialization, opens file and reads up to first event.

    Subsequent calls to read_next_event reads the next event and returns it
    in a raw form as self.raw_event.  parse_raw_event creates hits
    from the events in the .ht and .ia lines, returning self.parsed_event

    make_event_tree and show_event_tree are for diagnostics

    Units from Megablib: cm and keV
    """

    def __init__(self, full_file_name):
        """ Opens a sim file and reads to the first event.
            Returns file "handle" f

            full_file_name includes path, but not extension

            TODO: parse and return header information
            8/20 TS
        """
        #   Open
        f = open(full_file_name + '.sim')

        #   Skip to first event - text_line = 'SE'
        text_line = f.readline().strip('\n')

        while text_line=="" or text_line[0:2]!='SE':
            text_line = f.readline().strip('\n')
            print(text_line)

        self.f = f

    def read_next_event(self):
        """
        Reads the next event in the file, and returns all of
        the raw information in the dictionary self.raw_event

        Assumes that the file has previously be read up to and including
        the 'SE' line demarking the start of this event, and reads the
        next event 'SE' before returning.  If no next 'SE' is encountered,
        thensim['end_of_file'] is set True.
        """

        import sys
        import numpy as np
        from collections import deque

        #   event_info is essentially meta data for the event
        event_info = {}

        #   Event ID
        text_line = self.f.readline().strip('\n')
        numstring = text_line[2:].split()
        event_info['triggered_event_id'] = int(numstring[0])
        event_info['simulated_event_id']= int(numstring[1])

        #   Event time
        text_line = self.f.readline().strip('\n')
        event_info['time'] = float(text_line[2:])

        #   Energies: active, escaped and passive energy
        text_line = self.f.readline().strip('\n')
        event_info['active_energy'] = float(text_line[2:])
        text_line = self.f.readline().strip('\n')
        event_info['escaped_energy'] = float(text_line[2:])
        text_line = self.f.readline().strip('\n')
        event_info['passive_energy'] = float(text_line[2:])

        #   Specific passive material - PM
        #   The type of material is recored, but I'm not trackign that
        #   yet
        text_line = self.f.readline().strip('\n')
        event_info['energy_special_passive'] = []
        while text_line[0:2]=='PM':
            splitline=text_line[2:].split()
            event_info['energy_special_passive'].append(float(splitline[1]))
            text_line = self.f.readline().strip('\n')

        #   The event is returned as the dictionary raw_event
        raw_event = {}
        raw_event['event_info'] = event_info

        #%%   IA - interaction lines.  First initialize things

        #   Dictionary of interaction codes
        interaction_type_code = {}
        interaction_type_code['INIT'] = 1
        interaction_type_code['PAIR'] = 2
        interaction_type_code['COMP'] = 3
        interaction_type_code['PHOT'] = 4
        interaction_type_code['BREM'] = 5
        interaction_type_code['ANNI'] = 6
        interaction_type_code['RAYL'] = 7
        interaction_type_code['IONI'] = 8
        interaction_type_code['INEL'] = 9
        interaction_type_code['CAPT'] = 10
        interaction_type_code['DECA'] = 11
        interaction_type_code['ESCP'] = 12
        interaction_type_code['ENTR'] = 13
        interaction_type_code['EXIT'] = 14
        interaction_type_code['BLAK'] = 15

        ia = {}

        #   Initialize values - appended in loop below
        ia['interaction_type'] = []
        ia['interaction'] = deque()
        ia['parent_interaction'] = deque()
        ia['sub_detector'] = deque()
        ia['time'] = deque()
        ia['rx'] = deque()
        ia['ry'] = deque()
        ia['rz'] = deque()
        ia['particle'] = deque()
        ia['s_primary_x'] = deque()
        ia['s_primary_y'] = deque()
        ia['s_primary_z'] = deque()
        ia['polarization_in_x'] = deque()
        ia['polarization_in_y'] = deque()
        ia['polarization_in_z'] = deque()
        ia['energies_primary'] = deque()
        ia['particle_out'] = deque()
        ia['s_secondary_x'] = deque()
        ia['s_secondary_y'] = deque()
        ia['s_secondary_z'] = deque()
        ia['polarization_out_x'] = deque()
        ia['polarization_out_y'] = deque()
        ia['polarization_out_z'] = deque()
        ia['energies_secondary'] = deque()

        interaction_code = deque()

        #%%   Read IA lines and put into ia structure
        while text_line[0:2]=='IA':

            ia['interaction_type'].append(text_line[3:7])
            splitline = text_line[8:].split(';')

            ia['interaction'].append(int(splitline[0]))

            ia['parent_interaction'].append(int(splitline[1]))
            ia['sub_detector'].append(int(splitline[2]))
            ia['time'].append(float(splitline[3]))
            ia['rx'].append(float(splitline[4]))
            ia['ry'].append(float(splitline[5]))
            ia['rz'].append(float(splitline[6]))
            ia['particle'].append(int(splitline[7]))
            ia['s_primary_x'].append(float(splitline[8]))
            ia['s_primary_y'].append(float(splitline[9]))
            ia['s_primary_z'].append(float(splitline[10]))
            ia['polarization_in_x'].append(float(splitline[11]))
            ia['polarization_in_y'].append(float(splitline[12]))
            ia['polarization_in_z'].append(float(splitline[13]))
            ia['energies_primary'].append(float(splitline[14]))
            ia['particle_out'].append(int(splitline[15]))
            ia['s_secondary_x'].append(float(splitline[16]))
            ia['s_secondary_y'].append(float(splitline[17]))
            ia['s_secondary_z'].append(float(splitline[18]))
            ia['polarization_out_x'].append(float(splitline[19]))
            ia['polarization_out_y'].append(float(splitline[20]))
            ia['polarization_out_z'].append(float(splitline[21]))
            ia['energies_secondary'].append(float(splitline[22]))

            #   Interaction code
            interaction_code.append(
                interaction_type_code[ia['interaction_type'][-1]]
                )

            #   Read next text_line
            text_line = self.f.readline().strip('\n')

            #   If eof or other file error, break
            if text_line=='':
                break
        # Convert deques to numpy arrays
        for key in ia:
            ia[key] = np.array(ia[key])

        raw_event['ia'] = ia
        raw_event['interaction_code'] = np.array(interaction_code)

        #%%   Read HT (hit) lines and put into ht structure

        #   If eof or other error, return
        if text_line != '':

            ht = {}

            #   Initialize values - appended in loop below
            ht['detector'] = deque()
            ht['cell'] = deque()
            ht['rx'] = deque()
            ht['ry'] = deque()
            ht['rz'] = deque()
            ht['energy'] = deque()
            ht['time'] = deque()
            ht['interaction'] = deque()

            front_acd = False
            back_acd = False

            #   Loop over HTsim lines
            n = 0
            while text_line[0:2] == 'HT':

                n += 1

                #   Split line into text sections
                splitline = text_line[6:].split(';')

                #   Detector ID
                detector_id = int(splitline[0])
                ht['detector'].append(detector_id)

                #   Some ugliness here for Oliver's hack (documented at:
                #   https://confluence.slac.stanford.edu
                #       /display/GTPC/MEGALib+notes)
                #   which added a volume tag to the .sim file.  Here
                #   Check for volume tag added to our version of the .sim
                #   file by O. Hitchcock, as documented here:

                #   If detector_id is 4, then treat as ACD, and ignore
                #   what should be the volume tag but is currently (6/23)
                #   broken.   This all should be fixed and generalized.
                #   For now, do not record any details about hits in
                #   ACDS - energy, time, different hits in front and back
                if detector_id==4:
                    if float(splitline[4])>0:
                        front_acd = True
                    else:
                        back_acd = True

                #   TPC cells
                elif detector_id==5:
                    ht['cell'].append(int(splitline[1].split('_')[1]))
                    ht['rx'].append(float(splitline[2]))
                    ht['ry'].append(float(splitline[3]))
                    ht['rz'].append(float(splitline[4]))
                    ht['energy'].append(float(splitline[5]))
                    ht['time'].append(float(splitline[6]))
                    ht['interaction'].append(int(splitline[7]))

                else:
                    sys.exit('Error in read_next_event: bad detector_id')

                text_line = self.f.readline().strip('\n')
                if text_line == '':
                    break

            #   Convert ht structure to numpy
            for key in ht:
                ht[key] = np.array(ht[key])

            #   Assign to sim structure
            if len(ht['detector'])>0:
                raw_event['ht'] = ht
                raw_event['ht']['front_acd'] = front_acd
                raw_event['ht']['back_acd'] = back_acd

            #   Add incident energy from IA line
            raw_event['event_info']['incident_energy'] = \
                raw_event['ia']['energies_secondary'][0]

        #   Check for end of file.  If text_line is 'SE', then no
        if text_line[0:2] == 'SE':
            end_of_file = False
        else:
            end_of_file = True
        raw_event['end_of_file'] = end_of_file

        # meta data
        raw_event['meta'] = {}
        raw_event['meta']['interaction_codes'] = interaction_type_code

        # return raw_event
        self.raw_event = raw_event

    def parse_raw_event(self):

        """
        Parses the raw event information in self.raw_event, returning
        self.parsed_event.

        Adds up energy in HT lines, applying this to
        interactions in IA lines.  Currently adds one layer of
        clustering - identifies "x-rays" by energy only, and adds
        them to parent.  This requires special handling when x-rays
        created during bremstrahlung.
        Also generates event quality booleans

        Result is parsed_event, which is a list of only energy and initial
        locations of what ideally are electron tracks (with associated
        x-ray energies).  This is ready for vectoriztion to apply detector
        response

        Cell information:
            hit_cell - cell # for each hit, size(max_num_hits, num_events)
            cells_list - list of struck cells
            num_cells - number of sturck cells,
            hit_cell_index- hit index for cells_list
            cell_hit_count - number of hits in each cell

        Note the spatial information in HT lines and IA lines are not
            consistent - HT lines are in cell coordinates.
            TODO: Should clean this up.

        TODO: Many issues to clean up:
            "Compton" interactions that have an x-ray energy and occur at
                same location as previous - currently counted as separate
                hit, but should not be.
                Related - x-ray shows up at end of IA list and appears as
                separte hit at same location
            Events with inactive hits followed by energy deposition.
                Need to find active energy and assign to hit.
                any issues at cell boundaries?
                Add flag for first hit being like this - if calculating
                true angle?
            Two input parameters need to be handled: energy
                threshold for "x-rays" and distance for clustering.
                Probably not in params, but some other settings
            Clustering should probably also include check on
                distance from parent.  Will need to fix HT/IA lines
                position problem to do this.
            Mismatch of total energy and .sim active energy in ~10^-3 events
            The logic used here is a bit compliected, but still
                too simple, and not very satisfactory - a more
                general framework based on distances alone seems be better.
        """

        import numpy as np

        if not 'ht' in self.raw_event:
            print('*** Error in parse_raw_event - no ht lines ***')
            return

        #   THIS SHOULD BE SET ELSEWHERE
        clustering_energy_threshold = 4.0

        #   Convenient variable names.  These include the first
        #   ia element, which is not an interaction.
        #   This list is needed when referencing by parent_ia, since
        #   this will include the firt element.
        ia_all = self.raw_event['ia']['interaction']
        parent_ia_all = self.raw_event['ia']['parent_interaction']
        interaction_code_all = self.raw_event['interaction_code']
        sub_detector_all = self.raw_event['ia']['sub_detector']

        #   These start with the second ia element, and thus have the
        #   dimension of the event we are constructing to return
        #   Reference these directly
        #   TODO clean this comment up.  wtf?
        ia = ia_all[1:]
        parent_ia = parent_ia_all[1:]
        interaction_code = interaction_code_all[1:]
        sub_detector = sub_detector_all[1:]

        r = np.zeros((3, len(self.raw_event['ia']['rx'])-1))
        r[0, :] = self.raw_event['ia']['rx'][1:]
        r[1, :] = self.raw_event['ia']['ry'][1:]
        r[2, :] = self.raw_event['ia']['rz'][1:]

        s_primary = np.zeros((3, len(self.raw_event['ia']['rx'])-1))
        s_primary[0, :] = self.raw_event['ia']['s_primary_x'][1:]
        s_primary[1, :] = self.raw_event['ia']['s_primary_y'][1:]
        s_primary[2, :] = self.raw_event['ia']['s_primary_z'][1:]

        s_secondary = np.zeros((3, len(self.raw_event['ia']['rx'])-1))
        s_secondary[0, :] = self.raw_event['ia']['s_secondary_x'][1:]
        s_secondary[1, :] = self.raw_event['ia']['s_secondary_y'][1:]
        s_secondary[2, :] = self.raw_event['ia']['s_secondary_z'][1:]

        #   Number of Rayleigh scatters, and number of bremstrahlung
        num_rayleigh = np.sum(interaction_code==7)
        num_brem = np.sum(interaction_code==5)

        #   Loop over interactions and:
        #   + Assign ht energy to parent interactions
        #   + Assign cells to interactions
        energy = np.zeros(len(interaction_code))
        hit_cell = np.zeros(len(interaction_code),dtype=int)
        for ni in range(len(ia)):
            hits = self.raw_event['ht']['interaction']==ia[ni]
            if np.any(hits):
                all_cells = self.raw_event['ht']['cell'][hits]
                if len(np.unique(all_cells))>1:
                    print(f'Error: multiple cells, ia = {ia[ni]:2d}')
                hit_cell[ni] = all_cells[0]
                energy[ni] = sum(self.raw_event['ht']['energy'][hits])

        #   Deem a hit an "x-ray" if:
        #      + energy less than pre-defined threshold
        #      + photoabsorption
        #      + interaction is not in passive material
        #   In principle, should also add a check that the daughter is
        #       close to the location of parent, but not implemented yet
        x_ray = \
            (energy<clustering_energy_threshold) \
            & (interaction_code==4) \
            & (sub_detector!=0)

        #   A correction when x_ray with parent==1 - i.e.,
        #   the incoming gamma.
        #   This is a sub-class of events that have two consecutive photo
        #   interactions with parent = 1,
        #   the second of which is an apparent x-ray by energy
        #   Assign the previous interaction as the parent of these x-rays
        parent_ia[x_ray & (parent_ia==1)] = ia[x_ray & (parent_ia==1)] - 1

        #   A correction when x_ray with parent==5 - bremstrahlung.
        #   Assign the grandparent as the parent of these x-rays
        parent_ia[x_ray&(interaction_code_all[parent_ia-1]==5)] = \
            parent_ia_all[
                parent_ia[x_ray&(interaction_code_all[parent_ia-1]==5)]-1
                ]

        #   Throw up notice if x_ray, but parent passive.
        #   Can in principle happen if parent at edge of cell.
        #   Leave this for diagnositcs, could/should remove later.
        if any(x_ray
               & (sub_detector_all[parent_ia-1]==0)):
            print(
                'Note: x-ray in active region from passive parent hit'
                + ', event '
                + str(self.raw_event['event_info']['triggered_event_id'])
                )

        #   Add x-ray energy to parent, and assign x_ray cell to
        #   parent if parent = zero
        for ni in np.unique(parent_ia[x_ray])[::-1]:
            energy[ia==ni] += \
                np.sum(energy[x_ray&(parent_ia==ni)])

        #   Keep track of escaped energy depending on direction
        #   up or down (through or back-scatter, assuming planar geometry)
        escaped = interaction_code==12
        back = (
            self.raw_event['ia']['s_secondary_z'][0]
            * (s_primary[2][escaped]) < 0
            )
        escaped_back_energy = np.sum(
            self.raw_event['ia']['energies_primary'][1:][escaped][back]
            )
        escaped_through_energy = np.sum(
            self.raw_event['ia']['energies_primary'][1:][escaped][~back]
            )

        #   Flag "interactions" to be removed
        #       0  - inert material
        #       5  - bremstahlung
        #       6  - annihilation
        #       7  - rayleigh
        #       12 - escape
        remove = \
            x_ray \
            | (self.raw_event['ia']['sub_detector'][1:]==0) \
            | (interaction_code==5) \
            | (interaction_code==6) \
            | (interaction_code==7) \
            | (interaction_code==12) \

        #   Remove them, leaving only good "hits"
        ia = ia[~remove]
        interaction_code = interaction_code[~remove]
        energy = energy[~remove]
        hit_cell = hit_cell[~remove]
        r = r[:, ~remove]
        s_primary = s_primary[:, ~remove]
        s_secondary = s_secondary[:, ~remove]

        """
        THESE are old notes.  but should probably be thought about.
        To deal with here:
            comp - check distance, then add to predecesor (if any)
             - problem with event 1 - apparent x-ray is labeled as compton
            phot - check distance, then add to parent or not
            anni - should I check on energy = zero
        """

        #   List of cells with hits (length is <= number of hits)
        #   and index into this list for each hit (length = number of hits)
        cells_list, hit_cell_index, cell_hit_count \
            = np.unique(hit_cell, return_inverse=True, return_counts=True)

        #   Number of cells hit per event
        num_cells = cells_list.size

        #   Incident s is secondary s in first ia line
        s_incident = np.zeros((3,))
        s_incident[0] = self.raw_event['ia']['s_secondary_x'][0]
        s_incident[1] = self.raw_event['ia']['s_secondary_y'][0]
        s_incident[2] = self.raw_event['ia']['s_secondary_z'][0]

        #   Check on energy.  This is currently finding errors
        #   at the level of a few per 1000.  Those should be figured out.
        if not (self.raw_event['event_info']['active_energy']
                - sum(energy)) < 1e-3:
            print('Error: missing energy, event '
                  + str(self.raw_event['event_info']['triggered_event_id']))

        #   Event quality additional information

        #   Initial scatter in inert material or not
        if sub_detector[0] == 0:
            clean_entrance = False
            entrance_scatter_energy \
                = self.raw_event['ia']['energies_secondary'][0]
        else:
            clean_entrance = True
            entrance_scatter_energy = 0

        #   Missing energy - escaped or passive, above 1 meV
        if (self.raw_event['event_info']['escaped_energy']
            + self.raw_event['event_info']['passive_energy']) > 0.001:
            missing_energy = True
        else:
            missing_energy = False

        #   Missing energy after 1st scatter, independent of
        #   whether first scatter was in inert material (i.e., whether
        #   there was an  entrance scatter)
        if (np.any(sub_detector[1:]==0)
            | (self.raw_event['event_info']['escaped_energy'] > 0.001)):
            missing_energy_after_entrance = True
        else:
            missing_energy_after_entrance = False

        #   Package output
        parsed_event = {}

        #   Arrays of per hit information
        parsed_event['energy'] = energy
        parsed_event['r'] = r
        parsed_event['s_primary'] = s_primary
        parsed_event['s_secondary'] = s_secondary
        parsed_event['hit_cell'] = hit_cell
        parsed_event['hit_cell_index'] = hit_cell_index
        parsed_event['cells_list'] = cells_list
        parsed_event['num_cells'] = num_cells
        parsed_event['cell_hit_count'] = cell_hit_count
        parsed_event['interaction_code'] = interaction_code

        #   Per event information
        parsed_event['incident_energy'] \
            = self.raw_event['event_info']['incident_energy']
        parsed_event['s_incident'] = s_incident
        parsed_event['clean_entrance'] = clean_entrance
        parsed_event['entrance_scatter_energy'] = entrance_scatter_energy
        parsed_event['missing_energy'] = missing_energy
        parsed_event['missing_energy_after_entrance'] \
            = missing_energy_after_entrance
        parsed_event['num_rayleigh'] = num_rayleigh
        parsed_event['num_brem'] = num_brem
        parsed_event['escaped_back_energy'] = escaped_back_energy
        parsed_event['escaped_through_energy'] = escaped_through_energy

        self.parsed_event = parsed_event

    def make_event_tree(self):
        """
        Constructs a tree diagram of event

        8/24/2020   TS
        """

        import numpy as np

        tree = {}

        if not 'ht' in self.raw_event:
            print('*** Warning in DiagramEvent:'
                  + 'no HT lines, event skipped ***')
            self.tree = None
            return

        #   Add hits to interactions.
        tree['num_hits'] = np.empty(
            len(self.raw_event['ia']['interaction']),dtype=int)
        tree['hit_energy'] \
            = np.empty(len(self.raw_event['ia']['interaction']))
        for ni in range(len(self.raw_event['ia']['interaction'])):
            hits = self.raw_event['ht']['interaction']==(ni+1)
            tree['num_hits'][ni] = sum(hits)
            tree['hit_energy'][ni]= sum(self.raw_event['ht']['energy'][hits])

        #   These next section creates all_hit_energies and all_secondries.
        #   For each interaction these are lists of all secondary energies
        #   and interactions, including the interaction itself

        #   Start with list of interactions and associated parents
        interactions = self.raw_event['ia']['interaction'].tolist()
        parents = self.raw_event['ia']['parent_interaction'].tolist()

        #   Prepopulate the "all" fields with interactions/energies for this
        #   interaction
        tree['all_hit_energies'] = []
        tree['all_secondries'] = []
        for ni in range(len(interactions)):
            tree['all_secondries'].append([ni+1])
            tree['all_hit_energies'].append([tree['hit_energy'][ni]])

        #   This iterative procedure starts at the "bottom" of the event  -
        #   that is, interactions with no secondary, and works its way up,
        #   at each iteraction removing the bottom set of interactions,
        #   and so working back to the primary
        while len(interactions)>1:

            parent_list=set(parents)

            #   "bottom" of the chain is all interactions that are not
            #   a parent
            bottoms=[]
            for ni in range(len(interactions)):
                if set({interactions[ni]}).isdisjoint(parent_list):
                    bottoms.append(interactions[ni])

            #   Assign these bottom interactions to their parents
            these_parents = \
                self.raw_event['ia']['parent_interaction'] \
                    [np.array(bottoms)-1]
            bottom_cut = np.array(np.zeros(len(interactions)),dtype=bool)
            for nb in range(len(bottoms)):
                for nnb in range(len(tree['all_secondries'][bottoms[nb]-1])):
                    tree['all_secondries'][these_parents[nb]-1].append(
                        tree['all_secondries'][bottoms[nb]-1][nnb])
                bottom_cut[np.array(interactions)==(bottoms[nb])] = True

            #   Remove current bottom interactions from both
            #   interactions and parents
            interactions = \
                np.ndarray.tolist(np.array(interactions)[~bottom_cut])
            parents = \
                np.ndarray.tolist(np.array(parents)[~bottom_cut])


        #   Sort list of all interactions for convenience, then assign
        #   hit energies associated with each secondary
        for ni in range(len(tree['all_secondries'])):
            tree['all_secondries'][ni].sort()
            tree['all_hit_energies'][ni] = tree['hit_energy'][
                np.array(tree['all_secondries'][ni])-1]

        #   Create direct_secondaries, and direct_hit_energies.
        #   For each interaction,
        #   this is a list of only those interactions and eneriges that are
        #   direct secondaries of that interaction.
        tree['direct_secondaries'] = \
            [[] for i in range(len(self.raw_event['ia']['interaction']))]
        tree['direct_hit_energies'] = \
            [[] for i in range(len(self.raw_event['ia']['interaction']))]
        for ni in range(len(self.raw_event['ia']['interaction'])):
            if self.raw_event['ia']['parent_interaction'][ni]>0:
                tree['direct_secondaries'][
                    self.raw_event['ia']['parent_interaction'][ni]-1
                    ].append(ni+1)
                tree['direct_hit_energies'][
                    self.raw_event['ia']['parent_interaction'][ni]-1
                    ].append(tree['hit_energy'][ni])

        self.tree = tree

    def show_event_tree(self, startindex=0):

        """ Displays tree to command line """

        #   Recursive display of direct secondary information
        def dig_down(event, startindex, indentcounter):
            """ Called by show_event_tree, displays tree information.
            Calls itself recursively """

            c = ' ' * 3 * indentcounter

            for nn in range(len(
                    event.tree['direct_secondaries'][startindex])):

                ns=event.tree['direct_secondaries'][startindex][nn]
                nsi = ns - 1

                #   Need this tag
                if event.sim['ia']['sub_detector'][nsi]==5:
                    activetag='active '
                else:
                    activetag='passive'

                # if ia['type{nsi}=='COMP'
                #     or ia['type{nsi}=='BREM'
                #     or ia['type{nsi=='ANNI':
                #                iaetag='IA E_secondary = ' ...
                #         sprintf('#4.1f',ia['energies.out(nsi)) ', ']
                # else
                #     iaetag=[]

                print(
                    '%s'
                    'Int. %d. '
                    '%s, '
                    '%s, '
                    'E_hits: Tot = %4.2f, '
                    'This int = %4.2f'
                    % (
                        c,
                        ns-1,
                        event.sim['ia']['interaction_type'][nsi],

                        activetag,
                        # iaetag ...
                        sum(event.tree['all_hit_energies'][nsi]),
                        event.tree['hit_energy'][nsi]
                        )
                    )

                dig_down(event, nsi, indentcounter+1)

        #   Make tree if it doesn't exist
        if not 'tree' in self.raw_event['event_info']:
            self.make_event_tree()

        #   Bail if still no tree.  Should really have
        #   set error in make_event_tree
        if not self.tree:
            print('Bad event - no tree created')
            return

        #   Event num and event-level summed energies
        print('Event %d, keV: '
              'incident =  %4.1f, '
              'tot =  %4.1f, '
              'escaped =  %4.1f, '
              'passive = %4.1f'
               % (
               self.raw_event['event_info']['triggered_event_id'],
               self.raw_event['event_info']['incident_energy'],
               self.raw_event['event_info']['active_energy'],
               self.raw_event['event_info']['escaped_energy'],
               self.raw_event['event_info']['passive_energy']
               ))

        #   Other summary infomration
        # print(
        #     '  IA %d '
        #     '%s, '
        #     'Hit energies: this interaction = %4.2f, '
        #     'this + secondaries = %4.2f'
        #     % (
        #         self.raw_event['ia['interaction[startindex],
        #         self.raw_event['ia['interaction_type[startindex],
        #         self.tree['hit_energy[startindex],
        #         sum(self.tree['all_hit_energies[startindex])
        #         ))

        #   Recursively blab about secondaries
        dig_down(self, startindex, 1)

def read_events_from_sim_file(full_file_name,
                              geo_params,
                              num_events_to_read=1e10
                              ):
    """
    Opens .sim file and reads and parses all events to create hits, which
    are (ideally) separate electron recoil tracks.

    These are then put into a "flattened" set of arrays for vectorized
    calculations.  This is returned as the diciontaries 'truth', which is
    quantities of (mostly) dimension [num_events], and 'truth_hits',
    which is quanities of dimension [num_hits, num_events].
    %TODO clean up the poorly thought out use of copy, and replace the
        matrix representation with awkward structures
    """

    from . import geometry_tools
    from .math_tools import dot

    import numpy as np
    import copy

    #   Prep lists
    time_list = []
    triggered_id_list = []
    incident_energy_list = []
    s_incident_list = []
    active_energy_list = []
    escaped_energy_list = []
    passive_energy_list = []
    clean_entrance_list = []
    entrance_scatter_energy_list = []
    missing_energy_list = []
    missing_energy_after_entrance_list = []
    sum_hits_energy_list = []
    entrance_scatter_energy_list = []
    energy_list = []
    r_list = []
    s_primary_list = []
    s_secondary_list = []
    hit_cell_list = []
    hit_cell_index_list = []
    cells_list_list = []
    cell_hit_count_list = []
    num_cells_list = []
    num_hits_list = []
    interaction_code_list = []
    num_rayleigh_list = []
    num_brem_list = []
    escaped_back_energy_list = []
    escaped_through_energy_list = []

    #   Open file
    sim_file = Sim_File(full_file_name)

    # Read all events, put hit variables for each into long list
    for n in range(round(num_events_to_read)):

        #   Load next event
        sim_file.read_next_event()

        #   This skips events with no active energy deposited.
        #   Probably need to do something with skipped events
        if 'ht' in sim_file.raw_event:

            #   Generate truth_hits for this event
            sim_file.parse_raw_event()

            #   Add things to long list
            time_list.append(
                copy.copy(sim_file.raw_event['event_info']['time']))
            triggered_id_list.append(
                copy.copy(
                    sim_file.raw_event['event_info']['triggered_event_id']))
            incident_energy_list.append(
                copy.copy(
                    sim_file.raw_event['event_info']['incident_energy']))
            s_incident_list.append(
                copy.copy(sim_file.parsed_event['s_incident']))
            active_energy_list.append(
                copy.copy(sim_file.raw_event['event_info']['active_energy']))
            escaped_energy_list.append(
                copy.copy(sim_file.raw_event['event_info']['escaped_energy']))
            passive_energy_list.append(
                copy.copy(sim_file.raw_event['event_info']['passive_energy']))
            clean_entrance_list.append(
                copy.copy(sim_file.parsed_event['clean_entrance']))
            entrance_scatter_energy_list.append(
                copy.copy(sim_file.parsed_event['entrance_scatter_energy']))
            missing_energy_list.append(
                copy.copy(sim_file.parsed_event['missing_energy']))
            missing_energy_after_entrance_list.append(
                copy.copy(sim_file.parsed_event[
                    'missing_energy_after_entrance'
                    ]))
            sum_hits_energy_list.append(
                np.sum(sim_file.raw_event['ht']['energy']))
            energy_list.append(
                copy.copy(sim_file.parsed_event['energy']))
            r_list.append(
                copy.copy(sim_file.parsed_event['r']))
            s_primary_list.append(
                copy.copy(sim_file.parsed_event['s_primary']))
            s_secondary_list.append(
                copy.copy(sim_file.parsed_event['s_secondary']))
            hit_cell_list.append(
                copy.copy(sim_file.parsed_event['hit_cell']))
            hit_cell_index_list.append(
                copy.copy(sim_file.parsed_event['hit_cell_index']))
            cells_list_list.append(
                copy.copy(sim_file.parsed_event['cells_list']))
            cell_hit_count_list.append(
                copy.copy(sim_file.parsed_event['cell_hit_count']))
            num_cells_list.append(
                copy.copy(sim_file.parsed_event['num_cells']))
            num_hits_list.append(
                len(sim_file.parsed_event['energy']))
            interaction_code_list.append(
                copy.copy(sim_file.parsed_event['interaction_code']))
            num_rayleigh_list.append(
                copy.copy(sim_file.parsed_event['num_rayleigh']))
            num_brem_list.append(
                copy.copy(sim_file.parsed_event['num_brem']))
            escaped_back_energy_list.append(
                copy.copy(sim_file.parsed_event['escaped_back_energy']))
            escaped_through_energy_list.append(
                copy.copy(sim_file.parsed_event['escaped_through_energy']))

        if sim_file.raw_event['end_of_file']:
            break

    #   Check that anything is present
    if not num_hits_list:
        print('Error in parse_and_flatten_sim_file: no valid events')
        return

    #   These have fixed number of values per event - only need to
    #   create array from lists
    num_hits = np.array(num_hits_list)
    num_cells = np.array(num_cells_list)
    time = np.array(time_list)
    triggered_id = np.array(triggered_id_list)
    incident_energy = np.array(incident_energy_list)
    s_incident = np.array(s_incident_list).T
    active_energy = np.array(active_energy_list)
    escaped_energy = np.array(escaped_energy_list)
    passive_energy = np.array(passive_energy_list)
    clean_entrance = np.array(clean_entrance_list)
    entrance_scatter_energy = np.array(entrance_scatter_energy_list)
    missing_energy = np.array(missing_energy_list)
    missing_energy_after_entrance \
        = np.array(missing_energy_after_entrance_list)
    entrance_scatter_energy = np.array(entrance_scatter_energy_list)
    sum_hits_energy = np.array(sum_hits_energy_list)
    num_rayleigh = np.array(num_rayleigh_list)
    num_brem = np.array(num_brem_list)
    escaped_back_energy = np.array(escaped_back_energy_list)
    escaped_through_energy = np.array(escaped_through_energy_list)

    #   Maximum nubmer of hits and cells used below
    max_hits = max(num_hits)
    max_cells = max(num_cells)

    #   Hit or cell based data with multiple entries per event
    #   need more work.  First initialize arrays

    energy = np.zeros((max_hits, len(time)))
    r = np.zeros((3, max_hits, len(time)))
    s_primary = np.zeros((3, max_hits, len(time)))
    s_secondary = np.zeros((3, max_hits, len(time)))
    hit_cell = np.zeros((max_hits, len(time)), dtype=int)
    hit_cell_index = np.zeros((max_hits, len(time)), dtype=int)
    cells_list = np.zeros((max_cells, len(time)), dtype=int)
    cell_hit_count = np.zeros((max_cells, len(time)), dtype=int)
    alive  = np.zeros((max_hits, len(time)), dtype=bool)
    interaction_code = np.zeros((max_hits, len(time)), dtype=int)

    #   Loop over events, adding to "active" part of array.
    #   Is there a vectorized way to do this?
    for ne in range(len(energy_list)):
        energy[:num_hits[ne], ne] = energy_list[ne]
        r[:, :num_hits[ne], ne] = r_list[ne]
        s_primary[:, :num_hits[ne], ne] = s_primary_list[ne]
        s_secondary[:, :num_hits[ne], ne] = s_secondary_list[ne]
        hit_cell[:num_hits[ne], ne] = hit_cell_list[ne]
        hit_cell_index[:num_hits[ne], ne] = hit_cell_index_list[ne]
        cells_list[:num_cells[ne], ne] = cells_list_list[ne]
        cell_hit_count[:num_cells[ne], ne] = cell_hit_count_list[ne]
        alive[:num_hits[ne], ne] = True
        interaction_code[:num_hits[ne], ne] = interaction_code_list[ne]

    #   Change units for r from cm to m.
    r = r / 100

    #   Construct drift distance - front cells should have positive z,
    #   back cells have negative.
    z_drift = (
        (geo_params.cells['height'] / 2
         +  geo_params.z_centers['front_cells']
         - r[2, :, :]) * geo_params.cells['front_layer'][hit_cell - 1]
        + (geo_params.cells['height'] / 2
         -  geo_params.z_centers['back_cells']
         + r[2, :, :]) * geo_params.cells['back_layer'][hit_cell - 1]
        )
    z_drift[~alive] = 0

    #   Calculate "true" theta - from geometry of incident and first
    #   scattered vectors
    theta = np.arccos(dot(s_incident, s_primary[:, 0, :]))

    #   Put everything into truth and truth_hits dictionaries, paying
    #   attention to order
    truth = {}
    truth['time'] = time
    truth['triggered_id'] = triggered_id
    truth['incident_energy'] = incident_energy
    truth['s_incident'] = s_incident
    truth['active_energy'] = active_energy
    truth['escaped_energy'] = escaped_energy
    truth['passive_energy'] = passive_energy
    truth['missing_energy'] = missing_energy
    truth['missing_energy_after_entrance'] \
        = missing_energy_after_entrance
    truth['clean_entrance'] = clean_entrance
    truth['entrance_scatter_energy'] = entrance_scatter_energy
    truth['num_hits'] = num_hits
    truth['num_rayleigh'] = num_rayleigh
    truth['num_brem'] = num_brem
    truth['escaped_back_energy'] = escaped_back_energy
    truth['escaped_through_energy'] = escaped_through_energy
    truth['theta_geometry'] = theta

    truth_hits = {}
    truth_hits['energy'] = energy
    truth_hits['total_energy'] = sum_hits_energy
    truth_hits['r'] = r
    truth_hits['z_drift'] = z_drift
    truth_hits['num_hits'] = num_hits
    truth_hits['alive'] = alive
    truth_hits['s_primary'] = s_primary
    truth_hits['s_secondary'] = s_secondary
    truth_hits['hit_cell'] = hit_cell
    truth_hits['hit_cell_index'] = hit_cell_index
    truth_hits['cells_list'] = cells_list
    truth_hits['cell_hit_count'] = cell_hit_count
    truth_hits['num_cells'] = num_cells
    truth_hits['interaction_code'] = interaction_code

    #   Generate locations in cell coordinates
    truth_hits['r_cell'] = geometry_tools.global_to_cell_coordinates(
        truth_hits['r'],
        truth_hits['hit_cell'],
        geo_params,
        alive = truth_hits['alive']
        )

    #   Meta data from the last read sim events
    meta = sim_file.raw_event['meta']

    return truth, truth_hits, meta

def write_events_file(events, full_sim_file_name):
        """
        Saves events t disk, in two .h5s and a  .picckle.

        full_file_stub includes path but not extension

        11/12/21 TS
        """

        import os
        import sys
        import deepdish as dd
        import pickle

        #   Check that file name is same as stored in meta data
        path, file_name = os.path.split(full_sim_file_name)
        if file_name!=events.meta['sim_file_name']:
            sys.exit('Error in write_events_file - output name '
                     + "does not match meta['sim_file_name']")

        #   Save meta structure as pickle
        with open(os.path.join(
                full_sim_file_name
                + '.meta' + '.pickle'
                ),
                'wb') as f:
            pickle.dump(events.meta, f)

        #   Save events.truth to .h5 file
        dd.io.save(
            os.path.join(
                full_sim_file_name
                + '.truth' + '.h5'
                ),
            events.truth)

        #   Save events.truth to .h5 file
        dd.io.save(
            os.path.join(
                full_sim_file_name
                + '.truth_hits' + '.h5'
                ),
            events.truth_hits)

def read_events_file(full_file_name):
    """
    Loads events from disk

    full_file_stub includes path but not extension

    11/12/21 TS
    """

    import deepdish as dd
    import pickle

    #   Load meta data
    with open(full_file_name + '.meta' + '.pickle', 'rb') as f:
        meta = pickle.load(f)

    #   Load events.truth
    truth = dd.io.load(full_file_name + '.truth' + '.h5')

    #   Load events.truth_hits
    truth_hits = dd.io.load(full_file_name + '.truth_hits' + '.h5')

    return meta, truth, truth_hits

def write_evta_file(events, version='100'):
    """
    Writes events in events['measured_hits'] to an .evta file
    """

    import numpy as np

    #   Get evta file name
    events.meta['file_names'] = add_evta_file_names(
        events.meta['file_names'],
        events
        )

    #   Open file, write header
    f = open(events.meta['file_names']['path_evta'] + '.evta', 'w')

    f.write('Version ' + version + '\n')
    f.write('Type EVTA\n')
    f.write('\n')

    if version=='100':

        for ne in range(len(events.measured_hits['num_hits'])):

            f.write('SE\n')
            f.write(f'ID {events.truth["triggered_id"][ne]:1.0f}\n')
            f.write(f'TI {events.measured_hits["time"][ne]:20.9f}\n')

            for nh in range(events.measured_hits['num_hits'][ne]):
                z = (
                    events.measured_hits['r'][2,nh,ne]
                    + events.meta['params'].cells['geomaga_reference_zo']
                    )
                f.write(
                    'HT 5;'
                    + f'{events.measured_hits["r"][0,nh,ne]*100:10.7f};'
                    + f'{events.measured_hits["r"][1,nh,ne]*100:10.7f};'
                    + f'{z*100:10.7f};'
                    + f'{events.measured_hits["energy"][nh,ne]:10.7f}\n'
                    )

    elif version=='200':

        #   THIS IGNORES RECOMBINATION FLUCTUATIONS FOR MULTIPLE CELL HITS
        #   AND IN ANY CASE THIS CALCULAITON SHOULD HAPPEN
        #   IN response_tools

        #   Start with convenient varibles
        sigmaxy = events.meta['params'] \
            .constants['spatial_resolution']['sigma_xy']
        sigmaz = events.meta['params'] \
            .constants['spatial_resolution']['sigma_z']
        sigma_o = events.meta['params'] \
            .constants['simple_energy_resolution']['sigma_o']
        energy_o = events.meta['params'] \
            .constants['simple_energy_resolution']['energy_o']
        sigma_f = events.meta['params'] \
            .constants['simple_energy_resolution']['sigma_f']
        sigma_energy = np.sqrt(
            (sigma_o * energy_o)**2
            * events.measured_hits['energy']
            / energy_o
            + sigma_f**2
            )

        for ne in range(len(events.measured_hits['num_hits'])):

            f.write('SE\n')
            f.write(f'ID {events.truth["triggered_id"][ne]:1.0f}\n')
            f.write(f'TI {events.measured_hits["time"][ne]:20.9f}\n')

            for nh in range(events.measured_hits['num_hits'][ne]):
                z = (events.measured_hits['r'][2,nh,ne]
                     + events.meta['params']. \
                         cells['geomaga_reference_zo'])
                f.write(
                    'HT 5; '
                    + f'{events.measured_hits["r"][0,nh,ne]*100:10.7f}; '
                    + f'{events.measured_hits["r"][1,nh,ne]*100:10.7f}; '
                    + f'{z*100:10.7f}; '
                    + f'{events.measured_hits["energy"][nh,ne]:10.7f}; '
                    + f'{sigmaxy*100:10.7f}; '
                    + f'{sigmaxy*100:10.7f}; '
                    + f'{sigmaz*100:10.7f}; '
                    + f'{sigma_energy[nh,ne]:10.7f}\n'
                    )

def fix_sim_file_ht_lines(full_sim_file_name_in,
                          full_geo_file_name,
                          full_sim_file_name_out
                          ):

    """
    Fixes problem: ht lines in .sim file have positions based
    on a sub-detecor coordinate system, but should be global coordinates

    This routine not well checked.

    BIG KNOWN PROBLEM: LAST FEW LINES OF FILE ARE NOT CORRECTLY HANDLED

    @author: tshutt
    """
    import numpy as np
    import pickle

    from . import response_tools
    from . import params_tools

    #   Load geo_params that were generated for Cosima, then
    #   with these generate default response params
    with open(full_geo_file_name + '.pickle', 'rb') as f:
        geo_params = pickle.load(f)
    params = params_tools.ResponseParams(geo_params=geo_params)

    #   Open file_names
    f_in = open(full_sim_file_name_in + '.sim')
    f_out = open(full_sim_file_name_out + '.sim', 'w')

    #   Skip to first event - text_line = 'SE'
    text_line = f_in.readline().strip('\n')
    f_out.write(text_line + '\n')
    while text_line=="" or text_line[0:2]!='SE':
        text_line = f_in.readline().strip('\n')
        if text_line[0:2]!='SE':
            print(text_line)
        f_out.write(text_line + '\n')

    r = np.zeros(3)

    #   Read to end of data, denoted by a blank line
    nl = 0
    while text_line!="":

        nl += 1

        if (nl % 1000)==0:
            print('Event ' + str(nl))

        #   Read to HT line
        while text_line[0:2]!='HT':
            text_line = f_in.readline().strip('\n')
            if text_line[0:2]!='HT':
                f_out.write(text_line + '\n')

        #   Read all HT lines, use response_tools cell <-> global
        #   coordinate translation.  For this need cell of each hit. Also,
        #   that tool is in mks, so change r units
        while text_line[0:2]=='HT':

            splitline = text_line[6:].split(';')

            cell = np.array([int(splitline[1].split('_')[1])])

            r[0] = float(splitline[2]) / 100
            r[1] = float(splitline[3]) / 100
            r[2] = float(splitline[4]) / 100

            r_global =response_tools.global_to_cell_coordinates(
                r,
                cell,
                params,
                reverse = True
                ) * 100

            #   reconstruct text line and write out
            text_line_out = text_line[0:6] + ''.join([
                ';'.join(splitline[0:2]),
                ';', f'{r_global[0]:10.5f}',
                ';', f'{r_global[1]:10.5f}',
                ';', f'{r_global[2]:10.5f}',
                ';',
                ';'.join(splitline[5:])
                ])
            f_out.write(text_line_out + '\n')

            #   Read next text_line
            text_line = f_in.readline().strip('\n')

        f_out.write(text_line + '\n')

        #   Blank line is end of events
        if text_line=='':
            f_out.write(text_line + '\n')
            break

    #   End lines - last line starts with "TS"
    while 1:
        text_line = f_in.readline().strip('\n')
        f_out.write(text_line + '\n')
        if text_line[0:2]=='TS':
            break

    #   Apparently file open/close here is bad practice
    f_in.close()
    f_out.close()

def get_geo_tag(topology_id, values_id):
    """ Standard geometry tag for file names """
    geo_tag = f'_GeoT{topology_id:02d}v{values_id:02d}'
    return geo_tag

def get_theta_tag(cos_theta):
    """ Standard cos_theta tag for .sim and related file names """
    theta_tag = f'_Cos{cos_theta:2.1f}'
    return theta_tag

def get_sim_file_name(beam,
                       gamma_energy,
                       topology_id,
                       values_id=0,
                       cos_theta=None,
                       ):
    """ Base sim file name, without inc or id #s or extension """

    import sys

    #   Theta tag, if point source.
    if beam=='FarFieldPointSource':
        theta_tag = get_theta_tag(cos_theta)
    elif beam=='FarFieldIsotropic':
        theta_tag = ''
    else:
        sys.exit('Error: unsupported beam type in get_sim_file_name')

    #   Geo tag
    geo_tag = get_geo_tag(topology_id, values_id)

    #   Base and geo names
    file_name = (
        beam
        + f'_{gamma_energy/1000:5.3f}MeV'
        + theta_tag
        + geo_tag
        )

    return file_name

def get_geo_file_name(topology_id, values_id):
    """ Base geometry file name, without extension """

    #   Geo tag
    geo_tag = get_geo_tag(topology_id, values_id)

    #   Base and geo names
    file_name = 'GammaTPC' + geo_tag + '.geo'

    return file_name

def write_geo_files(path, params, values_id=0):
    """
    Writes geometry files - .setup, and copy of params into
        folder set by path

    values_id - integer tags for geometry values constants

    returns file_names
    """

    import sys
    import os
    import pickle

    #   Detector geometry must be 'geomega'
    if params.detector_geometry!='geomega':
        sys.exit('Error in write_geo_files: detector_geometry not geomega')

    #   File names
    file_name = get_geo_file_name(params.topology_id, values_id)

    #   Calculate params to update file image
    params.calculate()

    #   Write .setup from file image in params
    with open(os.path.join(path, file_name + '.setup'), 'w') as f:
        for line in params.setup_file_lines:
            f.write(line)

    #   Write params
    with open(os.path.join(path, file_name + '.pickle'), 'wb') as f:
            pickle.dump(params, f)

    return file_name

def write_source_file(data_path,
                      geo_full_file_name,
                      num_triggers,
                      beam,
                      energy,
                      cos_theta=None,
                      ):

    import os
    from math import acos

    topology_id, values_id = geo_full_file_name.split('_')[1].split('.')[0] \
        .split('T')[1].split('v')
    topology_id = int(topology_id)
    values_id = int(values_id)

    sim_file_name \
        = get_sim_file_name(
            beam,
            energy,
            topology_id,
            values_id,
            cos_theta
            )

    if beam=='FarFieldPointSource':
        theta = acos(cos_theta)
        beam_tag = 'FFPS'
        beam_values = f' {theta:7.5f}   0'
    elif beam=='FarFieldIsotropic':
        beam_tag = 'FFI'
        beam_values = ' '

    lines = ['']
    lines.append('Version          1 \n')
    lines.append('Geometry         ' + geo_full_file_name + '.setup\n')
    lines.append('CheckForOverlaps 1000 0.01 \n')
    lines.append('PhysicsListEM    Livermore \n')
    lines.append('\n')
    lines.append('StoreCalibrate                 true\n')
    lines.append('StoreSimulationInfo            true\n')
    lines.append('StoreOnlyEventsWithEnergyLoss  true  '
                 +'// Only relevant if no trigger criteria is given!\n')
    lines.append('DiscretizeHits                 true\n')
    lines.append('PreTriggerMode                 everyeventwithhits\n')
    lines.append('\n')
    lines.append('Run ' + beam_tag + '\n')
    lines.append(beam_tag + '.FileName           ' + sim_file_name + '\n')
    lines.append(beam_tag + '.NTriggers          '
                 + f'{num_triggers:5.0f} ' + '\n')
    lines.append('\n')
    lines.append('\n')
    lines.append(beam_tag + '.Source One \n')
    lines.append('One.ParticleType        1 \n')
    lines.append('One.Beam                ' + beam + beam_values + '\n')
    lines.append('One.Spectrum            Mono  ' + f'{energy:5.1f}' + '\n')
    lines.append('One.Flux                1000.0')

    #   Write .source file
    with open(os.path.join(data_path, sim_file_name + '.source'), 'w') as f:
        for line in lines:
            f.write(line)

def add_evta_file_names(file_names, events):
    """ First attempt at unified handling for file names.  Work
        in progress. Also returns meta information
    """

    import os

    #   Study tag, if study present
    study_tag = ''
    if 'study' in events.meta['params'].meta:
        # study_tag = events.meta['params']. \
        #     meta['study'].labels['study_tag']
        case = events.meta['params'].meta['case']
        case_tag = events.meta['params']. \
            meta['study'].labels['case_tag'][case]
        study_tag = '.s' + case_tag

    start, end = file_names['base'].split('.inc')

    file_names['evta'] = start + study_tag + '.inc' + end

    file_names['path_evta'] = os.path.join(
        file_names['paths']['data'],
        file_names['evta']
        )

    return file_names



