#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 21:08:02 2022

@author: tshutt
"""

class Events:
    """"
    Flattend representation of events in file.  Arrays padded in the
    scatters dimension. Dimensions are [scatters, event] or
    [space, scatters, event]

    When first called, returns truth, truth_hits and meta

    TODO: When relevant, add selection of charge readout.  Follow
        approach in electron_track_tools, including defined
        method for resetting params with different

    TODO: Review, extend various performance routines: pointing,
        pointing error, CKD ordering, etc.  These neihter well
        checked nor reviewed

    Units mks and keV, requiring a conversion from Megalib cm and keV
    """
    def __init__(self,
                 full_sim_file_name,
                 full_geo_file_name,
                 num_events_to_read=1e10,
                 read_sim_file = False,
                 write_events_files = True
                 ):
        """
        Reads events from .h5 file if found, otherwise from
        .sim file and writes h5 files unless directed not to.
        """

        import os
        import pickle

        from . import file_tools
        from . import params_tools

        #   If .h5 files present, read those
        if os.path.isfile(full_sim_file_name + '.truth' + '.h5') \
            and not read_sim_file:
            print('Reading .h5 files')
            self.meta, self.truth, self.truth_hits \
                = file_tools.read_events_file(full_sim_file_name)

            #   Trim arrays to only contain requested # o events
            if num_events_to_read < self.meta['num_events']:
                self = trim_events(self, num_events_to_read)

            elif num_events_to_read > self.meta['num_events']:
                print('*** Warning: saved files have '
                      + str(self.meta['num_events'])
                      + ' of ' + str(num_events_to_read) + ' requested'
                      )

        #   If not, read geo_params from pickle, then
        #   read .sim file, and write .h5 files unless told not to
        elif os.path.isfile(full_sim_file_name + '.sim'):

            print('Reading .sim file')

            #   Load geo_params that were generated for Cosima
            with open(full_geo_file_name + '.pickle', 'rb') as f:
                geo_params = pickle.load(f)

            #   This flag prevents updating geo_params
            geo_params.live = False

            #   Now read and parse .sim file
            self.truth, self.truth_hits, self.meta \
                = file_tools.read_events_from_sim_file(
                    full_sim_file_name,
                    geo_params,
                    num_events_to_read,
                    )

            #   File names and number of events added to meta data
            data_path, self.meta['sim_file_name'] \
                = os.path.split(full_sim_file_name)
            _, self.meta['geo_file_name'] \
                = os.path.split(full_geo_file_name)
            self.meta['num_events'] = self.truth['num_hits'].size
            self.meta['geo_params'] = geo_params

            #   Write .h5 files
            if write_events_files:
                print('Writing .h5 files')
                file_tools.write_events_file(self, full_sim_file_name)

        #   Otherwise file not found - error
        else:
            print('*** ERROR: ' + full_sim_file_name + '.sim not found')

        #   Generate default response params, and assign to events
        self.params = params_tools.ResponseParams(
            geo_params=self.meta['geo_params']
            )

    def reset_params(self, charge_readout='GAMPixG'):
        """ reset params with new charge readout scheme,
        removing measured things.
        NOTE: this here only as placeholder until charge readout
        is implemented for events
        """

        import params_tools

        #   Find fresh params
        self.params = params_tools.ResponseParams(
            charge_readout,
            self.geo_params,
            )

        #   Remove attributes that depend on params
        kill_keys = ['events_measured', 'cones']
        for kill_key in kill_keys:
            if hasattr(self, kill_key):
                delattr(self, kill_key)

    def apply_detector_response(self, params):
        """ Applies response defined by params to truth_hits,
        creating measure_hits """

        from . import response_tools

        #   Apply response
        self = response_tools.apply_detector_response(self, params)

        #   Save parameters used in meta data
        self.meta['params'] = params

    def write_evta_file(self, evta_version='100'):
        """ Writes events structure into evta file """

        from . import file_tools

        file_tools.write_evta_file(self, evta_version)

    def calculate_order(self,
                           num_hits_list,
                           cut=None,
                           use_truth_hits=False
                           ):
        """ Calculates hit order using ckd method only.
            If use_truth_hits, then calcualtes from MC truth
        """

        import numpy as np

        from . import reconstruction_tools

        if cut is None:
            cut = np.ones(self.truth['num_hits'].shape, dtype=bool)
        #   Calculate from truth hits

        #   Calculate from truth hits
        self.truth_hits['order'] = reconstruction_tools.calculate_order(
            self.truth_hits,
            num_hits_list,
            cut=cut
            )

        #   Calculate from measured hits, if present
        if 'measured_hits' in self:
            self.measured_hits['order'] \
                = reconstruction_tools.calculate_order(
                    self.measured_hits,
                    num_hits_list,
                    cut=cut
                    )

    def calculate_pointing(self):
        """
        Computes theta from energy of 1st Compton scatter

        Acts on truth_hits, and, if it exists, measured_hits

        Computes based on true order (1st hit is 1st scatter), and,
        if present, also on measured order

        9/5/22 TS  Overhaul of previous versions
        """

        import reconstruction_tools

        #   Pointing from truth hits
        reconstruction_tools.construct_pointing_truth_order(self.truth_hits)
        if 'order' in self.truth_hits:
            reconstruction_tools.construct_pointing_measured_order(
                self.truth_hits
                )

        #   Pointing from neasured hits
        if hasattr(self, 'measured_hits'):
            reconstruction_tools.construct_pointing_truth_order(
                self.measured_hits
                )
            if 'order' in self.measured_hits:
                reconstruction_tools.construct_pointing_measured_order(
                    self.measured_hits
                    )

    def calculate_pointing_error(self):
        """
        Computes pointing from vertices and response, and compares these to
        truth from vertices.  Uses truth order, and excludes events
        with missing energy.

        cones.
            dtheta.geometry - opening angle between true and pointing.ray,
                the measured ray
            dtheta.energy - difference between pointing.theta, the
                energy determined cone angle, and the true angle
            dtheta.combined - quadrature sum of these
            ray_length - length of measured ray
            theta - energetically measured cone angle
            energy - compoton scatter measured energy

        Note of length of fields: pointing.theta, etc., have length
        sum(pointing.success).  The length of cones fields is smaller,
        because of the additional requirement of no missing energy

        3/15/20     TS
        3/22 - python port, TS
        """

        import numpy as np

        import reconstruction_tools

        #   dot product
        def dot(a, b):
            return  np.sum(a * b, axis=0)

        #   Exclude missing energy event
        cut = ~self.truth['missing_energy']

        #   Measured pointing
        measured_theta = np.zeros(cut.shape)
        measured_ray = np.zeros((3, cut.size))
        measured_energy = np.zeros(cut.shape)
        measured_total_energy = np.zeros(cut.shape)
        measured_pointing \
            = reconstruction_tools.construct_pointing_truth_order(
            self.measured_hits
            )
        measured_theta[measured_pointing['success']] \
            = measured_pointing['theta']
        measured_ray[:, measured_pointing['success']] \
            = measured_pointing['ray']
        measured_energy[measured_pointing['success']] \
            = measured_pointing['energy']['compton']
        measured_total_energy[measured_pointing['success']] \
            = measured_pointing['energy']['total']

        #   True ray and theta
        truth_theta = np.zeros(cut.shape)
        truth_ray = np.zeros((3, cut.size))
        truth_pointing \
            = reconstruction_tools.construct_pointing_truth_order(
                self.truth_hits
                )
        truth_theta[truth_pointing['success']] = truth_pointing['theta']
        truth_ray[:, truth_pointing['success']] = truth_pointing['ray']

        #   Update cut to include successful measure and truth pointing
        cut = cut & truth_pointing['success'] & measured_pointing['success']

        #   Calculate cones
        cones = {}
        cones['dtheta'] = {}

        #   Geometry error is opening angle between true and measured rays,
        #   measured ray is pointing.ray
        cones['dtheta']['geometry'] = np.arccos(
            dot(measured_ray[:, cut], truth_ray[:, cut]) \
                / np.sqrt(dot(measured_ray[:, cut], measured_ray[:, cut]))
                / np.sqrt(dot(truth_ray[:, cut], truth_ray[:, cut]))
                )

        #   Energy error is difference between true and measured theta,
        #   measured is in pointing, while true is in vertices
        cones['dtheta']['energy'] =  np.abs(
            truth_theta[cut] - measured_theta[cut]
            )

        #   Combined cones is quadrature sum
        cones['dtheta']['combined']= np.sqrt(
            cones['dtheta']['geometry']**2 \
                + cones['dtheta']['energy']**2
            )

        cones['ray_length'] = np.sqrt(
            dot(measured_ray[:, cut], measured_ray[:, cut])
            )
        cones['theta'] = measured_theta[cut]
        cones['energy'] = measured_energy[cut]
        cones['total_energy'] = measured_total_energy[cut]
        cones['theta_true'] = truth_theta[cut]

        self.cones = cones

    def calculate_stats(self,
                        cut=None,
                        indices=None,
                        num_bins=None,
                        all_stats=False
                        ):
        """
        Calculate statistics, with options:
            cut : boolean cut to apply. If not supplied, all events used.
            Distributions calculated if indices and num_bins supplied:
                indices : integer index per event based on some
                    digitization (or other binning) done separtely
                num_bins : length of disribution, must span range of indices
                    Note this is needed because we can't reliably infer range
                    of indices from supplied indices
            all_stats : boolean.  Calculates extended statistics list.
                Default is false, will be speed + memory hit for
                distributions

        stats_sums : returned scalars or distributions
        """

        import numpy as np

        #   Calculate stats on subset of data defined by cut
        def calculate_stats_sums_with_cut(self, this_cut):

            stats_sums = {}

            stats_sums['full_energy'] = np.sum(
                ~self.truth['missing_energy'][this_cut]
                )
            stats_sums['clean_entrance_missing_energy'] =  np.sum(
                self.truth['clean_entrance'][this_cut]
                & ((self.truth['passive_energy'][this_cut]>0.001)
                | (self.truth['escaped_energy'][this_cut]>0.001))
                )
            stats_sums['bad_entrance_missing_energy'] \
                =  np.sum(
                    ~self.truth['clean_entrance'][this_cut]
                    &((self.truth['passive_energy'][this_cut]>0.001)
                      |(self.truth['escaped_energy'][this_cut]>0.001))
                    )

            if  hasattr(self, 'order'):
                stats_sums['full_energy_ordered'] = np.sum(
                    ~self.truth['missing_energy'][this_cut]
                    &self.order['tried'][this_cut]
                    &self.order['ordered'][this_cut]
                    )
                stats_sums['full_energy_disordered'] = np.sum(
                    ~self.truth['missing_energy'][this_cut]
                    &self.order['tried'][this_cut]
                    &~self.order['ordered'][this_cut]
                    )
                stats_sums['clean_entrance_missing_energy_ordered']\
                    =  np.sum(
                        self.truth['clean_entrance'][this_cut]
                        &((self.truth['passive_energy'][this_cut]>0.001)
                          |(self.truth['escaped_energy'][this_cut]>0.001))
                        &self.order['tried'][this_cut]
                        &self.order['ordered'][this_cut]
                        )
                stats_sums['clean_entrance_missing_energy_disordered'] \
                    =  np.sum(
                        self.truth['clean_entrance'][this_cut]
                        &((self.truth['passive_energy'][this_cut]>0.001)
                          |(self.truth['escaped_energy'][this_cut]>0.001))
                        &self.order['tried'][this_cut]
                        &~self.order['ordered'][this_cut]
                        )
                stats_sums['bad_entrance_missing_energy_ordered']\
                    =  np.sum(
                        ~self.truth['clean_entrance'][this_cut]
                        &((self.truth['passive_energy'][this_cut]>0.001)
                          |(self.truth['escaped_energy'][this_cut]>0.001))
                        &self.order['tried'][this_cut]
                        &self.order['ordered'][this_cut]
                        )
                stats_sums['bad_entrance_missing_energy_disordered'] \
                    =  np.sum(
                        ~self.truth['clean_entrance'][this_cut]
                        &((self.truth['passive_energy'][this_cut]>0.001)
                          |(self.truth['escaped_energy'][this_cut]>0.001))
                        &self.order['tried'][this_cut]
                        &~self.order['ordered'][this_cut]
                        )

            if all_stats:

                stats_sums['escaped'] = np.sum(
                    self.truth['escaped_energy'][this_cut]>0.001
                    )
                stats_sums['escaped_back'] = np.sum(
                    self.truth['escaped_back_energy'][this_cut]>0.001
                    )
                stats_sums['escaped_through'] = np.sum(
                    self.truth['escaped_through_energy'][this_cut]>0.001
                    )
                stats_sums['passive'] = np.sum(
                    self.truth['passive_energy'][this_cut]>0.001
                    )
                stats_sums['clean_entrance'] = np.sum(
                    self.truth['clean_entrance'][this_cut]
                    )
                stats_sums['clean_entrance_escaped'] = np.sum(
                    self.truth['clean_entrance'][this_cut]
                    & (self.truth['escaped_energy'][this_cut]>0.001)
                    )
                stats_sums['clean_entrance_passive'] = np.sum(
                    self.truth['clean_entrance'][this_cut]
                    & (self.truth['passive_energy'][this_cut]>0.001)
                    )
                stats_sums['clean_entrance_passive_and_escaped'] =  np.sum(
                    self.truth['clean_entrance'][this_cut]
                    & (self.truth['passive_energy'][this_cut]>0.001)
                    & (self.truth['escaped_energy'][this_cut]>0.001)
                    )
                stats_sums['clean_entrance_passive_or_escaped'] =  np.sum(
                    self.truth['clean_entrance'][this_cut]
                    & ((self.truth['passive_energy'][this_cut]>0.001)
                    | (self.truth['escaped_energy'][this_cut]>0.001))
                    )
                stats_sums['clean_entrance_no_escaped_back'] =  np.sum(
                    self.truth['clean_entrance'][this_cut]
                    & (self.truth['escaped_back_energy'][this_cut]<0.001)
                    )

            return stats_sums

        # #   Stats
        # pointing['stats'] = {}

        # pointing['stats']['initial'] = length(at_least_two_scatters)

        # pointing['stats']['success'] = sum(pointing['success)
        # pointing['stats']['success_no_missing_energy'] = sum( ...
        #     pointing['success& ...
        #     ~vertices.missingenergy)

        # pointing['stats']['at_least_two_scatters'] = sum( ...
        #     at_least_two_scatters)
        # pointing['stats']['at_least_two_scatters_no_missing_energy']
        #   = sum( ...
        #     at_least_two_scatters& ...
        #     ~vertices.missingenergy)

        # pointing['stats']['at_least_three_scatters'] = sum( ...
        #     hits['num_hits']>'] = 3)
        # pointing['stats']['at_least_three_scatters_no_missing_energy']
        #   = sum( ...
        #     hits['num_hits']>'] = 3& ...
        #     ~vertices.missingenergy)

        # pointing['stats']['at_least_four_scatters'] = sum( ...
        #     hits['num_hits']>'] = 4)
        # pointing['stats']['at_least_four_scatters_no_missing_energy']
        #   = sum( ...
        #     hits['num_hits']>'] = 4& ...
        #     ~vertices.missingenergy)

        # pointing['stats']['pair_after_compton'] = sum(sum( ...
        #     vertices.interaction(2:end,~vertices.entrancescatter)'] = '] =  ...
        #     vertices.interactioncode.pair))
        # pointing['stats']['pair_after_compton_no_missing_energy']
        #   = sum(sum( ...
        #     vertices.interaction(2:end,~vertices.entrancescatter&~ ...
        #     vertices.missingenergy)'] = '] =  ...
        #     vertices.interactioncode.pair))

        # #   Save stats from vertices
        # pointing['stats']['vertices'] = vertices.stats

        #   stats is dictionary
        self.stats = {}

        #   If no cut provided, is true for all events
        if cut is None:
            cut = np.ones(self.truth['time'].shape, dtype=bool)

        #   Find scalar sums over all events
        sums = calculate_stats_sums_with_cut(self, cut)
        self.stats['scalar_sums'] = {}
        for key in sums.keys():
            self.stats['scalar_sums'][key] = sums[key]
        self.stats['scalar_sums_sum'] = np.sum(cut)

        #   Disributions: calculate stats in num_bins bins with
        #   bin indices of all events provided
        if not indices is None and not num_bins is None:

            #  Find distribution sums, unfortunately looping over bins
            for bin_index in range(num_bins):

                #   Index from np.digitize starts with 1.
                index = bin_index + 1

                #   Get sums for events in this bin
                sums = calculate_stats_sums_with_cut(
                    self,
                    (indices==index) & cut
                     )

                #   For first bin, initialize distribution sums
                if index==1:
                    distribution_sums = {}
                    for key in sums.keys():
                        distribution_sums[key] = []
                    distribution_all_sums = []

                #   Append current bin to distribution_sums
                distribution_all_sums.append(
                    self.truth['time'][(indices==index) & cut].size
                    )
                for key in distribution_sums.keys():
                    if distribution_all_sums[bin_index] > 0:
                        distribution_sums[key].append(
                            sums[key]
                            )
                    else:
                        distribution_sums[key].append(0)

            #   Final results are dictionay of numpy arrays
            self.stats['distribution_sums'] = {}
            for key in distribution_sums.keys():
                self.stats['distribution_sums'][key] = np.array(
                    distribution_sums[key]
                    )
            self.stats['distribution_all_sums'] = distribution_all_sums

def trim_events(events, num_trimmed_events):
    """ Used when creating events instance, trims events structure to
    contain only first num_trimmed_events
        """

    import numpy as np

    cut = np.zeros(events.truth['num_hits'].shape, dtype=bool)
    if num_trimmed_events > len(cut):
        print('ERROR: num trimmed events exceeds number of events')
        return

    cut[0:num_trimmed_events] = True

    max_hits = np.max(events.truth['num_hits'][cut])

    for key in events.truth.keys():
        if events.truth[key].ndim==1:
            events.truth[key] \
                = events.truth[key][cut]
        elif events.truth[key].ndim==2:
            events.truth[key] \
                = events.truth[key][0:max_hits, cut]
        elif events.truth[key].ndim==3:
            events.truth[key] \
                = events.truth[key][:, 0:max_hits, cut]

    for key in events.truth_hits.keys():
        if events.truth_hits[key].ndim==1:
            events.truth_hits[key] \
                = events.truth_hits[key][cut]
        elif events.truth_hits[key].ndim==2:
            events.truth_hits[key] \
                = events.truth_hits[key][0:max_hits, cut]
        elif events.truth_hits[key].ndim==3:
            events.truth_hits[key] \
                = events.truth_hits[key][:, 0:max_hits,  cut]

    events.meta['num_events'] = num_trimmed_events

    return events





