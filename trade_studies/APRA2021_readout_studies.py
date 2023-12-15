class PixelPitchDriftDistance:
    """
      Separately scan:
          - drift distance
          - pixels pitch

      null is optional, default is False.  If True, then uses uses 5 cm drft
      and default pixel pitch - that is, nullifies the study

      As in all 2 parameter scans, the results can be thought of in two
      basis sets: (a) in the linear space of cases, or (b) in the 2d space
      of [drift_distance][pitch]
      Tools to access both basis sets are supplied in study.kit

      10/20 phython port   TS
    """
    def __init__(self, null=False):

        import numpy as np

        #   These are the values to study
        if not null:
            drift_distances = [0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, \
                                0.075, 0.1, 0.2, 0.3]
            pitches = [0.2e-3, 0.3e-3, 0.4e-3, 0.5e-3, 0.75e-3, 1e-3]
            select_drift_distance_indices = [2, 6, 10]
            select_pitch_indices = [0, 2, 5]
            default_drift_distance_index = 6
            default_pitch_index= 2
        else:
            import params_tools
            params = params_tools.Params()
            drift_distances = [0.05]
            pitches = [params.inputs['pixels']['pitch']]
            select_drift_distance_indices = [0]
            select_pitch_indices = [0]
            default_drift_distance_index = 0
            default_pitch_index = 0

        #  Kit contains useful meta information - including tools to switch
        #  between bases of all cases or the two variables
        self.kit = {}
        self.kit['drift_distance'] = []
        self.kit['pitch'] = []
        self.kit['case'] = \
            [[0 for i in range(len(pitches))] \
             for j in range(len(drift_distances))]
        self.kit['num_pitches'] = len(pitches)
        self.kit['num_drift_distances'] = len(drift_distances)
        self.kit['num_cases'] = \
            self.kit['num_pitches'] \
            * self.kit['num_drift_distances']
        self.kit['pitches'] = np.unique(pitches)
        self.kit['drift_distances'] = np.unique(drift_distances)
        self.kit['select_drift_distance_indices'] = \
            select_drift_distance_indices
        self.kit['select_pitch_indices'] = \
            select_pitch_indices
        self.kit['default_drift_distance_index'] = \
            default_drift_distance_index
        self.kit['default_pitch_index'] = default_pitch_index

        #   Initialize labels
        self.labels = {}

        self.labels['study_name'] = 'APRA2021 Pixel Pitch Drift'
        self.labels['study_tag'] = 'a21pix'

        self.labels['case'] = []
        self.labels['case_tag'] = []
        self.labels['drift_distance'] = []
        self.labels['pitch'] = []

        #   Initiali main arrays
        self.fields = []
        self.sub_fields = []
        self.values = []

        #   now assign values. Loop separately over pitch and drift, though
        #   most things are indexed by case
        nc = 0
        for nd in range(len(drift_distances)):
            for npp in range(len(pitches)):

                #   kit saves information both by pithc/drift and
                #   by case
                self.kit['pitch'].append(npp)
                self.kit['drift_distance'].append(nd)
                self.kit['case'][nd][npp] = nc
                nc += 1

                #   labels and tags
                self.labels['case'].append(
                    f'{drift_distances[nd]*100:3.1f} cm drift'
                    + ', {pitches[npp]*1e6:3.0f}'
                    + r' $\mu$' + ' pitch'
                    )
                self.labels['case_tag'].append(
                    'D' + str(nd) + 'P' + str(npp))
                self.labels['drift_distance'].append(
                    f'{drift_distances[nd]*100:3.1f} cm drift'
                    )
                self.labels['pitch'].append(
                    f'{pitches[npp]*1e6:4.0f}'
                    + r' $\mu$m' + ' pitch'
                    )

                #   Save params fields and sub-fields, and values
                self.fields.append(['data_params', 'pixels'])
                self.sub_fields.append(['drift_distance', 'pitch'])
                self.values.append([drift_distances[nd], pitches[npp]])

class SpatialResolution:
    """
      Scan over spatial resolution

      null is optional, default is False.  If True, then uses default

      Note that drift_distance is added to charge_drift, but really does
      not belong there.  This is a bit ugly, but works.

      11/20 python port   TS
    """
    def __init__(self, null=False):

        #   Spatial resolutions to change, and default
        if not null:
            resolution_xy = [0.2e-3, 0.3e-3, 0.4e-3, 0.5e-3, 0.75e-3, 1e-3]
            resolution_z = resolution_xy
        else:
            import response_definition
            response = response_definition.Response()
            resolution_xy = response.spatial_resolution.sigma_xy
            resolution_z = response.spatial_resolution.sigma_z

        #   Build study structure
        self.fields = \
            [['spatial_resolution', 'spatial_resolution']] \
            * len(resolution_xy)
        self.sub_fields =\
            [['sigma_xy', 'sigma_z']] \
            * len(resolution_xy)
        self.values = []
        for nc in range(len(resolution_xy)):
            self.values.append(
                [resolution_xy[nc], resolution_z[nc]]
                )

        #  Kit contains useful meta information
        self.kit = {}
        self.kit['num_cases'] = len(self.values)

        #   Study labels
        self.labels = {}

        self.labels['study_name'] = 'APRA2021 Spatial Resolution'
        self.labels['study_tag'] = 'apra21sr'

        self.labels['case'] = []
        self.labels['case_tag'] = []
        for nc in range(len(resolution_xy)):
            self.labels['case'].append(
                r'$\sigma_{xzy}$:'
                + f' {resolution_xy[nc]*1e6:3.0f}'
                + r' $\um$'
                )
            self.labels['case_tag'].append(
                'sxyz'
                + f'{resolution_xy[nc]*1e6:04.0f}'
                )

class FullResponse:
    """
      Vary several parameters to create optimistic, nominal and pessimistic
      values for both energy and spatial response

      12/9 TS
    """
    def __init__(self, null=False):

        self.fields = []
        self.sub_fields = []
        self.values = []

        self.labels = {}

        self.labels['study_name'] = 'APRA2021 Full Repsponse'
        self.labels['study_tag'] = 'a21fr'

        self.labels['case'] = []
        self.labels['case_tag'] = []

        #   Optimistic
        self.labels['case'].append('Optimistic')
        self.labels['case_tag'].append('Opt')
        self.fields.append([
            'spatial_resolution',
            'spatial_resolution',
            'coarse_grids',
            'light',
            'material'
            ])
        self.sub_fields.append([
            'sigma_xy',
            'sigma_z',
            'excess_noise_factor',
            'collection',
            'sigma_p'
            ])
        self.values.append([
            200e-6,
            200e-6,
            1,
            0.3,
            0.04
            ])

        #   Nominal
        self.labels['case'].append('Nominal')
        self.labels['case_tag'].append('Nom')
        self.fields.append([
            'spatial_resolution',
            'spatial_resolution',
            'coarse_grids',
            'light',
            'material'
            ])
        self.sub_fields.append([
            'sigma_xy',
            'sigma_z',
            'excess_noise_factor',
            'collection',
            'sigma_p'
            ])
        self.values.append([
            400e-6,
            400e-6,
            2,
            0.1,
            0.06
            ])

        #   Pessimistic
        self.labels['case'].append('Pessimistic')
        self.labels['case_tag'].append('Pes')
        self.fields.append([
            'spatial_resolution',
            'spatial_resolution',
            'coarse_grids',
            'light',
            'material'
            ])
        self.sub_fields.append([
            'sigma_xy',
            'sigma_z',
            'excess_noise_factor',
            'collection',
            'sigma_p'
            ])
        self.values.append([
            750e-6,
            750e-6,
            5,
            0.05,
            0.06
            ])

        #  Kit contains useful meta information
        self.kit = {}
        self.kit['num_cases'] = len(self.fields)





