import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import tqdm
import h5py
import scipy.stats as st

from GAMPixTools import electron_track_tools

# edep-sim generation volume is set so that z < 0, so z == 0 defines the anode
depth = 0 # include no additional shifting

# Find files - get list of files in folder with single energy tracks.
# input_edepsim_file = '/home/dan/studies/GAMPix/muon_sample/dumpTree_single_particle_1c2374f4-aaeb-4004-983f-213dcbc63871.h5'
# input_edepsim_file = '/home/dan/studies/GAMPix/neutrino_sample/joined.genie.edepsim.04de6916-8bdc-427c-8517-123000e9a743.h5'
input_edepsim_dir = '/home/dan/studies/GAMPix/neutrino_sample'

primE = []
intCode = []
trueChargeDep = []
chargeCollGampix = []
chargeCollTiles = []
chargeCollLarpix = []

# evID = 6
# for evID in range(2,10)
for infile in os.listdir(input_edepsim_dir):
    input_edepsim_file = os.path.join(input_edepsim_dir, infile)
    with h5py.File(input_edepsim_file) as f:
        evIDs = np.unique(f['primaries']['eventID'])
    # for evID in tqdm.tqdm(evIDs[:50]):
    # for evID in tqdm.tqdm(range(100)):
    for evID in tqdm.tqdm(evIDs):

        try:
            track = electron_track_tools.Track(input_edepsim_file,
                                               evID,
                                               input_format = 'dumpTree')
        except IndexError:
            continue
            
        track.reset_params(charge_readout='GAMPixD')
        track.readout_charge(depth)

        trackLP = electron_track_tools.Track(input_edepsim_file,
                                             evID,
                                             input_format = 'dumpTree')    
        trackLP.reset_params(charge_readout='LArPix')
        trackLP.readout_charge(depth)

        if np.sum(track.pixel_samples["samples_triggered"]) != 0:
            # print (np.sum(track.raw_track["num_e"]))
            # print (np.sum(track.pixel_samples["samples_triggered"]))
            # print (np.sum(track.coarse_tiles_samples["samples_triggered"]))
            # print (np.sum(trackLP.pixel_samples["samples_triggered"]))
            
            trueChargeDep.append(np.sum(track.raw_track["num_e"]))
            chargeCollGampix.append(np.sum(track.pixel_samples["samples_triggered"]))
            chargeCollTiles.append(np.sum(track.coarse_tiles_samples["samples_triggered"]))
            chargeCollLarpix.append(np.sum(trackLP.pixel_samples["samples_triggered"]))
            
            with h5py.File(input_edepsim_file) as f:
                evMask = f['primaries']['eventID'] == evID
                primaryE = f['primaries']['E'][evMask]
                code = f['primaries']['code'][evMask]
                primE.append(primaryE[0])
                intCode.append(code[0])

trueChargeDep = np.array(trueChargeDep)
chargeCollGampix = np.array(chargeCollGampix)
chargeCollTiles = np.array(chargeCollTiles)
chargeCollLarpix = np.array(chargeCollLarpix)
primE = np.array(primE)
intCode = np.array(intCode)

# fig = plt.figure()
# ax = fig.gca()
# ax.scatter(trueChargeDep, chargeCollGampix)

# -------- scatter plot ---------

fig = plt.figure()
ax = fig.gca()
ax.scatter(trueChargeDep, chargeCollGampix/trueChargeDep)

# Qbins = np.linspace(2.8e7, 3.e7, 11) 
Qbins = np.linspace(0, 1.5e7, 21) 
QbinCenters = 0.5*(Qbins[1:] + Qbins[:-1]) 

# -------- binned stat - dep charge on x-axis ---------

fig = plt.figure()
ax = fig.gca()
eff = chargeCollGampix/trueChargeDep
bs_pix_med = st.binned_statistic(trueChargeDep,
                                 eff,
                                 bins = Qbins,
                                 statistic = 'median')
bs_pix_lq = st.binned_statistic(trueChargeDep,
                                eff,
                                bins = Qbins,
                                statistic = lambda x: np.quantile(x, 0.16))
bs_pix_uq = st.binned_statistic(trueChargeDep,
                                eff,
                                bins = Qbins,
                                statistic = lambda x: np.quantile(x, 0.84))
ax.errorbar(QbinCenters,
            bs_pix_med[0],
            yerr = (bs_pix_med[0] - bs_pix_lq[0],
                    bs_pix_uq[0] - bs_pix_med[0]),
            fmt = 'o',
            label = 'GAMPix Pixels')

eff = chargeCollTiles/trueChargeDep
bs_tile_med = st.binned_statistic(trueChargeDep,
                                  eff,
                                  bins = Qbins,
                                  statistic = 'median')
bs_tile_lq = st.binned_statistic(trueChargeDep,
                                 eff,
                                 bins = Qbins,
                                 statistic = lambda x: np.quantile(x, 0.16))
bs_tile_uq = st.binned_statistic(trueChargeDep,
                                 eff,
                                 bins = Qbins,
                                 statistic = lambda x: np.quantile(x, 0.84))
ax.errorbar(QbinCenters,
            bs_tile_med[0],
            yerr = (bs_tile_med[0] - bs_tile_lq[0],
                    bs_tile_uq[0] - bs_tile_med[0]),
            fmt = 'o',
            label = 'GAMPix Tiles')

eff = chargeCollLarpix/trueChargeDep
bs_lp_med = st.binned_statistic(trueChargeDep,
                                eff,
                                bins = Qbins,
                                statistic = 'median')
bs_lp_lq = st.binned_statistic(trueChargeDep,
                               eff,
                               bins = Qbins,
                               statistic = lambda x: np.quantile(x, 0.16))
bs_lp_uq = st.binned_statistic(trueChargeDep,
                               eff,
                               bins = Qbins,
                               statistic = lambda x: np.quantile(x, 0.84))
ax.errorbar(QbinCenters,
            bs_lp_med[0],
            yerr = (bs_lp_med[0] - bs_lp_lq[0],
                    bs_lp_uq[0] - bs_lp_med[0]),
            fmt = 'o',
            label = 'LArPix-like Pixels')

plt.legend()
ax.set_xlabel(r'$Q_{\mathrm{deposited}}$ [e]')
ax.set_ylabel(r'$Q_{\mathrm{collected}}/Q_{\mathrm{deposited}}$')

Ebins = np.linspace(1, 4, 21) 
EbinCenters = 0.5*(Ebins[1:] + Ebins[:-1]) 

# -------- binned stat - CC - primary E on x-axis ---------

fig = plt.figure()
ax = fig.gca()

CCmask = intCode == b'CC'
NCmask = ~CCmask

eff = chargeCollGampix/trueChargeDep
bs_pix_med = st.binned_statistic(primE[CCmask],
                                 eff[CCmask],
                                 bins = Ebins,
                                 statistic = 'median')
bs_pix_lq = st.binned_statistic(primE[CCmask],
                                eff[CCmask],
                                bins = Ebins,
                                statistic = lambda x: np.quantile(x, 0.16))
bs_pix_uq = st.binned_statistic(primE[CCmask],
                                eff[CCmask],
                                bins = Ebins,
                                statistic = lambda x: np.quantile(x, 0.84))
ax.errorbar(EbinCenters,
            bs_pix_med[0],
            yerr = (bs_pix_med[0] - bs_pix_lq[0],
                    bs_pix_uq[0] - bs_pix_med[0]),
            fmt = 'o',
            label = 'GAMPix Pixels')

eff = chargeCollTiles/trueChargeDep
bs_tile_med = st.binned_statistic(primE[CCmask],
                                  eff[CCmask],
                                  bins = Ebins,
                                  statistic = 'median')
bs_tile_lq = st.binned_statistic(primE[CCmask],
                                 eff[CCmask],
                                 bins = Ebins,
                                 statistic = lambda x: np.quantile(x, 0.16))
bs_tile_uq = st.binned_statistic(primE[CCmask],
                                 eff[CCmask],
                                 bins = Ebins,
                                 statistic = lambda x: np.quantile(x, 0.84))
ax.errorbar(EbinCenters,
            bs_tile_med[0],
            yerr = (bs_tile_med[0] - bs_tile_lq[0],
                    bs_tile_uq[0] - bs_tile_med[0]),
            fmt = 'o',
            label = 'GAMPix Tiles')

eff = chargeCollLarpix/trueChargeDep
bs_lp_med = st.binned_statistic(primE[CCmask],
                                eff[CCmask],
                                bins = Ebins,
                                statistic = 'median')
bs_lp_lq = st.binned_statistic(primE[CCmask],
                               eff[CCmask],
                               bins = Ebins,
                               statistic = lambda x: np.quantile(x, 0.16))
bs_lp_uq = st.binned_statistic(primE[CCmask],
                               eff[CCmask],
                               bins = Ebins,
                               statistic = lambda x: np.quantile(x, 0.84))
ax.errorbar(EbinCenters,
            bs_lp_med[0],
            yerr = (bs_lp_med[0] - bs_lp_lq[0],
                    bs_lp_uq[0] - bs_lp_med[0]),
            fmt = 'o',
            label = 'LArPix-like Pixels')

plt.legend()
plt.title(r'Charged-Current Interactions')
ax.set_xlabel(r'Primary Energy [GeV]')
ax.set_ylabel(r'$Q_{\mathrm{collected}}/Q_{\mathrm{deposited}}$')

# -------- binned stat - NC - primary E on x-axis ---------

fig = plt.figure()
ax = fig.gca()

eff = chargeCollGampix/trueChargeDep
bs_pix_med = st.binned_statistic(primE[NCmask],
                                 eff[NCmask],
                                 bins = Ebins,
                                 statistic = 'median')
bs_pix_lq = st.binned_statistic(primE[NCmask],
                                eff[NCmask],
                                bins = Ebins,
                                statistic = lambda x: np.quantile(x, 0.16))
bs_pix_uq = st.binned_statistic(primE[NCmask],
                                eff[NCmask],
                                bins = Ebins,
                                statistic = lambda x: np.quantile(x, 0.84))
ax.errorbar(EbinCenters,
            bs_pix_med[0],
            yerr = (bs_pix_med[0] - bs_pix_lq[0],
                    bs_pix_uq[0] - bs_pix_med[0]),
            fmt = 'o',
            label = 'GAMPix Pixels')

eff = chargeCollTiles/trueChargeDep
bs_tile_med = st.binned_statistic(primE[NCmask],
                                  eff[NCmask],
                                  bins = Ebins,
                                  statistic = 'median')
bs_tile_lq = st.binned_statistic(primE[NCmask],
                                 eff[NCmask],
                                 bins = Ebins,
                                 statistic = lambda x: np.quantile(x, 0.16))
bs_tile_uq = st.binned_statistic(primE[NCmask],
                                 eff[NCmask],
                                 bins = Ebins,
                                 statistic = lambda x: np.quantile(x, 0.84))
ax.errorbar(EbinCenters,
            bs_tile_med[0],
            yerr = (bs_tile_med[0] - bs_tile_lq[0],
                    bs_tile_uq[0] - bs_tile_med[0]),
            fmt = 'o',
            label = 'GAMPix Tiles')

eff = chargeCollLarpix/trueChargeDep
bs_lp_med = st.binned_statistic(primE[NCmask],
                                eff[NCmask],
                                bins = Ebins,
                                statistic = 'median')
bs_lp_lq = st.binned_statistic(primE[NCmask],
                               eff[NCmask],
                               bins = Ebins,
                               statistic = lambda x: np.quantile(x, 0.16))
bs_lp_uq = st.binned_statistic(primE[NCmask],
                               eff[NCmask],
                               bins = Ebins,
                               statistic = lambda x: np.quantile(x, 0.84))
ax.errorbar(EbinCenters,
            bs_lp_med[0],
            yerr = (bs_lp_med[0] - bs_lp_lq[0],
                    bs_lp_uq[0] - bs_lp_med[0]),
            fmt = 'o',
            label = 'LArPix-like Pixels')

plt.legend()
plt.title(r'Neutral-Current Interactions')
ax.set_xlabel(r'Primary Energy [GeV]')
ax.set_ylabel(r'$Q_{\mathrm{collected}}/Q_{\mathrm{deposited}}$')

# fig = plt.figure()
# ax = fig.gca()
# ax.scatter(trueChargeDep, chargeCollTiles)

fig = plt.figure()
ax = fig.gca()
ax.scatter(trueChargeDep, np.array(chargeCollTiles)/np.array(trueChargeDep))

plt.show()
