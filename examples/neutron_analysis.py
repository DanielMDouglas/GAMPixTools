import sys
import matplotlib.pyplot as plt
import h5py
from GAMPixTools import electron_track_tools


file = sys.argv[1]
f = h5py.File(file)
depth = 0
m_n = 940.6

eff_evs = []
eff_drifted_evs = []
pos_z_evs = []
true_pos_z_evs = []

E_in_pts = []
true_pts = []
drifted_pts = []
detected_pts = []
eff_pts = []

for traj in f['trajectories']:
    # if this is the priamry
    if traj['traj_id'] == 0:
        ev_n = traj['event_id'] # note: events will be in order
        KE = traj['E_start'] - m_n 
        E_in_pts.append(KE)

        print('----------------------------------------------------')
        print(f'Event {ev_n}')
        
        track = electron_track_tools.Track(file, ev_n,
                                           input_format = 'dumpTree',
                                           origin_shift = [0, 0, -5])

        track.reset_params(charge_readout = 'LArPix')
        track.readout_charge(depth)
    
        
        pos_z = 0
        true = track.raw_track['num_e'].sum()
        true_pts.append(true)
        for z in track.raw_track['r'][2, :]:
            if z > 0:
                pos_z += 1
                print(z)
        if pos_z > 0:
            true_pos_z_evs.append(ev_n)
            pos_z_evs.append(ev_n)
                
        print(f'true: {true}')
        print(f'# of pos z: {pos_z}')
        
        
        pos_z = 0
        drifted = track.pixel_samples['samples_raw'].sum()
        drifted_pts.append(drifted)
        for z in track.pixel_samples['r_raw'][2, :]:
            if z > 0:
                pos_z += 1
                print(z)
        if pos_z > 0:
            pos_z_evs.append(ev_n)
                
        print(f'drifted: {drifted}')
        print(f'# of pos z: {pos_z}')


        pos_z = 0
        detected = track.pixel_samples['samples_triggered'].sum()
        detected_pts.append(detected)
        for z in track.pixel_samples['r_triggered'][2, :]:
            if z > 0:
                pos_z += 1
                print(z)
        if pos_z > 0:
            pos_z_evs.append(ev_n)
                
        print(f'detected: {detected}')
        print(f'# of pos z: {pos_z}')

        
        drift_eff = drifted / true
        eff = detected / true
        eff_pts.append(eff)
        if drift_eff > 1:
            print(f'DRIFT EFF > 1: {drift_eff}')
            eff_drifted_evs.append(ev_n)
        if eff > 1:
            print(f'EFF > 1: {eff}')
            eff_evs.append(ev_n)

            
print('Summary:')
print(f'evs w eff > 1: {eff_evs}')
print(f'evs w drifted eff > 1: {eff_drifted_evs}')
print(f'evs w +z: {pos_z_evs}')
#print(f'evs w true +z: {true_pos_z_evs}')
            
plt.style.use('dark_background')

plt.figure(1)
plt.plot(E_in_pts, true_pts, '.')
plt.title('Ionization Yield vs. Neutron Energy')
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Ionized Electrons (number of electrons)')
plt.savefig('ionized.png')
plt.close()

plt.figure(2)
plt.plot(E_in_pts, drifted_pts, '.')
plt.title('Drifted Electrons vs. Neutron Energy')
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Drifted Electrons (number of electrons)')
plt.savefig('drifted.png')
plt.close()


plt.figure(3)
plt.plot(E_in_pts, detected_pts, '.')
plt.title('Detected Electrons vs. Neutron Energy')
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Detected Electrons (number of electrons)')
plt.savefig('detected.png')
plt.close()

plt.figure(4)
plt.plot(E_in_pts, eff_pts, '.')
plt.title('Detection Efficiency vs. Neutron Energy')
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Efficiency (detected/ionized electrons)')
plt.savefig('eff.png')
plt.close()

