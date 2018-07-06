from __future__ import division
from __future__ import print_function


import sys, io, os
import numpy as np
from scipy import interpolate
from skimage.transform import resize
import pandas as pd
import h5py
import copy
import pickle
sys.path.insert(0,'../../eht-imaging') # or wherever you've installed eht-imaging
import ehtim as eh
import warnings
warnings.filterwarnings("ignore") # to suppress matplotlib backend warnings from ehtim



def open_hdf5(file_name):
    # from lia
    file = h5py.File(file_name, 'r')
    name = list(file.keys())[0]
    d = file[str(name)]
    result = d[:]
    file.close()
    return(result)


def shift_x(array, pix):
    '''
    The hotspot and sim movies weren't centered the same,
    so this shifts the center so they are aligned,
    from lia
    '''
    ss = np.shape(array)
    new = np.zeros((ss[0],ss[1]))
    for j in range(ss[0]):
        # i,j = y,x, want to move in x, so keep y the same, so loop in x
        new[:,j+pix] = array[:,j]
    return(new)


def make_lower_res(sim, dimensions, shift=-15):
    '''
    Brings a higher res movie matrix to lower, interpolating if needed
    using skimage.transform.resize
    This takes the center of the higher res image in a clumsy way
    from lia, but modified
    '''
    # dimesions is np.shape() of new sim
    new_sim = np.zeros(dimensions)
    time, xdim, ydim = dimensions
    for i in np.arange(time):
        frame = sim[i,:,:]
        frame_shift = shift_x(frame,shift)
        frame_resize = resize(frame_shift[112:512-112, 112:512-112], (xdim,ydim), mode='constant')
        new_sim[i,:,:] = frame_resize
    return(new_sim)

def interpolate_hotspot(N_frames, N_ave=5, N_movie_frames=1024, image_size = 128):
    """
    requires pixel_dict.pkl to be created & saved in ../data/interm
    """
    N_ave = int(N_ave)
    # original had 1 period in 50 frames
    x = np.linspace(0,50, N_frames*N_ave)
    # only calculate up to the # of frames we need
    x = x[:int(N_movie_frames*N_ave)]
    # we resized the movie to 128x128
    interp_hotspot = np.zeros((len(x), image_size, image_size))
    pixel_arr = np.arange(image_size)
    with open('../data/interim/pixel_dict.pkl', 'rb') as f:
        pixel_dict = pickle.load(f)
    for i in pixel_arr:
        for j in pixel_arr:
            t, c, k = pixel_dict[(i, j)]
            spline = interpolate.BSpline(t, c, k, extrapolate='periodic')
            interp_hotspot[:,i,j] = spline(x)
    # frame averaging        
    interp_hotspot = np.array([np.mean(interp_hotspot[i:i+N_ave, :, :], axis=0) for i in range(0, len(interp_hotspot), N_ave)])
    return(interp_hotspot)

def generate_hotspot_frames(hotspot_period, movie_duration=60, N_ave=5, N_movie_frames=1024, image_size = 128):
    """
    indexes and/or tiles to resize hotspot movie 
    
    hotspot_period    time in seconds for hotspot to complete 1 cycle
    frame_dur         time in seconds for a frame to last
    movie_duration    time in seconds for the entire movie
    """
    frame_rate = N_movie_frames/movie_duration
    print('generating {0:.2f} sec movie at frame rate of {1:.2f} /sec'.format(movie_duration,frame_rate))
    N_hotspot_frames = int(hotspot_period * frame_rate)
    interp_hotspot = interpolate_hotspot(N_hotspot_frames, N_ave, N_movie_frames, image_size)
    n_loops, remainder = divmod(N_movie_frames,N_hotspot_frames)
    if n_loops == 1 and remainder == 0: 
        return(interp_hotspot)
    elif n_loops ==0 and remainder !=0:
        resized_hotspot = interp_hotspot[:int(remainder),:,:]
    else:
        resized_hotspot = np.tile(interp_hotspot, (int(n_loops),1,1))
        if remainder != 0:
            resized_hotspot = np.concatenate([np.tile(interp_hotspot, (int(n_loops),1,1)), interp_hotspot[:int(remainder),:,:]])
    return(resized_hotspot)


def roll_frames(frames, offset):
    """
    Starts movie at a different frame
    like a phase offset

    offset (float from 0.0 to 1.0)   approx. fraction of the cycle you want to start
    """
    if offset > 1.0 or offset < 0.0:
        raise ValueError('offset should be between 0.0 and 1.0')
    start_frame = int(len(frames)*offset)
    rolled_frames = np.roll(frames, start_frame, axis=0)
    return(rolled_frames)

def combine_frames(sim_frames, hotspot_frames, hotspot_ratio):
    if hotspot_ratio > 1.0 or hotspot_ratio < 0.0:
        raise ValueError('Hotspot ratio should be between 0.0 and 1.0')
    combined_frames = (1 - hotspot_ratio)*sim_frames + hotspot_ratio*hotspot_frames
    return(combined_frames)

    
def add_noise(no_noise_obs, noisy_obs, noise_factor=20.0):
    with_noise_obs = copy.deepcopy(no_noise_obs)
    pertinent_labels = ['vis', 'qvis', 'uvis', 'vvis', 'sigma', 'qsigma', 'usigma', 'vsigma'] 
    for i in np.arange(len(no_noise_obs.data)):
        for label in pertinent_labels:
            with_noise_obs.data[i][label] += noise_factor*noisy_obs.data[i][label]
    return(with_noise_obs)

def make_vex_continuous(vex):
    '''
    Takes a .vex file and edits the start times & duration
    to make it a continuous observation
    '''
    test_vex = copy.deepcopy(vex)

    # first, get rid of slew times
    for i in np.arange(len(test_vex.sched)-1):
        start_hr = test_vex.sched[i]['start_hr']
        end_hr = test_vex.sched[i]['start_hr']+test_vex.sched[i]['scan'][0]['scan_sec']/3600.
        # replace the next obs start hour with last obs end hour
        test_vex.sched[i+1]['start_hr'] = end_hr

    for i in np.arange(len(test_vex.sched)):
        if test_vex.sched[i]['source'] != 'SGRA':
            start_hr = test_vex.sched[i]['start_hr']
            scan_sec = test_vex.sched[i]['scan'][0]['scan_sec']
            #data_size = test_vex.sched[i]['scan'][0]['data_size']

            # take laxt obs of SgrA* and replace timestamps with off source one
            last_sgrA = copy.deepcopy(test_vex.sched[i-1])
            last_sgrA['start_hr'] = start_hr
            for j in np.arange(len(last_sgrA['scan'])):
                last_sgrA['scan'][j]['scan_sec'] = scan_sec
                #last_sgrA['scan'][j]['data_size'] = data_size

            # replace off source obs with this new one
            test_vex.sched[i] = last_sgrA

    return(test_vex)

def scale_vex_obstime(continuous_vex, movie_duration):
    """
    only works with continuous vex files, currently
    """
    start_hr = continuous_vex.sched[0]['start_hr']
    scale = movie_duration/np.sum([entry['scan'][0]['scan_sec'] for entry in continuous_vex.sched])
    for i in np.arange(len(continuous_vex.sched)):
        continuous_vex.sched[i]['start_hr'] = start_hr
        for j in np.arange(len(continuous_vex.sched[i]['scan'])):
            new_scan_sec =  continuous_vex.sched[i]['scan'][j]['scan_sec'] * scale
            continuous_vex.sched[i]['scan'][j]['scan_sec'] = new_scan_sec
        start_hr+=new_scan_sec/3600.
    return(continuous_vex)

def remove_telescope(vex, telescope_string):
    test_vex = copy.deepcopy(vex)
    for i in np.arange(len(test_vex.sched)):
        for j in np.arange(len(test_vex.sched[i]['scan'])):
            if test_vex.sched[i]['scan'][j]['site'] == telescope_string:
                arr = np.arange(len(test_vex.sched[i]['scan']))
                new_scan = {k:test_vex.sched[i]['scan'][k] for k in np.delete(arr, j)}
                test_vex.sched[i]['scan'] = new_scan
    return(test_vex)


def scale_frames(mean_flux, movie_frames):
    #time_arr = np.arange(np.shape(combined_movie.frames)[0])
    #frames = np.array([combined_movie.frames[t].reshape((combined_movie.ydim, combined_movie.xdim)) for t in time_arr])
    time_arr = np.arange(np.shape(movie_frames)[0])
    total_flux_per_frame = np.array([np.sum(movie_frames[t,:,:]) for t in time_arr])
    scale = mean_flux/np.mean(total_flux_per_frame)
    movie_frames = scale * np.array(movie_frames)
    return(movie_frames)


def generate_toy_data(vex, export_name='../data/generated/test.uvfits', movie_duration=60, mean_flux=3.5,hotspot_period=15,hotspot_ratio=0.1, hotspot_offset=0.0, noise_factor=0.0, constant_sim = False, save_movie = False):
    '''
    Runs through the walkthrough in a single function, if you need to generate many test cases. 
    Requires '../data/interim/pixel_dict.pkl' and '../data/interim/trimmed_sim.h5' 
    both are generated in 0_data_prep.ipynb

    args:
    vex                vex observation file (see 1_walkthrough.ipynb)
    export_name        the uvfits file name to be saved
    movie_duration     in seconds
    mean_flux          in Jy
    hotspot_period     in seconds
    hotspot_ratio      the ratio of hotspot to quiescent flux; 0.0 gives no hotspot, 1.0 gives only hotspot
    hotspot_offset     float fraction of cycle to offset the phase of the hotspot (see 1_walkthrough.ipynb)
    noise_factor       multiplicative factor for combination of noise movie with hotspot+quiescent movie; 
                       uses ehtim's thermal noise from eht-imaging/arrays/SITES.txt SEFD values
    constant_sim       if True, uses quiescent flux simulation average frame for movie and observation
    save_movie         eg. 'test.mp4' if you want to save the pre-observation movie; defaults to False
    '''

    # if file already exists, rewrite it
    try:
        os.remove(export_name)
    except OSError:
        pass
    
    # set movie parameters
    M_unit = np.sqrt(27)* eh.RADPERUAS
    pixel = 128 * M_unit/512.
    ra = 17.761122472222223
    dec = -28.992189444444445
    framedur = movie_duration/1024. 

    # hotspot
    hotspot_frames = generate_hotspot_frames(hotspot_period, movie_duration)
    if hotspot_offset > 0.0: 
        hotspot_frames = roll_frames(hotspot_frames, hotspot_offset)
    hotspot_movie = eh.movie.Movie(hotspot_frames, framedur=framedur,psize=pixel,ra=ra, dec=dec)
    
    # quiescent flux
    sim_frames = open_hdf5('../data/interim/trimmed_sim.h5')
    if constant_sim:
        averaged_frame = np.mean(sim_frames, axis=0)
        for i in np.arange(len(sim_frames)):
            sim_frames[i] = averaged_frame 
    sim_movie = eh.movie.Movie(sim_frames, framedur=framedur, psize=pixel,ra=ra, dec=dec)
    
    # combined
    combined_frames = combine_frames(sim_frames, hotspot_frames, hotspot_ratio)
    combined_frames = scale_frames(mean_flux, combined_frames)
    combined_movie = eh.movie.Movie(combined_frames, framedur=framedur, psize=pixel,ra=ra, dec=dec)
    obs = combined_movie.observe_vex(vex, 'SGRA', t_int=framedur, synchronize_start=True, 
                              sgrscat=False, add_th_noise=False,
                              opacitycal=True, ampcal=True, phasecal=True, frcal=True, dcal=True,
                              jones=False, inv_jones=False)
    if save_movie:
        combined_movie.export_mp4('{0}.mp4'.format(export_name.split('.uvfits')[0]))

        
    # optional noise
    if noise_factor > 0.0:
        # dampen sim to very low signal, observe with noise params to get a noisy movie
        dampening_factor = 1e-10
        just_noise_movie = copy.deepcopy(sim_movie)
        just_noise_movie.frames = dampening_factor*np.array(sim_movie.frames)
        noisy_obs = just_noise_movie.observe_vex(vex, 'SGRA', t_int=framedur, synchronize_start=True, 
                                      sgrscat=False, add_th_noise=True,
                                      opacitycal=True, ampcal=True, phasecal=True, frcal=True, dcal=True,
                                      jones=False, inv_jones=False)

        obs = add_noise(obs, noisy_obs, noise_factor) 
        
    obs.save_uvfits(export_name) # export to uvfits 