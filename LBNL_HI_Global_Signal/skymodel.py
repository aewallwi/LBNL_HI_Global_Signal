import numpy as np
import healpy as hp
from .data import DATA_PATH
import os
import pickle
from scipy.interpolate import interp1d


def initialize_gsm(
    frequencies,
    nside_sky=256,
    save_cube=False,
    output_dir="./",
    clobber=False,
):
    """Initialize GSM.

    Parameters
    ----------
    frequencies: array-like
        1d-array of frequencies (float)
    nside_sky: int
        nsides of healpix sky-model
    save_cube: bool, optional
        if True, save data to a numpy array to save time.
    output_dir: str, optional

    Returns
    -------
    gsmcube: array-like
        (npix, nfreqs) array of healpix values in Jy / sr.
    """
    gsm_file = os.path.join(
        output_dir,
        f"gsm_cube_f0_{frequencies[0]/1e6:.1f}MHz_nf_{len(frequencies)}_df_{np.mean(np.diff(frequencies/1e3)):.1f}_kHz_nside_{nside_sky}.npz",
    )
    if not os.path.exists(gsm_file) or clobber:
        from pygdsm import GlobalSkyModel

        gsm = GlobalSkyModel(freq_unit="Hz")
        rot = hp.rotator.Rotator(coord=["G", "C"])
        gsmcube = np.zeros((len(frequencies), hp.nside2npix(nside_sky)))
        for fnum, f in enumerate(frequencies):
            mapslice = gsm.generate(f)
            mapslice = hp.ud_grade(mapslice, nside_sky)
            # convert from galactic to celestial
            gsmcube[fnum] = rot.rotate_map_pixel(mapslice)
        # convert gsm cube from K to Jy / Sr. multiplying by 2 k_b / lambda^2 * ([Joules / meter^2 / Jy] =1e26)
        gsmcube = 2 * gsmcube * 1.4e-23 / 1e-26 / (3e8 / frequencies[:, None]) ** 2
        np.savez(gsm_file, map=gsmcube)
    else:
        gsmcube = np.load(gsm_file)["map"]
    return gsmcube


def add_gleam(frequencies, hp_input, nsrcs=10000):
    """Add GLEAM sources to a map via nearest neighbor gridding.

    Parameters
    ----------
    frequencies: array-like
        1d-array of frequencies (float)
    hp_input: array-like
        Nfreqs x Npix healpix array (units of Jy / Sr) to add gleam sources to.

    Returns
    -------
    hp_input: array-like
        hp_input array with gleam sources added in.
    """
    npix = hp_input.shape[1]
    nside = hp.npix2nside(npix)
    pixarea = hp.nside2pixarea(nside)
    theta, phi = hp.pix2ang(nside, range(npix))
    gleam_srcs = np.loadtxt(os.path.join(DATA_PATH, "catalogs/gleam_bright.txt"), skiprows=44)[:, :nsrcs]
    for srcrow in gleam_srcs:
        ra = np.radians(srcrow[0])
        zen = np.pi / 2 - np.radians(srcrow[1])
        f200 = srcrow[-1]
        alpha = srcrow[-2]
        pixel = hp.ang2pix(nside, zen, ra)
        hp_input[:, pixel] += f200 * (frequencies / 200e6) ** alpha / pixarea
    return hp_input


def z_to_nu(redshift):
    """
    Redshift to frequency conversion for 21cm line
    """
    f0 = 1420405751.77  # rest freqeuncy in Hz
    return f0/(z+1)

def temp_to_jy(temperature, frequency, t_unit='mK', f_unit='Hz')
    """
    Convert brightness temperature to specific intensity (Jy/Sr)
    """
    c = 3e8  # speed of light in m/s
    k = 1.38e-23  # boltzmann constant in J/K

    if t_unit == 'mK':
        temperature *= 1e3  # convert to K
    elif t_unit != 'K':
        raise ValueError('Invalid temperature unit, must be mK or K')

    if f_unit == 'MHz':
        frequency *= 1e6  # convert to Hz
    elif f_unit == 'GHz':
        frequency *= 1e9  # convert to Hz
    elif t_unit != 'Hz':
        raise ValueError('Invalid frequency unit, must be Hz, MHz, or GHz')

    I = 2*k/c**2 * frequency**2 * temperature  # intensity in SI
    I *= 1e26  # convert to Jy/Sr
    return I

def add_global_signal(frequencies, hp_input):
    """
    This adds the avg global signal everywhere across the sky
    (i.e. sets the 21cm intensity to the all-sky avg everywhere)
    """
    with open(os.path.join(DATA_PATH, 'global_signal/test_21cm.pkl'), 'r') as f:
        global_sim = pickle.load(f)
    t21 = global_sim.history['dTb']
    redshifts = global_sim.history['z']
    freqs21 = z_to_nu(redshifts)
    # interpolate global signal to the foreground frequencies
    interpolator = interp1d(freqs21, t21, bounds_error=False, fill_value=0)
    t21_interp = interpolator(frequencies)
    intensity = temp_to_Jy(t21_interp, frequencies)
    hp_input += np.tile(intensity, (1, hp_input.shape[1]))
    return hp_input
