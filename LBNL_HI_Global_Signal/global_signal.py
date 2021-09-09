import ares

def get_global_signal_ARES():
    sim = ares.simulations.Global21cm()
    sim.run()
    temp = sim.history['dTb']
    return temp

def convert_temp_to_Jy(temp, frequencies):
    """
    Assuming temperature in mK, frequency in MHz
    """
    # SI units
    c = 3e8
    k = 1.38e-23

    f_Hz = frequencies * 1e6
    t_K = temperature/1e3

    # Rayleigh-Jeans
    I = 2*k/c**2 * f_Hz**2 * t_K
    I *= 1e26  # convert to Jy/Sr
    return I

def intensity_healpix(hp_arr, intensity):
    """
    Add global signal to healpix array with shape N_frequencies x N_pixels. Intensity has shape N_frequencies.
    """
    # this assumes that the 21cm signal is == the avg everywhere, should we instead generate a random field? how?
    # ARES doesn't return power spectrum
    return hp_arr + np.tile(intensity, (1, hp_arr.shape[1]))

def compute_visibilites(beam, 21cm_signal)
    """
    Integral of product of beam and 21cm intensity, both have shape N_frequencies * N_pixels
    """
    integrand = beam * 21cm_signal
    # missing a e^2pi*i factor? arxiv:1709.03984
    return np.sum(integrand, axis=1)  # discrete integral
