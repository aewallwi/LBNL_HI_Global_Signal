import glob
import argparse
import os
import yaml
import re
from pyuvdata import UVBeam
import numpy as np
import pandas as pd
import scipy

def convert_amp_phase_txt_to_uvbeam(
    beam_folder,
    beam_type="efield",
    telescope_name="EIGSEP",
    feed_name="hera_vivaldi",
    feed_version="v0",
    model_name="hera_phaseii_conductors_only",
    model_version="v0",
    history="",
    feed_pol="x",
    rotate_pol=True,
    x_orientation="north",
    efield_to_power=True,
    convert_to_healpix=True,
    nside=64,
    save=False,
    output_file=None,
    clobber=False,
    interpolation_function='az_za_simple',
    funits='GHz'
):
    """
    Convert list of txt files to uvbeam object that can be parsed by simulator.
    """
    # generate feed yaml
    filename_list = glob.glob(beam_folder + "*.txt")
    # extract frequencies
    if funits == 'GHz':
        re_freq = re.compile("f=0.[0-9]{1,3}")
        frequencies = [float(re_freq.findall(fname)[0].split("=")[-1]) * 1e9 for fname in filename_list]
    elif funits == 'MHz':
        re_freq = re.compile("f=[0-9]{1,3}")
        frequencies = [float(re_freq.findall(fname)[0].split("=")[-1]) * 1e6 for fname in filename_list]

    # sort filenames by frequencies
    #print(filename_list)
    #print(frequencies)
    filename_list = sorted(filename_list, key=lambda x: frequencies[filename_list.index(x)])
    frequencies = sorted(frequencies)
    #print(filename_list)
    #print(frequencies)
    uvb = UVBeam()
    uvb.read_cst_beam(
        filename=filename_list,
        beam_type=beam_type,
        feed_pol=feed_pol,
        rotate_pol=True,
        frequency=frequencies,
        model_name=model_name,
        model_version=model_version,
        telescope_name=telescope_name,
        feed_name=feed_name,
        feed_version=feed_version,
        history=history,
        x_orientation=x_orientation,
    )
    uvb.interpolation_function = interpolation_function
    if beam_type == "efield" and efield_to_power:
        uvb.efield_to_power()
    if convert_to_healpix:
        uvb.to_healpix(nside=nside)
        uvb.interpolation_function='healpix_simple'
    if save:
        uvb.write_beamfits(output_file, clobber=clobber)

    return uvb


def convert_realimag_ffs_to_amp_phase_txt(folder_path):
    """
    Convert Real/Imag in FFS format to Mag/Phase for theta and phi components of E-fields
    Parameters:
    ----------
    """

    folder_path_ffs = glob.glob(os.path.join(folder_path, "ffs"))[0]
    outdir = glob.glob(os.path.join(folder_path, "txt"))[0]
    filenames = sorted(glob.glob(os.path.join(folder_path_ffs, "*ffs")))
    freqs = np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+", filename)[-2]) for filename in filenames]) * 1e-3
    order = np.argsort(freqs)
    freqs = freqs[order]
    filenames = np.array(filenames)[order]
    for fi, filename in enumerate(filenames):
        print("Converting {}".format(filename))
        table = pd.read_table(filename)
        idx = np.where(table == "// >> Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi): ")[0][0]
        # Export E-field from FFS format
        data = table[idx + 1 :].reset_index(drop=True)
        phi_beam = np.zeros(data.size)
        theta_beam = np.zeros(data.size)
        E_theta_real = np.zeros(data.size)
        E_theta_imag = np.zeros(data.size)
        E_phi_real = np.zeros(data.size)
        E_phi_imag = np.zeros(data.size)
        for i in range(data.size):
            phi_beam[i], theta_beam[i], E_theta_real[i], E_theta_imag[i], E_phi_real[i], E_phi_imag[i] = [
                float(val) for val in data.loc[i][0].split()
            ]
        E_theta = E_theta_real + 1j * E_theta_imag
        E_phi = E_phi_real + 1j * E_phi_imag
        wave = scipy.constants.c / (freqs[fi] * 1e6)
        k = 2 * np.pi / wave
        # Convert farfield patterns (V) to farfields (V/m)
        # with a reference distance of 1m which is the default of CST
        r = 1
        E_theta = E_theta * np.exp(-1j * k * r) / r
        E_phi = E_phi * np.exp(-1j * k * r) / r
        # Calculate magnitude of E-field (4th & 6th columns)
        E_theta_mag = np.abs(E_theta)
        E_phi_mag = np.abs(E_phi)
        # Calculate phase of E-field (5th & 7th columns)
        E_theta_phs = np.rad2deg(np.angle(E_theta))
        E_theta_phs = np.where(E_theta_phs < 0, E_theta_phs + 360, E_theta_phs)
        E_phi_phs = np.rad2deg(np.angle(E_phi))
        E_phi_phs = np.where(E_phi_phs < 0, E_phi_phs + 360, E_phi_phs)
        # Calculate the 3rd column
        A_x = E_theta_mag
        A_y = E_phi_mag
        delta = E_theta_phs - E_phi_phs
        A_x_sq = A_x ** 2
        A_y_sq = A_y ** 2
        delta_rad = np.deg2rad(delta)
        semi_major = np.sqrt(
            (A_x_sq + A_y_sq + np.sqrt((A_x_sq - A_y_sq) ** 2 + 4 * A_x_sq * A_y_sq * np.cos(delta_rad) ** 2)) / 2
        )
        E_mag = semi_major
        # Calculate the axial ratio (last column)
        E1 = E_theta
        E2 = E_phi
        axial_ratio = np.sqrt(
            (np.abs(E1) ** 2 + np.abs(E2) ** 2 + np.abs(E1 ** 2 + E2 ** 2))
            / (np.abs(E1) ** 2 + np.abs(E2) ** 2 - np.abs(E1 ** 2 + E2 ** 2))
        )
        with open(os.path.join(outdir, "farfield_(f={:g})_[1]_e_field.txt".format(freqs[fi])), "w") as f:
            f.write(
                "Theta [deg.]  Phi   [deg.]  Abs(E   )[V/m   ]   Abs(Theta)[V/m   ]  Phase(Theta)[deg.]  Abs(Phi  )[V/m   ]  Phase(Phi  )[deg.]  Ax.Ratio[      ]\n"
            )
            f.write(
                "------------------------------------------------------------------------------------------------------------------------------------------------------\n"
            )
            for i in range(theta_beam.size):
                if phi_beam[i] != 360:
                    line = [
                        "{:>7,.3f}        {:>7,.3f}       {:>10,.8e}     {:>10,.8e}     {:>10,.8e}      {:>10,.8e}      {:>10,.8e}       {:>10,.8e}\n".format(
                            theta_beam[i],
                            phi_beam[i],
                            E_mag[i],
                            E_theta_mag[i],
                            E_theta_phs[i],
                            E_phi_mag[i],
                            E_phi_phs[i],
                            axial_ratio[i],
                        )
                    ]
                    f.writelines(line)
                else:
                    continue
