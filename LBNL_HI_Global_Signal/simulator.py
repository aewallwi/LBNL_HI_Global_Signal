import numpy as np
import itertools
import scipy.special as sp
import copy
import healpy as hp
import os
import yaml
from pyuvdata import UVData
from hera_sim.visibilities import vis_cpu
from pyuvsim.simsetup import initialize_uvdata_from_params, _complete_uvdata




def compute_global_signal(
    obs_yaml,
    basename,
    output_dir="./",
    nside_sky=256,
    clobber=False,
    include_gsm=True,
    include_gleam=True,
    nsrcs_gleam=10000,
):
    """Compute visibilities for global sky-model with white noise EoR.

    Simulate visibilities at a single time for a Golomb array of antennas located at the HERA site.
    Uses the Global Sky Model (GSM) to compute foregrounds and simulates EoR signal as a white noise
    healpix map. Antenna configuration is saved to

    Parameters
    ----------
    obs_yaml: str
        path to pyuvsim observation yaml specifying array layout etc...
        can be generated with garrays.initialize_telescope_yamls
    basename: str
        basename for outputs.
    output_dir: str, optional
        path to directory to output simulation products
        deault is './'
    nside_sky: int, optional
        healpix nside setting the resolution of the simulated sky.
        default is 256.
    clobber: bool, optional
        Overwrite existing UVData files.
        Default is False. If False, read any existing files and return them
        rather then simulating them.
    output_dir: str, optional
        directory to write simulation config files.
    clobber: bool, optional
        overwrite existing config files.
    include_gsm: bool, optional.
        include desourced gsm in sky model.
        default is True.
    include_gleam: bool, optional.
        include gleam point sources in sky model.
        default is True.
    nsrcs_gleam: int, optional
        number of brightest gleam sources to include in sky model
        default is 10000

    Returns
    -------
    uvd_fg: UVData object
        UVData with visibilites of foreground emission.

    """
    # only perform simulation if clobber is true and fg_file_name does not exist and eor_file_name does not exist:
    # generate GSM cube
    fg_file_name = os.path.join(
        output_dir,
        basename
        + f"compressed_fg_{include_gsm}_gleam_{include_gleam}_nsrc_{nsrcs_gleam}.uvh5",
    )
    if not os.path.exists(fg_file_name) or clobber:
        from . import skymodel

        uvdata, beams, beam_ids = garrays.initialize_uvdata(
            output_dir=output_dir,
            clobber=clobber,
            obs_param_yaml_name=obs_yaml,
        )
        if include_gsm:
            fgcube = skymodel.initialize_gsm(uvdata.freq_array[0], nside_sky=nside_sky, output_dir=output_dir)
        else:
            fgcube = np.zeros((len(uvdata.freq_array[0], hp.nside2npix(nside_sky))))
        if include_gleam:
            fgcube = skymodel.add_gleam(uvdata.freq_array[0], fgcube, nsrcs=nsrcs_gleam)
        fg_simulator = vis_cpu.VisCPU(
            uvdata=uvdata,
            sky_freqs=uvdata.freq_array[0],
            beams=beams,
            beam_ids=beam_ids,
            sky_intensity=fgcube,
        )
        fg_simulator.simulate()
        fg_simulator.uvdata.vis_units = "Jy"
        uvd_fg = fg_simulator.uvdata
        if compress_by_redundancy:
            # compress with quarter wavelength tolerance.
            uvd_fg.compress_by_redundancy(tol=0.25 * 3e8 / uvd_fg.freq_array.max())
        if not include_autos:
            uvd_fg.select(bls=[ap for ap in uvd_fg.get_antpairs() if ap[0] != ap[1]], inplace=True)
        uvd_fg.write_uvh5(fg_file_name, clobber=True)
    else:
        uvd_fg = UVData()
        uvd_fg.read(fg_file_name)
    # only do eor cube if file does not exist.
    return uvd_fg
