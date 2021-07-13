import numpy as np
import yaml
import os
from pyuvsim.simsetup import _complete_uvdata, initialize_uvdata_from_params


def initialize_telescope_yamls(
    basename,
    beam_path,
    output_dir="./",
    clobber=False,
    df=1000e3,
    nf=200,
    f0=50e6,
    start_time=2459122.5835133335,  # .25108 + 8. / 24.,
    integration_time=100,
    Ntimes=864,
    polarization_array=[-5, -6, -7, -8],
    telescope_location=(39.2543951, -113.3613616, 1908.0),  # lat, lon, and elevation above sealevel for marjum pass.
    antenna_positions=[(0.0, 0.0, 0.0)],
    telescope_name="EIGSEP",
):
    """Initialize observing yaml files for simulation.

    Parameters
    ----------
    antenna_positions: array-like
        Nants x 3 array of ENU antenna positions.
    basename: str
        identifying string for yamls
    output_dir: str, optional
        directory to write simulation config files.
        default is current dir ('./')
    clobber: bool, optional
        overwrite existing config files.
        default is False
    df: float, optional
        frequency channel width (Hz)
        default is 400e3
    nf: integer, optional
        number of frequency channels to simulate
        Default is 200
    f0: float, optional
        minimum frequency to simulate (Hz)
        default is 120e6
    start_time: float, optional
        JD starting time of observations (units of days)
        default is 2459122.58351335 which corresponds roughly to
        4 hours LST at the HERA site in South Africa.
    integration_time: float, optional
        Duration of each integration (units of seconds)
        default is 100 seconds
    Ntimes: int, optional
        number of time samples.
    polarization_array: list, optional
        list of polarizations
        default is [-5, -6, -7, -8] ('xx', 'xy', 'yx', 'yy')
    Returns
    -------
    obs_param_yaml_name: str
        path to obs_param yaml file that can be fed into
        pyuvsim.simsetup.initialize_uvdata_from_params()
    telescope_yaml_name: str
        path to telescope yaml file. Referenced in obs_param_yaml file.
    csv_file: str
        path to csv file containing antenna locations.
    """

    csv_name = os.path.join(output_dir, f"{basename}_antenna_layout.csv")
    telescope_yaml_name = os.path.join(output_dir, f"{basename}_telescope_defaults.yaml")

    telescope_yaml_dict = {
        "beam_paths": {0: beam_path},
        "telescope_location": f"{telescope_location}",
        "telescope_name": f"{telescope_name}",
        "x_orientation": "north",
    }
    obs_param_dict = {
        "freq": {
            "Nfreqs": int(nf),
            "bandwidth": float(nf * df),
            "start_freq": float(f0),
        },
        "telescope": {
            "array_layout": csv_name,
            "telescope_config_name": telescope_yaml_name,
        },
        "time": {
            "Ntimes": Ntimes,
            "duration_days": integration_time * Ntimes / (24 * 3600.0),
            "integration_time": integration_time,
            "start_time": start_time,
        },
        "polarization_array": polarization_array,
    }
    if not os.path.exists(telescope_yaml_name) or clobber:
        with open(telescope_yaml_name, "w") as telescope_yaml_file:
            yaml.safe_dump(telescope_yaml_dict, telescope_yaml_file)
    # write csv file.
    lines = []
    lines.append("Name\tNumber\tBeamID\tE    \tN    \tU\n")
    for i, x in enumerate(antenna_positions):
        lines.append(f"ANT{i}\t{i}\t{i}\t{x[0]:.4f}\t{x[1]:.4f}\t{x[2]:.4f}\n")
    if not os.path.exists(csv_name) or clobber:
        with open(csv_name, "w") as csv_file:
            csv_file.writelines(lines)

    obs_param_yaml_name = os.path.join(output_dir, f"{basename}.yaml")
    if not os.path.exists(obs_param_yaml_name) or clobber:
        with open(obs_param_yaml_name, "w") as obs_param_yaml:
            yaml.safe_dump(obs_param_dict, obs_param_yaml)
    return obs_param_yaml_name, telescope_yaml_name, csv_name


def initialize_uvdata(
    obs_param_yaml_name,
    output_dir="./",
    clobber=False,
    compress_by_redundancy=False,
):
    """Prepare configuration files and UVData to run simulation.

    Parameters
    ----------
    output_dir: str, optional
        directory to write simulation config files.
    clobber: bool, optional
        overwrite existing config files.
    compress_by_redundancy: bool, optional
        If True, compress by redundancy. Makes incompatible with VisCPU for now.

    Returns
    -------
    uvd: UVData object
        blank uvdata file
    beams: UVBeam list
    beam_ids: list
        list of beam ids

    """
    uvdata, beams, beam_ids = initialize_uvdata_from_params(obs_param_yaml_name)

    beam_ids = list(beam_ids.values())
    beams.set_obj_mode()
    _complete_uvdata(uvdata, inplace=True)
    if compress_by_redundancy:
        uvdata.compress_by_redundancy(tol=0.25 * 3e8 / uvdata.freq_array.max(), inplace=True)
    return uvdata, beams, beam_ids
