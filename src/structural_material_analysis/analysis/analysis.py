import pandas as pd
import matplotlib.pyplot as plt
import openmc
import re

def extract_tbr_contributions(sp, tally_name):
    """
    Extracts TBR-related contributions from a given tally in an OpenMC statepoint.

    Parameters
    ----------
    sp : openmc.StatePoint
        The loaded statepoint object containing tallies.
    tally_name : str
        The name of the tally to analyze (e.g., "TBR", "TBR blanket", "TBR channel").

    Returns
    -------
    dict
        A dictionary with:
        - 'TBR': total tritium breeding ratio (sum of all mean values)
        - 'TBR std. dev.': sum of standard deviations
        - 'Li6': contribution to TBR from Li6
        - 'Li7': contribution to TBR from Li7
    """
    df = sp.get_tally(name=tally_name).get_pandas_dataframe()
    return {
        "TBR": df['mean'].sum(),
        "TBR std. dev.": df['std. dev.'].sum(),
        "Li6": df[df['nuclide'] == 'Li6']['mean'].sum(),
        "Li7": df[df['nuclide'] == 'Li7']['mean'].sum()
    }

def compute_region_heating(sp, tally_base_name, number_of_neutrons_per_second, volume_cm3):
    """
    Compute and return the power deposition on a specific region for different particle types
    from a statepoint object.

    Parameters:
        sp : openmc.StatePoint
            The statepoint object loaded from an OpenMC simulation.
        tally_base_name : str
            The base name of the tally (e.g., 'heating_fw', 'heating_cc', 'heating_blanket').
        number_of_neutrons_per_second : float
            The neutron source rate [neutrons/s].
        volume_cm3 : float
            The volume of the region in cm³.

    Returns:
        pd.DataFrame
            A table with particle type and power deposited in MW.
    """
    eV_to_Joules = 1.60218e-19
    volume_m3 = volume_cm3 * 1e-6  # Convert cm³ to m³

    particles = ['n', 'electron', 'positron', 'photon', '']  # '' for total
    data = []

    for particle in particles:
        tally_name = f"{tally_base_name}_{particle}" if particle else tally_base_name
        tally = sp.get_tally(name=tally_name)
        power_MW = eV_to_Joules * number_of_neutrons_per_second * tally.mean.sum() / 1e6
        power_density = power_MW / volume_m3  # MW/m³
        label = particle.capitalize() if particle else "Total"
        data.append((label, power_MW, power_density))

    df = pd.DataFrame(data, columns=["Particle", "Power [MW]", "Power Density [MW/m³]"])
    return df    

# The following commented-out code is an alternative way to compute the heating on the first wall.
    # tallies = {
    #     'neutron': 'heating_fw_n',
    #     'electron': 'heating_fw_electron',
    #     'positron': 'heating_fw_positron',
    #     'photon': 'heating_fw_photon',
    #     'total': 'heating_fw'
    # }

    # data = []
    # for particle, tally_name in tallies.items():
    #     tally = sp.get_tally(name=tally_name)
    #     power_MW  = eV_to_Joules * number_of_neutrons_per_second * tally.mean.sum() / 1e6
    #     power_density = power_MW / volume_m3  # MW/m³
    #     data.append((particle.capitalize(), power_MW, power_density))

    # df = pd.DataFrame(data, columns=["Particle", "Power [MW]", "Power Density [MW/m³]"])
    # return df

def plot_neutron_spectrum(tally, region_name, region_volume_cm3, number_of_neutrons_per_second):
    """
    Plot the neutron energy spectrum at a specific region from an OpenMC tally.

    Parameters:
        tally (openmc.Tally): Tally object with neutron flux and energy filter.
        region_name (str): Name of the region for which the spectrum is plotted.
        region_volume_cm3 (float): Volume of the region in cm³.
        number_of_neutrons_per_second (float): Source neutron rate in neutrons/second.
    """
    df = tally.get_pandas_dataframe()

    energy_bins = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
    energy_filter = openmc.EnergyFilter(energy_bins)

    flux_raw = df['mean'].values
    flux_err_raw = df['std. dev.'].values
    
    flux = (flux_raw / region_volume_cm3) * number_of_neutrons_per_second
    flux_err = (flux_err_raw / region_volume_cm3) * number_of_neutrons_per_second

    plt.figure(figsize=(8, 5))
    plt.step(energy_filter.values[:-1] / 1e6, flux, where='mid', label='Flux')
    plt.fill_between(energy_filter.values[:-1] / 1e6, flux - flux_err, flux + flux_err, alpha=0.3, label='± 1σ')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Energy [MeV]")
    plt.ylabel('Flux [n/cm²-s]')
    plt.title(f"Neutron spectrum in {region_name}")
    plt.xlim(2e-7, 30)
    plt.ylim(1e3, 1e14)
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.show()
    
    print(f"Mean neutron flux in {region_name}: {sum(flux):.2e} ± {sum(flux_err):.2e} n/cm²-s")

def compute_dpa_from_tally(df, number_of_neutrons_per_second, atomic_mass_u, density_g_per_cm3, volume_cm3, Ed_eV=55):
    """
    Compute the displacements per atom (DPA) from a damage-energy tally.

    Parameters
    ----------
    df : pandas.DataFrame
        The tally DataFrame from OpenMC containing damage energy in eV.
    number_of_neutrons_per_second : float
        Neutron production rate [n/s].
    atomic_mass_u : float
        Average atomic mass of the material [atomic mass units].
    density_g_per_cm3 : float
        Density of the material [g/cm^3].
    volume_cm3 : float
        Volume of the region considered [cm^3].
    Ed_eV : float, optional
        Displacement threshold energy (default is 55 eV).

    Returns
    -------
    dict
        Dictionary containing intermediate values and final DPA.
    """
    u_to_g = 1.66054e-24  # atomic mass unit in grams
    damage_energy_in_eV = df['mean'].sum()

    displacement_per_source_neutron = damage_energy_in_eV / (2 * Ed_eV)
    displacement_per_source_neutron_recomb = displacement_per_source_neutron * 0.8

    number_of_neutrons_per_year = number_of_neutrons_per_second * 60 * 60 * 24 * 365.25
    displacements_total = number_of_neutrons_per_year * displacement_per_source_neutron_recomb

    atom_mass_g = atomic_mass_u * u_to_g
    number_of_atoms = volume_cm3 * density_g_per_cm3 / atom_mass_g

    DPA = displacements_total / number_of_atoms

    return {
        "damage_energy_in_eV": damage_energy_in_eV,
        "displacement_per_source_neutron": displacement_per_source_neutron,
        "displacement_per_source_neutron_recomb": displacement_per_source_neutron_recomb,
        "number_of_neutrons_per_year": number_of_neutrons_per_year,
        "displacements_total": displacements_total,
        "number_of_atoms": number_of_atoms,
        "DPA": DPA
    }
    
def extract_leakage(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Search for the line containing the Leakage Fraction
    match = re.search(r'Leakage Fraction\s*=\s*([\d.eE+-]+)\s*\+/-\s*([\d.eE+-]+)', content)
    
    if match:
        leakage = float(match.group(1))
        error = float(match.group(2))
        print(f"Leakage Fraction: {leakage:.5e} ± {error:.5e}")
        return leakage, error
    else:
        print("Leakage Fraction not found in the file.")
        return None, None