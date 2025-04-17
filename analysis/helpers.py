import pandas as pd
import matplotlib.pyplot as plt
import openmc

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

def compute_fw_heating(sp, number_of_neutrons_per_second):
    """
    Compute and return the power deposition on the first wall for different particle types
    from a statepoint object.

    Parameters:
        sp : openmc.StatePoint
            The statepoint object loaded from an OpenMC simulation.
        number_of_neutrons_per_second : float
            The neutron source rate [neutrons/s].

    Returns:
        pd.DataFrame
            A table with particle type and power deposited in MW.
    """
    eV_to_Joules = 1.60218e-19

    tallies = {
        'neutron': 'heating_fw_n',
        'electron': 'heating_fw_electron',
        'positron': 'heating_fw_positron',
        'photon': 'heating_fw_photon',
        'total': 'heating_fw'
    }

    data = []
    for particle, tally_name in tallies.items():
        tally = sp.get_tally(name=tally_name)
        heat_MW = eV_to_Joules * number_of_neutrons_per_second * tally.mean.sum() / 1e6
        data.append((particle.capitalize(), heat_MW))

    df = pd.DataFrame(data, columns=["Particle", "Power [MW]"])
    return df

def plot_neutron_spectrum(tally, fw_volume, number_of_neutrons_per_second):
    """
    Plot the neutron energy spectrum at the first wall from an OpenMC tally.

    Parameters:
        tally (openmc.Tally): Tally object with neutron flux and energy filter.
        fw_volume (float): Volume of the first wall in cm³.
        number_of_neutrons_per_second (float): Source neutron rate in neutrons/second.
    """
    df = tally.get_pandas_dataframe()

    energy_bins = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
    energy_filter = openmc.EnergyFilter(energy_bins)

    openmc_flux_fw = df['mean'].values
    openmc_flux_err_fw = df['std. dev.'].values

    flux = (openmc_flux_fw / fw_volume) * number_of_neutrons_per_second
    flux_err = (openmc_flux_err_fw / fw_volume) * number_of_neutrons_per_second

    plt.figure(figsize=(8, 5))
    plt.step(energy_filter.values[:-1] / 1e6, flux, where='mid', label='Flux')
    plt.fill_between(energy_filter.values[:-1] / 1e6, flux - flux_err, flux + flux_err, alpha=0.3, label='± 1σ')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Energy [MeV]")
    plt.ylabel('Flux [n/cm²-s]')
    plt.title("Neutron spectrum at first wall")
    plt.xlim(2e-7, 30)
    plt.ylim(1e3, 1e14)
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()