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
