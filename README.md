# Structural Materials Analysis

OpenMC models and analysis tools for studying first wall and vacuum vessel
materials in an ARC-class fusion reactor. The repository includes constructive
solid geometry (CSG) models, material definitions, neutronics calculations, and
notebooks for activation, transmutation, and shutdown dose rate studies.

**Install the local Python package before running the scripts or notebooks.**
They import `structural_material_analysis` from the `src/` directory.

## Requirements and dependencies

| Dependency | Purpose |
|------------|---------|
| [OpenMC](https://docs.openmc.org/en/stable/usersguide/install.html) | Neutron and photon transport, nuclear data, and depletion calculations. Both the executable and the Python API are required. |
| [tokamak_radiation_environment](https://github.com/SteSeg/tokamak_radiation_environment.git) | Build the tokamak geometry components. |

The shutdown dose rate workflows require **OpenMC 0.15.1 or later**, which
provides [`openmc.deplete.d1s`](https://docs.openmc.org/en/stable/releasenotes/0.15.1.html).

## Installation

### 1. Create an environment with OpenMC

For more detailed instructions on configuring and installing OpenMC, see the [Installation and Configuration](https://docs.openmc.org/en/stable/usersguide/install.html#usersguide-install) section of the official OpenMC User's Manual.


### 2. Install the geometry dependency and this package

```bash
git clone https://github.com/davidepettinari/structural-materials-analysis.git
cd structural-materials-analysis
pip install -e .
```

## Repository structure

```text
src/structural_material_analysis/
    materials/              Material compositions and properties
    geometry/               CSG reactor geometry
    analysis/               Tally processing, spectra, heating, and DPA
analysis/
    firstwall/              First wall models and postprocessing
    vacuumvessel/
        1_neutronics/       Vacuum vessel neutronics
        2_activation_transmutation_/  Activation and transmutation notebooks
        3_shut_down_dose_rate/        D1S shutdown dose rate studies
```

## Usage

After installation and nuclear data configuration, run a model from its own
directory. For example, starting at the repository root:

```bash
cd analysis/firstwall
python openmc_model.py -m tungsten
```

For the vacuum vessel, starting at the repository root:

```bash
cd analysis/vacuumvessel/1_neutronics
python openmc_model.py -m inconel718
```

Use `python openmc_model.py --help` in either directory to view the command-line
options. Review `settings.particles` and `settings.batches` before a simulation:
the first wall model currently uses 100 batches of 100 million particles each.

The scripts write a statepoint named `<material>.h5` in the working directory.
The postprocessing notebooks read these files from their `natural_flibe_results/` or `enr_flibe_results/`
directory; place new results there or update the notebook paths. The first wall
notebook also reads simulation stdout logs from `results/output/<material>.o`.

## Atomic displacement energies used to compute DPA

| Element | Ed (eV) | Element | Ed (eV) |
|---------|---------|---------|---------|
| Be      | 31      | Co      | 40      |
| C       | 31      | Ni      | 40      |
| Mg      | 25      | Cu      | 40      |
| Al      | 27      | Zr      | 40      |
| Si      | 25      | Nb      | 40      |
| Ca      | 40      | Mo      | 60      |
| Ti      | 40      | Ag      | 60      |
| V       | 40      | Ta      | 90      |
| Cr      | 40      | W       | 55      |
| Mn      | 40      | Au      | 30      |
| Fe      | 40      | Pb      | 25      |

The displacement threshold energy, E_d, is 25 eV for all other elements.

[The FISPACT-II User Manual](https://fispact.ukaea.uk/wp-content/uploads/2021/05/FISPACT_manual_old.pdf)


## Vacuum vessel materials

The neutronics script currently accepts `inconel718`, `ss316l`,
`silicon_carbide`, `vcrti`, `ods_eurofer`, `ht9`, and `rusfer`.

The broader list of candidate materials is:

- Inconel 718
- Eurofer
- Hiperfer
- Rusfer
- CLAM
- F82H
- XM19
- SS316
- Steel 660
- ODS alloys
- G91 T1
- G91 T2
- Silicon Carbide
