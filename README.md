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
| Python 3 | Run the package, scripts, and notebooks; the installation example below uses Python 3.11. |
| [OpenMC](https://docs.openmc.org/en/stable/usersguide/install.html) | Neutron and photon transport, nuclear data, and depletion calculations. Both the executable and the Python API are required. |
| [tokamak_radiation_environment](https://github.com/davidepettinari/tokamak_radiation_environment) | Build the tokamak geometry components. |
| NumPy | Numerical operations in the simulation scripts and notebooks. |
| pandas | Process tally results and tabular data. |
| Matplotlib | Plot spectra and analysis results. |
| JupyterLab or Jupyter Notebook | Required only to open and run the notebooks. |
| Git | Clone this repository and install the geometry dependency from GitHub. |
| pip and setuptools | Install the local Python packages; setuptools is the build backend. |

The shutdown dose rate workflows require **OpenMC 0.15.1 or later**, which
provides [`openmc.deplete.d1s`](https://docs.openmc.org/en/stable/releasenotes/0.15.1.html).
OpenMC's own dependencies, including h5py and SciPy, are installed by its package
manager. Nuclear data libraries must be downloaded and configured separately.

## Installation

### 1. Create an environment with OpenMC

With Conda installed, create and activate an environment on Linux or macOS:

```bash
conda create --name structural-materials --channel conda-forge --strict-channel-priority python=3.11 "openmc>=0.15.1" numpy pandas matplotlib pip
conda activate structural-materials
```

For other installation methods, including Docker, follow the
[OpenMC installation guide](https://docs.openmc.org/en/stable/usersguide/install.html).
Use the same Python environment for all the following steps.

### 2. Install the geometry dependency and this package

```bash
python -m pip install "git+https://github.com/davidepettinari/tokamak_radiation_environment.git"
git clone https://github.com/davidepettinari/structural-materials-analysis.git
cd structural-materials-analysis
python -m pip install -e .
```

If you already cloned the repository, run `python -m pip install -e .` from its
root directory. The editable installation makes changes under `src/` available
without reinstalling. For a regular installation, use `python -m pip install .`.

The current `pyproject.toml` does not declare runtime dependencies, so install
the dependencies above before using the package.

To work with the notebooks, also install JupyterLab in this environment:

```bash
python -m pip install jupyterlab
```

### 3. Configure nuclear data

Download OpenMC-compatible HDF5 data from the
[OpenMC nuclear data collection](https://openmc.org/data/).
The models enable photon transport, so the library must include both neutron
and photon interaction data for the selected materials. Set the path to its
`cross_sections.xml` file:

```bash
export OPENMC_CROSS_SECTIONS="/absolute/path/to/nuclear_data/cross_sections.xml"
```

Activation, depletion, and D1S calculations also require a compatible depletion
chain XML file, including decay photon data for D1S:

```bash
export OPENMC_CHAIN_FILE="/absolute/path/to/chain_endfb80_sfr.xml"
```

Replace these example paths with your own. See the
[OpenMC data configuration guide](https://docs.openmc.org/en/stable/usersguide/data.html)
for details.

Some activation and shutdown dose rate scripts and notebooks currently set
`openmc.config['chain_file']` or `openmc.config['cross_sections']` to absolute
paths on the author's machine. Update those assignments to your local paths,
or remove them to use the environment variables above: explicit assignments
override the environment configuration.

### 4. Check the installation

```bash
openmc --version
python -c "import openmc; from openmc.deplete import d1s; import tokamak_radiation_environment; import structural_material_analysis; print('Imports successful')"
```

These checks verify the executable and package imports. They do not run a
transport calculation or validate the nuclear data configuration.

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
The postprocessing notebooks read these files from their sibling `results/`
directory; place new results there or update the notebook paths. The first wall
notebook also reads simulation stdout logs from `results/output/<material>.o`.

To open the notebooks, run `jupyter lab` from the repository root and select the
installed environment as the kernel.

## Development directions

- Extend the current CSG neutronics and activation workflows.
- Develop a DAGMC model.
- Explore CAD conversion and meshing with GEOUNED, StellarMesh, cad2openmc, and
  Attila4MC.

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

## Materials

### First wall materials

- [x] tinbc
- [x] tungsten
- [x] zrnbc
- [x] vmoc
- [x] vtac
- [x] zrhfc
- [x] zrtac
- [x] nbhfc
- [x] hftac

### Vacuum vessel materials

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
- other

## Analysis

| Function | Metric | Goal |
|---------|--------|------|
| Breed tritium | TBR | ≥ 1.1 |
| Capture neutron energy | EM | > 1 |
| Shield magnets | Magnet lifetime | ≥ 5 FPY |
| Activation analysis | Dose rate | < 100 years |
| Neutron damage | DPA | ≥ 5 years |
