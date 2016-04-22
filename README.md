######################################################################
# Polarized synthesis of molecular lines including the Zeeman effect
######################################################################

## Installation
Prior to running the code, two libraries have to be compiled. This accelerates the
computation of the Stokes profiles. You need to make sure that the gfortran is installed
in your system. If not, it should be pretty straightforward to install it.

1. Enter the directory ``dust`` and type
    
        python setup.py build_ext --inplace

   This generates the pydust library that computes the dust opacity.

2. Enter the directory ``formal`` and type
    
        python setup.py build_ext --inplace

   This generates the pydelopar library that computes the routines for radiative
   transfer.

## Running the code

The code is run from the command line for a configuration file ``configuration.ini`` with

    python lte.py configuration.ini

Additionally, I recommend the users to run it from an interactive ``IPython`` session
because they will have access to all calculations in the class ``lte``:

    In [1]: run lte.py configuration.ini
    In [2]: lte.stokesTotal[k][i][:,:] -> array of size 5xNLambda containing the frequency ([0,:]) and Stokes parameters ([1:5,:], I, Q, U, V)
                                        for impact parameter k and region i

## Dependencies
This code depends on the following set of Python packages, that can be easily installed
using ``conda``, ``pip`` or any packaging system:

    - numpy
    - matplotlib
    - configparser

The code has been tested in Python 3.4+, but I don't discard it can work in Python 2.7+.

## Configuration file
All the configuration is controlled from a configuration file that is human-readable.
An example is copied in the following

    # Configuration File

    #***********************************
    # General information
    #***********************************
    [General]
    File with linelist = "data/soFinal.linelist"
    File with energies = "data/soFinal.energy"
    File with model atmosphere = "atmos/model1G.atmos"
    Output file = "figs/model1G_1kms"
    Microturbulence = 1.0
    Mass = 48.065
    Field topology = 'radial'
    Type = 'star'
    Distance = 250.0
    Impact parameters = 0.0, 0.03, 0.06, 0.09, 0.12

    #***********************************
    # Information about regions in GHz
    #***********************************
    [Wavelength regions]
    Number of regions = 1
    Region 0 = 219.93, 219.97, 500

The ordering of the keywords is free, provided those of the ``general`` section are not
mixed with those of the ``Wavelength regions`` ones.

1. File with linelist = "data/soFinal.linelist" -> file that gives the linelist of the potential lines to be included
2. File with energies = "data/soFinal.energy" -> file with the energies, which are used to compute the partition function
3. File with model atmosphere = "atmos/model1G.atmos" -> model atmosphere
4. Output file = "figs/model1G_1kms" -> a pdf figure will be generated at the end
5. Microturbulence = 1.0 -> microturbulence velocity in km/s
6. Mass = 48.065 -> mass of the molecule in AMU
7. Field topology = 'radial' -> 'radial' (field pointing in the radial direction) or 'perpendicular' (always perpendicular to the line of sight)
8. Type = 'star' -> 'star' (there is a star in the center, which defines the boundary condition) or 'cloud' (no central star)
9. Distance = 250.0 -> distance to the source in pc
10. Impact parameters = 0.0, 0.03, 0.06, 0.09, 0.12 -> impact parameters to compute in arcsec

11. Number of regions = 1 -> number of spectral regions to compute. They can contain several lines
12. Region 0 = 219.93, 219.97, 500 -> definition of each region, with the limits in GHz and the number of steps in the spectral direction

## Model atmosphere
The model atmosphere is defined in a file with the following structure. The directory ``atmos`` contains a few
programs (``genModel.py``) to generate parametric model atmospheres in the appropriate format.

    r [cm]      n[cm^-3]     A(mol)   Tk [K]    Tdust[K]     v[km/s]    B[G]
    100
    2.500e+13   3.380e+10   1.000e-06    2200.000    1667.300       5.000       0.010
    2.710e+13   2.877e+10   1.000e-06    2096.000    1614.300       5.000       0.010
    2.938e+13   2.448e+10   1.000e-06    1996.900    1563.000       5.000       0.010
    3.185e+13   2.083e+10   1.000e-06    1902.500    1513.400       5.000       0.010

## Energies
The file with the energies of the molecule has the following format. The directory ``data`` contains a few
programs (``genModel.py``) to read external files and convert them into the appropriate format.

    ind      E [Hz]               E [K]       g
    0    0.000000e+00           0.00000       1
    1    3.000281e+10           1.43991       3
    2    9.293710e+10           4.46027       5
    3    1.922411e+11           9.22611       7

## Linelist
The file with the linelist for the molecule has the following format. The directory ``data`` contains a few
programs (``genModel.py``) to read external linelists and convert them into the appropriate format.
    
    up   low     nu [Hz]          Elow [Hz]       Aul [Hz]      gu        gl     Ju  Jl
    6   5    13043808370.000 316354896983.24384  2.911e-08   -0.12500    1.87500 1.0 1.0 
    2   1    30001544580.000            0.00000  2.361e-07    2.87500    0.00000 1.0 0.0