# Iversion code

The inversion code was developed for the FLEX DISC project. Inversion conciders only FLOX low-resolution reflectance (400-900 nm) and can be run in Octave.

# Components

All path are in relation to the main folder SCOPE (`.`)

`invertFlox.m` - main script to be run with MATLAB or Octave (> 8.0)
- reads FLOX data (two files, up and down radiance) from `input/dataset inversion`
    - strict format
        - first column - wl
        - header line - XHH_MM_SS
- reads necessary auxiliary files
    - most of them are SCOPE default files (atmosphere, FLUSPECT SACs)
    - **the key file** `input/inversion_input_data.csv` regulates which parameters to retrieve
- runs codes from `src/inversion`
- computes Jacobian, uncertainty and fAPAR and fAPAR_Chl

