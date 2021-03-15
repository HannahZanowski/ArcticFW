# ArcticFW
Code for computing Arctic freshwater storage and fluxes in 7 CMIP6 models
for the historical simulation (1950-2015) and SSP1-2.6 and SSP5-8.5 (2015-2100).
The code is currently set up to work with CMIP6 output on the NCAR Cheyenne
supercomputer and Casper analysis cluster, so all directory names will
need to be changed accordingly, as well as those directories that point
to my own local files.

Paper citation: Zanowski, H., A. Jahn, and M.M. Holland, 2021: Arctic Ocean Freshwater
in CMIP6 Ensembles: Declining Sea Ice, Increasing Ocean Storage and Export. JGR Oceans.

1. Files of the form Arctic_external_fw_fluxes_*.py compute the river 
and P-E fluxes into the Arctic Ocean for each model

2. Files of the form Arctic_fw_gateways_*.py compute the solid and liquid
freshwater fluxes through Bering, Nares, Barrow, Fram straits + the Barents
Sea Opening as well as the solid and liquid storage within the boundary defined
by these passageways for each model

3. Files of the form Arctic_fw_davis_*.py compute the solid and liquid freshwater
fluxes through Davis Strait as well as the solid and liquid storage as defined by
the boundaries of Bering, Davis, Fram straits + the Barents Sea Opening

4. Files of the form Arctic_liq_fw_salt_and_vol_decomp_*.py decompose the liquid
freshwater fluxes into their salinity and volume-driven changes in the Bering, Nares, Barrow, Fram,
and Davis straits + the Barents Sea Opening for each model. These decompositions span
2000-2100 and so cover part of the historical simulation as well as all of the future
simulations for both SSP1-2.6 and SSP5-8.5

5. get_models.py is a function that the above codes use to make getting some
specific model and directory info on Casper easier

6. Time series produced by this code can be found at the NSF Arctic Data Center
