# toolbelt
USGS Design Values tool
- Calls USGS API to retrieve: [design values](https://earthquake.usgs.gov/ws/designmaps/) and [hazard curves](https://earthquake.usgs.gov/ws/nshmp/).
- Calls are made in groups of sites called "site sets" in the "sites.csv" file, using the "USGS_pull_hazard_values.m" matlab function. Results are saved into matlab datastructure. The function is quite slow as the USGS API calls are rate limited.
- The matlab script "USGS_calc_FREr.m" does additional calculations on the data saved in the .mat file to generate specific risk targeted design parameters and design specturm (e.g., for the Functional Recovery Earthquake: FREr)

Spectra Tool
- Converts ground motion files from either NGA formatting, 5 column formatting, or single signal formatting, into single signal format in units of g with dt as the first entry
- Calculates spectra of ground motion from 0 to 5 seconds for damping ratio's of .01, .02, .03, and .05
- Saves and plots spectra
