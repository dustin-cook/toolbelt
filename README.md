# toolbelt
USGS Design Values tool
- Calls USGS API to retrieve design values (https://earthquake.usgs.gov/ws/designmaps/)
- Inputs are in the form of the csv file titled "inputs" and requried input fields for each site to run
- Required fields are: lat, lng, site_class, risk_category, and reference_doc
- Run matlab script to call USGS api for each site in the inputs csv and save results to csv titled "outputs"

Spectra Tool
- Converts ground motion files from either NGA formatting, 5 column formatting, or single signal formatting, into single signal format in units of g with dt as the first entry
- Calculates spectra of ground motion from 0 to 5 seconds for damping ratio's of .01, .02, .03, and .05
- Saves and plots spectra
