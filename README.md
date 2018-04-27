# toolbelt
USGS Design Values tool
- Calls USGS API to retrieve design values (https://earthquake.usgs.gov/ws/designmaps/)
- Inputs are in the form of the csv file titled "inputs" and requried input fields for each site to run
- Required fields are: lat, lng, site_class, risk_category, and reference_doc
- Run matlab script to call USGS api for each site in the inputs csv and save results to csv titled "outputs"
