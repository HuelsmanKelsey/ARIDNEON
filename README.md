# ARIDNEON
Programmatically access and open NEON's HSI AOP and TOS vegetation data in drylands.

Step 1: GetSpectra

Step 1.5: Analyze / Map Spectra 

Step 2: GetVeg

Step 2.5 Analyze / Map Veg

Step 3: Combine them

Step 3.5 Analyze them together


Overview
Get Spectra

Get Veg
Choose a site/location and filter the vegetation dataset
Crosswalk / interpret the data
Visualize and understand the distributions 

Whole plot, subplot size-specific analyses of % cover

relative frequency for plot and larger subplots

interpretation caveats and how to present them

NEON in situ (vegetation) data workflow: capturing vegetation “vibe” using PA or % cover:

a PCA that captures the variability within and among domains, sites, plots, subplots and subplots (matrix rows) in terms of cover and composition, which is then used to characterize each organizational level: 
[a diagram of subplot, plot, site, domain, all NEON sites]





I’m envisioning output that shares the loadings of each species in each principle component AND that provides a score in each component of each subplot. 

Ideally, there will be an output that shows the clustering of plots in the component space and demonstrates how those plots are positioned in terms of species loading values in each principle component. 


These values will be merged with the information extracted in the next step from remotely sensed imagery, so consider how they translate to a spatial pattern/grid/map.


Details: Get Spectra
1. shapefile containing the polygons of each subplot, find the NEON site containing those (pretending I don’t know) OR if a shapefile is not provided accepts the 4-letter code of a NEON site (e.g., ‘BLAN’), which is used to download a hyperspectral imagery tile (DP…) using NEON’s API. 
2. Opens the tile once it is downloaded, and extracts the relevant information (reflectance, bands and associated wavelengths, and other crucial metadata) from the hdf5 file for analysis but does not yet output a dataframe. From the extracted information, 
Summarize into two outputs: 1) a 5-band .tif of the entire tile, and 2) a dataframe of summary information from subsections of that tile. 

By default the shapefile is NEON’s subplots, where they conduct vegetation surveys. 
These subplots are 100m2, 10m2, and 1m2. 
From each of those shapefiles, I want to extract an additional 2 pixel buffer around them.

3. The summary information should include: 
a mean reflectance (across all wavelengths) and a CV (standard deviation/mean) of the reflectance among all pixels within the space. 
For the 100m2 subplot I would also like a summary of the CV among groups of 4 pixels (i.e., the 25 2x2 sections in the 10x10m subplot).

4. I want to pair each of these outputs with the outputs from step 5. 
Each plot and subplot should have spectral summary metrics and principle component scores. 
I want to use the reflectance summary metrics to predict principle component scores. 
