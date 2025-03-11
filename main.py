from v4.EVI import *

# GET EVIs FOR NATURA 2000 SITE ########################################################################################
# TODO: From the data, it seems that the site DE6308301 has 11 polygons. The function get_EVI() uses get_satellite() to get the EVI for these 11 polygons. Points that I am not sure about:
# In my case, I get exactly the same EVI for these 11 polygons (so maybe they are all the same?)
# In the function get_EVI(), I did not find where these 11 EVIs then get “averaged” to obtain 1 EVI for the site DE6308301
# I am not sure I understand evi_pilot_sites = set(evi_df[(evi_df["Start year"] == self.start) & (evi_df["End year"] == self.end)][self.site_code_variable].unique())  (current line 352). First self.start  is a dictionary and not a date, so maybe you wanted to write self.start[self.site_codes]  or similar? But also, this would work if site_codes is 1, and for more site codes you’d probably need a loop. Second, self.end  may be np.nan , and when it is this expression will always give you false. As this expression is right now, in my example gives empty when it should give something, and the same site (DE6308301) gets re-computed every time I run and re-added to my output/natura2000/evi_natura2000.csv.
# Bottom line: since you’ll change the code next week, I’ll let you do that, and when you are done you should test it on /v3/test/test.py. As you do, make sure that I have put the correct directories at the beginning (i.e., that shp_path_natura, gdf_path_natura, image_path_natura, site_path_natura, param_evi_natura and evi_path_natura are all correct). I think the goal here is to get the EVI for DE6308301 once, and make sure that it does not recompute it a second time. On my side, I will generate the inputs I need to write the code for the matching and write it on Monday.
# TODO: Edo has re-made requirements based on mine, and added scikit-learn without version for now, and pushed it. I should pull and see if on my computer it allows me to install scikit-learn.
### INPUTS #############################################################################################################

# Original shape file
shp_path_natura = r'../data/natura2000/eea_v_3035_100_k_natura2000_p_2022_v01_r00/SHP'

# Treated shape file
gdf_path_natura = r'../data/natura2000/shp'

# Read the existing GeoDataFrame
gdf = gpd.read_file(gdf_path_natura)

# Folder for satellite images
image_path_natura = r'ouput\natura2000\images'

# File with site identifiers
site_path_natura = r'../data/natura2000/sites/natura2000sites_forest_treated.csv'
site_file_natura = pd.read_csv(site_path_natura)[["SITECODE", "DATE_SAC"]].drop_duplicates()

# Ensure DATE_SAC is converted to datetime
site_file_natura['DATE_SAC'] = pd.to_datetime(site_file_natura['DATE_SAC'], errors='coerce')

# Set index and extract year
site_dict = site_file_natura.set_index('SITECODE')['DATE_SAC'].dt.year.to_dict()

# Convert the values from float to int
site_dict = {k: int(v) for k, v in site_dict.items()}

# Parameters for EVI
param_evi_natura = {
    "project": "ee-e4s-biodiversity-loss",
    "crs": 'EPSG:4326',
    "satellite": 'LANDSAT/LE07/C02/T1',
    "scale": 30,  # Landsat 7 has a resolution of 30m
    "start": site_dict,  # Landsat 7 data goes from Dec 1999
    "end": np.nan,  # Landsat 7 data goes until Jan 2024
    "variable": "SITECODE"
}

summer_natura = True  # if true only takes images from May 1 to Sept 30, else all year considered
image_natura = False

# Existing EVI file
evi_path_natura = r'../ouput/natura2000/evi_natura2000.csv'


### PROCESS ############################################################################################################

# Read the existing EVI file to get already processed sites, if it exists
try:
    evi_data = pd.read_csv(evi_path_natura)
    processed_sites = evi_data[(evi_data["End year"] == param_evi_natura.get("end")) | (evi_data["End year"] == param_evi_natura.get("start"))]["SITECODE"].drop_duplicates().values
except FileNotFoundError:
    # If evi_natura2000.csv does not exist, consider no sites processed
    processed_sites = []

# Add sites that do not have images
processed_sites = np.append(processed_sites, ["NL9801023", "ES1110016", "FI0100005"])

# Filter out sites that are already processed
unprocessed_sites = site_file_natura[(~site_file_natura["SITECODE"].isin(processed_sites)) & (site_file_natura["SITECODE"].str.startswith("SE"))].drop_duplicates(subset="SITECODE")

# Define the chunk size
chunk_size = 100  # Adjust as needed

# Process chunks of unprocessed_sites
if not unprocessed_sites.empty:
    for i, chunk in enumerate(np.array_split(unprocessed_sites, max(1, len(unprocessed_sites) // chunk_size))):
        if chunk.empty:  # Skip empty chunks
            continue
        try:
            # Process the current chunk
            evi = EVI(chunk, shp_path_natura, gdf_path_natura, image_path_natura, evi_path_natura,
                      param_evi_natura, summer_natura, image_natura)
            evi.get_evi()
            print(f"Chunk {i + 1} processed successfully.")
        except Exception as e:
            # Log errors and optionally save problematic chunks
            print(f"Error processing chunk {i + 1}: {e}")
else:
    print("No unprocessed sites to process.")

# GET EVIs FOR MINING SITES ############################################################################################

### INPUTS #############################################################################################################

# Original shape file
shp_path_mine = r'data/mining\gdf_treated_wo_ownership.gpkg'

# Treated shape file
gdf_path_mine = r'data/mining\gdf_treated_wo_ownership'

# Folder for satellite images
image_path_mine = r'output\mining\images'

# File with site identifiers
site_path_mine = r'data/mining\gdf_treated_wo_ownership.csv'
site_file_mine = pd.read_csv(site_path_mine)

mine_dict = site_file_mine.set_index('facility_id')['production_start'].to_dict()
# Convert the values from float to int
mine_dict = {k: int(v) for k, v in mine_dict.items()}

# Parameters for EVI
param_evi_mine = {
    "project": "ee-e4s-biodiversity-loss",
    "crs": 'EPSG:4326',
    "satellite": 'LANDSAT/LE07/C02/T1',
    "scale": 30,  # Landsat 7 has a resolution of 30m
    "start": mine_dict,  # Landsat 7 data goes from Dec 1999
    "end": 2023, # Landsat 7 data goes until Jan 2024
    "variable": "facility_id"
}

summer_mine = True  # if true only takes images from May 1 to Sept 30, else all year considered
image_mine = False

# Existing EVI file
evi_path_mine = r'../ouput/mining/evi_mines_wo_ownership.csv'

### PROCESS ############################################################################################################

# Read the existing EVI file to get already processed sites, if it exists
try:
    evi_data_mine = pd.read_csv(evi_path_mine)
    processed_sites = evi_data_mine["facility_id"].drop_duplicates().values
except FileNotFoundError:
    # If evi_mines.csv does not exist, consider no sites processed
    processed_sites = []

# Add sites that have problems
processed_sites = np.append(processed_sites, [])

# Filter out sites that are already processed
unprocessed_sites = site_file_mine[~site_file_mine["facility_id"].isin(processed_sites)]

# Get the EVI of unprocessed sites
if not unprocessed_sites.empty:
    evi = EVI(unprocessed_sites, shp_path_mine, gdf_path_mine, image_path_mine, evi_path_mine,
              param_evi_mine, summer_mine, image_mine)
    evi.get_evi()
