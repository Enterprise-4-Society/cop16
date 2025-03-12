import os
from tarfile import data_filter
import pandas as pd
from utils import load_data, get_polygons

###### LOAD DATA ####################################################################################################
# Load the shape files of Natura 2000 sites
# gdf_path = 'data/natura2000/shp'
gdf_path = 'data/natura/eea_v_3035_100_k_natura2000_p_2021_v12_r01/SHP'

# Load Natura 2000 sites data
if os.path.exists(''):
    data_path = ''   # TODO: Add appropriate path name
else:
    data_path = 'data/natura/eea_v_3035_100_k_natura2000_p_2021_v12_r01/TABULAR/MDB/Natura2000_end2021_rev1.mdb'
    # Available here: https://www.eea.europa.eu/en/datahub/datahubitem-view/6fc8ad2d-195d-40f4-bdec-576e7d1268e4?activeAccordion=1091667%2C1084066

# Define variables of interest
dict_variables = {'NATURA2000SITES': ["SITECODE",        # Unique code which forms the key-item within the database.
                                      "SITENAME",        # Site name in the local language.
                                      "COUNTRY_CODE",    # Two digit country code the site belongs to.
                                      # "SITETYPE",        # Type of classification for the site:
                                                          # A: SPAs (Special Protection Areas - sites designated under the Birds Directive);
                                                          # B: SCIs and SACs (Sites of Community Importance and Special Areas of Conservation - sites designated under the Habitats Directive);
                                                          # C: where SPAs and SCIs/SACs boundaries are identical (sites designated under both directives).
                                      # "DATE_SPA",        # Date site classified as SPA.
                                      # "DATE_PROP_SCI",   # Date site proposed as eligible for identification as a Site of Community Importance (SCI).
                                      # "DATE_CONF_SCI",   # Date site has been confirmed as a Site of Community Importance (SCI).
                                      "DATE_SAC",        # Date site designated as SAC.
                                      "AREAHA",          # Surface area of a site in hectares. Although it is an obligatory field, the value -99 is given for sites for which the areas are unknown.
                                                         # A value of 0 can be correct if the site is a cave or cliff. In this case, the field 2.3 is obligatory.
                                      #"LATITUDE",
                                      #"LONGITUDE"
                                      ],
                  "HABITATS": ["SITECODE",          # Unique code which forms the key-item within the database.
                               "HABITATCODE",       # Code for the habitat type listed in Annex I of Directive 92/43/EEC.
                               "DESCRIPTION",       # Name of the habitat type listed in Annex I of Directive 92/43/EEC.
                               "COVER_HA",          # Area cover in hectares (of this habitat on the site of interest).
                               # "RELSURFACE",        # Area of the site covered by the natural habitat type in relation to the total area covered by that natural habitat type within the national territory.
                               # "REPRESENTATIVITY"   # Degree of representativity of the habitat type on the site.
                        ],
                  "BIOREGION": ["SITECODE",             # Unique code which forms the key-item within the database.
                                "BIOGEOGRAPHICREG",
                                "PERCENTAGE"]
                  }

###### DEFINE PARAMETERS ###########################################################################################

# Define the forest habitats
classification_path = 'data/natura/classification/habitat_classification.xlsx'
leave_type = "Broad-leaved"

###### CLEAN DATA ##################################################################################################




###### FILTER FOR RELEVANT SITES ###################################################################################
### Filter based on habitats #######################################################################################
def filter_by_habitat(df, forest_path=classification_path, leave_type=leave_type):
    """
    Filters sites of df based on habitat classification and returns a list of site codes that contain the specified
    habitat type.

    Parameters:
        df (pd.DataFrame): The main DataFrame containing habitat data with columns incl. "HABITATCODE" and "SITECODE".
        forest_path (str): Path to the habitat classification file.
                           Defaults to `classification_path`.
        leave_type (str): The specific habitat type to filter (e.g., "Broad-leaved", "Needle-leaved", "Mixed).
                          Defaults to `leave_type`.

    Returns: list: A list of unique site codes (`SITECODE`) that contain the specified habitat type.
    """
    habitat_class = load_data(forest_path)

    # Filter habitat by type and extract habitat code
    habitat_code = habitat_class[habitat_class["Leave type"] == leave_type]["HABITATCODE"].unique().tolist()

    # Identify sites with the relevant habitat codes
    sites = df[df["HABITATCODE"].isin(habitat_code)]["SITECODE"].unique().tolist()

    return sites


### Filter based on country ########################################################################################
def filter_by_country(df):
    """
    Filters sites of df based on country and returns a list of site codes that are within countries of interest
    according to the compute_cost.py.

    Parameters:
        df (pd.DataFrame): The main DataFrame containing habitat data with columns incl. "COUNTRY_CODE" and "SITECODE".

    Returns: list: A list of unique site codes (`SITECODE`) that are within countries of interest.
    """
    from compute_cost import get_country
    country = get_country()
    return df[df["COUNTRY_CODE"].isin(country)]["SITECODE"].unique()

### Filter for sites subject to reforestation and conservation measures #############################################
def filter_for_sac(df):
    """
    Filters sites of df that are Special Areas of Conservation and returns a list of site codes.

    Parameters:
        df (pd.DataFrame): The main DataFrame containing habitat data with columns incl. "DATE_SAC" and "SITECODE".

    Returns: list: A list of unique site codes (`SITECODE`) that are within countries of interest.
    """
    return df[df["DATE_SAC"].notna()]["SITECODE"].unique()


def get_natura_sites(dfs=None):

    if dfs is None:
        dfs = load_data(data_path, table_dict=dict_variables)

    # Define the databases for the filtering
    natura = dfs["NATURA2000SITES"]
    habitats = dfs["HABITATS"]

    # Filter sites
    filter_habitat_sites = filter_by_habitat(habitats)
    filter_country_sites = filter_by_country(natura)
    filter_sac_sites = filter_for_sac(natura)

    # Find the intersect sites
    sites_of_interest = list(set(filter_habitat_sites) & set(filter_country_sites) & set(filter_sac_sites))

    return sites_of_interest


###### CHECK DATA ##################################################################################################
def check_for_area(df):
    abnormal_sites = df[df['COVER_HA'] > df["AREAHA"]]
    if not abnormal_sites.empty:
        raise ValueError(f"Some sites have a habitat cover larger than total area in ha. These sites include:"
                         f">\n {abnormal_sites["SITECODE"].unique}")


###### MERGE DATA ##################################################################################################





###### PROCESS DATA ################################################################################################
def process_natura_sites(gdf_path=gdf_path, data_path=data_path, dict_variables=dict_variables):
    print("\nProcessing Natura 2000 sites\n")
    gdf = load_data(gdf_path)
    data = load_data(data_path, table_dict=dict_variables)

    # Define sites of interest
    sites_of_interest = get_natura_sites(dfs=data)

    # Filter data for sites of interest
    data_processed = {}
    for table_name, df in data.items():
        data_processed[table_name] = df[df["SITECODE"].isin(sites_of_interest)]

    # Merge AREAHA to habitats for data check
    data_processed["HABITATS"] = data_processed["HABITATS"].merge(data_processed["NATURA2000SITES"][["SITECODE", "AREAHA"]], how="left", on="SITECODE")

    # Check for inconsistencies
    # check_for_area(data_processed["HABITATS"])  # TODO: Solve inconsistencies (9 sites)

    print("\nNatura 2000 data processed\n")

    # Process gdf
    gdf_processed = get_polygons(gdf, "SITECODE", sites_of_interest)

    print("\nNatura 2000 polygons processed\n")

    return data_processed["NATURA2000SITES"], data_processed["HABITATS"], data_processed["BIOREGION"], gdf_processed

# natura, habitats, bioregion, gdf_processed = process_natura_sites()
# TODO: Check ee_geometry is usable for EVI

