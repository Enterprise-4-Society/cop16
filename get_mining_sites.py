import geopandas as gpd
import get_mine_owners
from compute_evi import get_satellite
from utils import load_data

###### LOAD DATA AND DEFINE PARAMETERS ##############################################################################
# Load the Jasanski geodataframe (main shape file)
gdf_path = 'data/mining/open_database_mine_production/data/facilities.gpkg'

# Load the Maus geodataframe (complementary shape file)
gdf2_path = 'data/mining/Maus-etal_2022_V2_allfiles/global_mining_polygons_v2.gpkg'

# Define filter
filter_ownership = False
by_parent_id = True

# Define evi output path
if not filter_ownership:
    evi_path = 'output/mining/no_ownership/'
else:
    evi_path = 'output/mining/with_ownership/'

###### CLEAN DATA ###################################################################################################

### Allocate main-site production start when subsite production start is missing ####################################
def allocate_missing_production_start(gdf):
    """Fills missing production start dates from the parent facility when available."""

    # Step 1: Identify rows without a 'production_start' (NaN values)
    missing_production_start = gdf[gdf['production_start'].isna()]

    # Step 2: Loop through rows without 'production_start' and find matching parent rows with '00' subsite
    for idx, row in missing_production_start.iterrows():

        # Find the parent facility with subsite_facility_id == '00' for this parent_id
        parent_row = gdf[(gdf['parent_facility_id'] == row['parent_facility_id']) & (gdf['subsite_facility_id'] == '00')]

        # If the parent has a 'production_start' value, copy it to the current row
        if not parent_row.empty and not parent_row['production_start'].isna().values[0]:
            # Update the original DataFrame with the parent's 'production_start'
            gdf.at[idx, 'production_start'] = parent_row['production_start'].values[0]

    # Step 3: Remove missing obs and convert date in int
    gdf = gdf[gdf['production_start'].notna()]  # TODO: Decide if we should change the assumption
    gdf["production_start"] = gdf["production_start"].astype(int)

    return gdf

### Allocate main-site geometry when subsites production start are missing ##########################################
def allocate_missing_geometry(gdf):
    """Fills missing geometries from the parent facility when no other subsite has a valid geometry."""
    # Step 1: Identify rows with empty geometry
    missing_geometry = gdf[gdf['geometry'].is_empty]

    # Step 2: Loop through rows with empty geometry and find matching parent rows with subsite_facility_id == '00' that have a valid geometry
    for idx, row in missing_geometry.iterrows():
        parent_id = row['parent_facility_id']

        # Check if there are any subsites (except current) with a valid geometry
        if not gdf[(gdf['parent_facility_id'] == parent_id) &
                   (gdf['subsite_facility_id'] != "00") &
                   (gdf['subsite_facility_id'] != row["subsite_facility_id"]) &
                   (~gdf['geometry'].is_empty)].empty:
            # Skip assigning geometry if another subsite has valid geometry
            continue

        # Find the parent facility with subsite_facility_id == '00' for this parent_id
        parent_rows = gdf[(gdf['parent_facility_id'] == parent_id) & (gdf['subsite_facility_id'] == '00')]

        # If we reach this point, it means the only matching rows for this parent have empty geometry
        # Replace the geometry of the current row with the geometry of the parent row
        if not parent_rows.empty:
            gdf.at[idx, 'geometry'] = parent_rows.iloc[0]['geometry']
    return gdf

def remove_duplicated_sites(gdf, filter_ownership=filter_ownership):
    """Remove duplicated sites from the gdf of interest. Reasons for dropping these sites are indicated in comments."""
    # Manual assignment of duplicates based on associated ownership structure in own_treated
    to_drop =["COM00272.02", "COM00962.02", "COM00583.02", "COM00205.06", "COM01019.03", "COM01016.03", "COM01016.04", "COM00850.02"] + [ # Reason 1: Subsites (0X) or main-sites (00) within one facility have the same geometry and same ownership structure: we keep only one subsite - the lowest number i.e. 01 if 01 and 02 exist or 00 for 00 and 04)
        "COM00824.00", "COM01282.00", "COM01134.00", "COM00969.00", "COM00489.01", "COM00489.02"] + [                                     # Reason 2: Subsites (0X) or main-sites (00) within different facilities have the same geometry and same ownership structure: keep one facility - the one with the lowest value in facility_id)

        "COM00774.01", "COM00774.02", "COM01114.00", "COM00624.02", "COM00662.00"]                                                        # Reason 3: Subsites (0X) or main-sites (00) within different facilities have the same geometry but different ownership structures: drop them altogether

    if not filter_ownership:
        to_drop = (to_drop +
                   ["COM01017.00", "COM01031.00", "COM01188.00","COM00955.01","COM01402.05","COM01402.03","COM01064.02","COM00938.00"])   # Reason 4: Subsites (0X) or main-sites (00) within different facilities have the same geometry: drop them altogether

    gdf = gdf[~gdf["facility_id"].isin(to_drop)]
    return gdf

###### FILTER SITES #################################################################################################
### Filter by production start ######################################################################################
def filter_by_production_start(gdf, year=2001):
    """Keeps only sites with available production start date and after a certain year - generally that matches
    satellite data availability."""
    gdf = gdf[(gdf["production_start"].notna()) & (gdf["production_start"] >= year)]
    return gdf

### Filter by ownership ############################################################################################
def filter_by_ownership(gdf, filter_ownership, by_parent_id):
    """Keeps only sites with ownership information."""
    if filter_ownership:
        own = get_mine_owners.process_mine_ownership()
        id_type = "parent_facility_id" if by_parent_id else "facility_id"
        facility_ids_own = own[id_type].unique()
        gdf = gdf[gdf[id_type].isin(facility_ids_own)]
    return gdf

### Filter for polygon data availability############################################################################

def filter_by_geometry(gdf, gdf2):
    """Keeps only sites with valid polygon data if gdf2 is provided.
    Note: The default behavior of geopandas.sjoin (with predicate="intersects") does not automatically check whether
    all points in a MultiPoint geometry are within a polygon. It checks only based on the entire MultiPoint geometry,
    rather than assessing individual points within it. This means that if any one of the points in a MultiPoint
    geometry intersects a polygon, the join will consider it a match. With "within" as a predicate, we ensure that we
    select only polygons that contains all points within a Multipoint geometry."""
    if gdf2 is not None:
        # Get the id of the facilities whose geometry is completely within a polygon of gdf2
        matched_gdf = gpd.sjoin(gdf, gdf2, how='inner', predicate='within')
        ids = matched_gdf["facility_id"].unique()

        # Add the polygons of the facilities whose geometry is completely within a polygon of gdf2 to the main dataframe
        return gpd.sjoin(gdf2, gdf[gdf["facility_id"].isin(ids)], how='inner', predicate='intersects')
    return gdf

def filter_for_subsites(gdf, filter_ownership):
    """Filter by only keeping sub sites when there is data availability."""
    # Step 1: Identify parent_facility_ids that have more than one unique subsite_facility_id
    # Create a Series counting the number of unique subsite_facility_id per parent_facility_id
    subsite_count_per_parent = gdf.groupby('parent_facility_id')['subsite_facility_id'].nunique()

    # Step 2: Identify parent_facility_ids with only one subsite_facility_id (i.e., only '00')
    multiple_subsite_parents = subsite_count_per_parent[subsite_count_per_parent > 1].index

    # Step 3: Drop rows with subsite_facility_id == '00' if they are in the list of single_subsite_parents
    gdf_filtered = gdf[~((gdf['subsite_facility_id'] == '00') & (gdf['parent_facility_id'].isin(multiple_subsite_parents)))]

    # Reset index after filtering
    if filter_ownership is True:
        return gdf_filtered.reset_index(drop=True).drop_duplicates(subset='geometry', keep='first')
    else:
        return gdf_filtered.reset_index(drop=True)


def filter_mining_sites(gdf, gdf2=None, filter_ownership=True, by_parent_id=True):
    """Applies multiple filtering steps to retain relevant mining sites."""
    gdf = filter_by_production_start(gdf)
    gdf = filter_by_ownership(gdf, filter_ownership, by_parent_id)
    gdf = filter_by_geometry(gdf, gdf2)
    gdf = filter_for_subsites(gdf, filter_ownership)
    check_for_duplicate_polygons(gdf)
    check_for_overlapping_polygons(gdf, gdf2)
    return gdf

###### CHECK DATA ###################################################################################################
def check_for_duplicate_polygons(gdf):
    """Checks for duplicate polygons and raises an error if found."""
    duplicates = gdf[gdf.duplicated(subset='geometry', keep=False)]
    if not duplicates.empty:
        raise ValueError("Duplicate polygons found in the filtered dataset.")

def check_for_overlapping_polygons(gdf, gdf2):
    """Checks for overlapping polygons and raises an error if found."""
    if gdf2 is not None:
        # Find overlaps within the GeoDataFrame by doing a self-overlay
        overlapping_areas = gpd.overlay(gdf2, gdf2, how='intersection')

        # Identify the geometry where there is an overlap in the Maus et al database
        overlapping_geometry = overlapping_areas.geometry
        if not gdf[gdf["geometry"].isin(overlapping_geometry)].empty:
            raise ValueError("Overlapping polygons found in the filtered dataset.")

###### PROCESS DATA #################################################################################################

def process_mining_sites(gdf_path, gdf2_path=gdf2_path, filter_ownership=filter_ownership, by_parent_id=by_parent_id):
    """Processes mining sites and returns a filtered geodataframe and corresponding ownership dataframe."""
    print("\nProcessing mining sites\n")
    # Generate gdf_treated
    gdf = load_data(gdf_path)
    gdf2 = load_data(gdf2_path) if gdf2_path else None

    # Clean gdf
    gdf = get_mine_owners.clean_facility_ids(gdf)
    gdf = allocate_missing_production_start(gdf)
    gdf = allocate_missing_geometry(gdf)
    gdf = remove_duplicated_sites(gdf)

    # Filter data for sites of interest
    gdf_processed = filter_mining_sites(gdf, gdf2, filter_ownership, by_parent_id)

    print("\nMining sites treated\n")

    # Process gdf
    evi = get_satellite(gdf_processed, "facility_id", "production_start", evi_path)

    # Select sites of interest
    evi = evi[evi["facility_id"].isin(gdf_processed["facility_id"])]

    # Generate ownership files based on filtered sites from gdf_treated
    own = get_mine_owners.process_mine_ownership()
    own_processed = own[own["facility_id"].isin(gdf_processed["facility_id"])].drop(columns="id")

    return gdf_processed, own_processed, evi

