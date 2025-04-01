###### DEFINE PARAMETERS #############################################################################################
post_date = 2023

# Load the Jasanski geodataframe (main shape file)
gdf_path = 'data/mining/open_database_mine_production/data/facilities.gpkg'

###### DEFINE TOOLS ##################################################################################################

def get_evi_diff(df, site_id):
    """Compute the difference in EVI scores between two time points for each site. This function assumes each site_id
    appears exactly twice (i.e., for two years), and computes the difference as the second EVI value minus the first.

    Parameters:
        df : pandas.DataFrame - The dataframe containing at least 'site_id' and 'EVI score' columns.
        site_id : str - The name of the column identifying unique sites.

    Returns:
        pandas.DataFrame - A dataframe with columns [site_id, 'EVI score difference'].
    """
    evi_diff = df.groupby(site_id).agg({
        'EVI score': lambda x: x.iloc[1] - x.iloc[0]
    }).rename(columns={'EVI score': 'EVI score difference'}).reset_index()

    return evi_diff

###### COMPUTE COST PER EVI OF NATURA2000 SITES ######################################################################

def compute_cost_evi():
    """
    Compute the cost per EVI gain for Natura 2000 sites. This function loads EVI data for Natura 2000 sites, computes
    the EVI gain, merges it with site-level cost data, and estimates the cost per EVI gained

    Returns:
        df : pandas.DataFrame - The dataframe containing per-site cost, EVI difference, and cost per EVI gain.
        average_cost_evi : float - The average cost per unit of EVI gain across all considered Natura 2000 sites.
    """
    from compute_cost import process_cost
    from get_natura_sites import process_natura_sites

    # Get data
    cost_natura = process_cost()[["SITECODE", "annual_cost_site"]].drop_duplicates()
    data_natura, gdf_natura, evi_natura = process_natura_sites()

    del data_natura, gdf_natura

    # Compute EVI difference
    evi_natura = evi_natura.sort_values(by=['SITECODE', 'End year'])

    # This assumes each site_id appears exactly twice
    evi_diff = get_evi_diff(evi_natura, "SITECODE")

    # Merge the EVI difference back to the original dataframe
    evi_natura = evi_natura.merge(evi_diff, on='SITECODE')

    # Compute cost per EVI gain
    df = cost_natura.merge(evi_natura[evi_natura["End year"] != 2023][["SITECODE", "EVI score difference", "End year"]].rename(columns={"End year": "Year"}), on="SITECODE", how="left")

    df = df[df["EVI score difference"] > 0]
    df["cost_per_evi"] = df["annual_cost_site"] * (post_date - df["Year"]) / df["EVI score difference"]

    average_cost_evi = df["cost_per_evi"].mean()

    return df, average_cost_evi

###### COMPUTE COST OF RESTORATION FOR MINES #########################################################################

def compute_cost_mining():
    """
    Estimate restoration cost equivalents for mining sites using average Natura EVI costs. This function retrieves EVI
    scores for mining sites, computes the EVI score difference per site, and estimates the environmental restoration
    cost equivalent based on the average cost per EVI point derived from Natura 2000 sites.

    Returns:
        evi_mine : pandas.DataFrame - A dataframe containing mining site IDs, EVI differences, production start year,
        and estimated restoration cost equivalents. If the EVI improved, the cost is set to zero.
        """
    # Get cost values
    cost_natura, average_cost_evi = compute_cost_evi()
    del cost_natura

    # Get mines EVI difference
    from get_mining_sites import process_mining_sites
    mine_gdf, mine_own, evi_mine = process_mining_sites(gdf_path)

    evi_diff = get_evi_diff(evi_mine, "facility_id")

    # Merge the EVI difference back to the original dataframe
    evi_mine = evi_mine.merge(evi_diff, on='facility_id')

    # Clean to keep only sites with start of the production
    evi_mine = evi_mine[evi_mine["End year"] != post_date].drop(columns=["Start year", "Summer only", "Number of images", "EVI score"]).rename(columns={"End year": "production_start"})

    # Compute the cost
    evi_mine["cost_evi"] = evi_mine["EVI score difference"] * average_cost_evi

    # Change cost to 0 if EVI improved
    evi_mine.loc[evi_mine["EVI score difference"] > 0, "cost_evi"] = 0

    return evi_mine

###### GET TREATED DATA #############################################################################################

natura, average = compute_cost_evi()
mine = compute_cost_mining()

natura.to_csv("output/natura/20250314/cost_per_site.csv", index=False)
mine.to_csv("output/mining/cost_per_site.csv", index=False)



