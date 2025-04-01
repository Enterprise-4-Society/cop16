from utils import load_data

###### LOAD DATA ####################################################################################################
# Load the PAF cost data
paf_path = 'data/natura/paf/20250314/paf_habitat_costs(data).csv'
encoding = "ISO-8859-1"

# TODO: Check the PAFs

###### CLEAN DATA ###################################################################################################
def clean_paf(df, country):
    # Remove missing information, consider only relevant countries and habitat level costs
    df_treated = df[(df["HABITATCODE"].notna()) & (df["annual_cost_abs[eur]"].notna()) &
                    (df["COUNTRY_CODE"].isin(country)) & (df["associated_sites"].isna())].astype(
        {"annual_cost_abs[eur]": 'int',
         "target":'int'})

    # Treat cost per ha data
    df_treated.loc[(df_treated["COUNTRY_CODE"] == "IE") &
                   (df_treated["measure_type"] == "one-off"), "annual_cost_ha"] = df_treated["annual_cost_abs[eur]"] * 7 / df_treated["target"]
    df_treated.loc[(df_treated["COUNTRY_CODE"] == "IE") &
                   (df_treated["measure_type"] == "recurring"), "annual_cost_ha"] = df_treated["annual_cost_abs[eur]"] / df_treated["target"]

    # The IE shows inconsistencies in the annual cost per ha computation across types of measure.
    # In the measure description it suggests that:
    # for recurring measures -  annual cost per ha =  annual_cost_abs[eur] / target
    # while for one-off measures - annual cost per ha = annual_cost_abs[eur] * 7 / target

    # Find non-numeric rows in costs
    def find_non_numeric_rows(df):
        mask = df.applymap(lambda x: not isinstance(x, (int, float)))
        return df[mask.any(axis=1)]

    rows_to_clean = find_non_numeric_rows(df_treated[["annual_cost_abs[eur]", "target"]])

    if not rows_to_clean.empty:
        print("Costs and hectares are non-numeric. Please fix.")

    return df_treated


###### EXTRACT DATA ##################################################################################################
def get_country(file_path=paf_path):
    """
    Get list of countries available in the PAF cost database.
    """
    df = load_data(file_path, delimiter=";")
    df.columns.values[0] = "COUNTRY_CODE"
    country_list = df[["COUNTRY_CODE"]].drop_duplicates().dropna().values.tolist()

    # Flatten the list of lists
    country_list = [item[0] for item in country_list]

    return country_list


###### PROCESS DATA #################################################################################################
def process_cost(file_path=paf_path):
    """
    Process PAF data and compute site-level annual cost estimates based on habitat and biogeographic composition of
    sites.
    """
    from get_natura_sites import process_data

    # Extract data
    paf = load_data(file_path, delimiter=";")
    dfs, sites_of_interest = process_data()

    country = dfs["NATURA2000SITES"]["COUNTRY_CODE"].drop_duplicates().tolist()
    # Clean paf
    paf = clean_paf(paf, country)

    # Cost by habitat and region type
    cost = paf[["COUNTRY_CODE", "HABITATCODE", "BIOGEOGRAPHICREG", 'measure_descrip', 'annual_cost_ha']]
    cost_hab = cost.groupby(
        ['COUNTRY_CODE', 'HABITATCODE', 'BIOGEOGRAPHICREG'],
        as_index=False
    ).agg({
        'annual_cost_ha': 'sum',
        'measure_descrip': lambda x: ' | '.join(x)
    })
    # Clean data_processed
    data = dfs["NATURA2000SITES"][["SITECODE", "COUNTRY_CODE"]].merge(
        dfs["HABITATS"][["SITECODE", "HABITATCODE", "COVER_HA", "AREAHA"]].merge(
            dfs["BIOREGION"][["SITECODE", "BIOGEOGRAPHICREG"]],
            on="SITECODE", how="left"),
        on="SITECODE", how="left")

    # Compute cost at the site level
    cost_site = data.merge(cost_hab, how="left", on=["COUNTRY_CODE", "HABITATCODE", "BIOGEOGRAPHICREG"])
    cost_site["annual_cost_hab"] = cost_site["annual_cost_ha"] * cost_site["COVER_HA"]
    cost_site = (cost_site.dropna(subset=["annual_cost_ha"])
                 .merge(
        cost_site[["SITECODE", "annual_cost_hab"]].groupby("SITECODE").sum(),
        on="SITECODE", how="left")
                 .rename(columns={"annual_cost_hab_y": "annual_cost_site",
                                                                                                                                                                                   "annual_cost_hab_x": "annual_cost_hab"}))
    return cost_site
