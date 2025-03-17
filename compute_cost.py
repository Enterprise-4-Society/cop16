from itertools import count

import pandas as pd
import os
from utils import load_data

###### LOAD DATA ####################################################################################################
# Load the PAF cost data
paf_path = 'data/natura/paf/20250314/paf_habitat_costs(data).csv'
encoding = "ISO-8859-1"

# TODO: Check the PAFs
###### GET DATA ####################################################################################################
def get_country():
    df = load_data(file_path=paf_path, encoding=encoding, delimiter=";")
    country_list = df[["COUNTRY_CODE"]].drop_duplicates().dropna().values.tolist()

    # Flatten the list of lists
    country_list = [item[0] for item in country_list]

    return country_list


###### PROCESS DATA #################################################################################################
def process_cost(file_path=paf_path):
    from get_natura_file import get_natura_sites

    df = load_data(file_path)

    sites_of_interest = get_natura_sites()

    return df













# # Paths
# paf_path = r'../data/cost/paf_habitat_costs_20241126.xlsx'
# site_path = r'../data/natura2000/sites/natura2000sites_forest.csv'
# output_path = r'../ouput/natura2000/cost_site.csv'
#
# # Ensure the output directory exists
# os.makedirs(os.path.dirname(output_path), exist_ok=True)
#
# # Open the selected sites and the PAF costs
# paf = pd.read_excel(paf_path, sheet_name="data")
# sites = pd.read_csv(site_path)
#
# # Select only the Ireland data from the PAF dataset & sites
# paf_ie = paf[paf["COUNTRY_CODE"] == "IE"]
# sites_ie = sites[sites["COUNTRY_CODE"] == "IE"]
#
#
# paf_ie["annual_cost_ha"] = paf_ie["annual_cost_abs[eur]"] / paf_ie["target"]
#
# # Compute the cost of country-habitat-geographic region
# habitat_costs = (paf_ie[paf_ie["associated_sites"].isna()][["COUNTRY_CODE", "HABITATCODE", "BIOGEOGRAPHICREG", "annual_cost_ha"]]
#                  .groupby(["COUNTRY_CODE", "HABITATCODE", "BIOGEOGRAPHICREG"])
#                  .sum()
#                  .reset_index())
#
# # Merge the sites data and the costs at the site level
# final = (sites_ie[["SITECODE", "COUNTRY_CODE", "HABITATCODE", "BIOGEOGRAPHICREG", "COVER_HA"]]
#          .merge(paf_ie[['COUNTRY_CODE', 'HABITATCODE', 'BIOGEOGRAPHICREG', 'annual_cost_ha']],
#                       left_on=['COUNTRY_CODE', 'HABITATCODE', 'BIOGEOGRAPHICREG'],
#                       right_on=['COUNTRY_CODE', 'HABITATCODE', 'BIOGEOGRAPHICREG'],
#                       how='left'))
#
# final["annual_cost"] = final["annual_cost_ha"] * final["COVER_HA"]
#
# # Export
#
# final.to_csv(output_path, index=False)