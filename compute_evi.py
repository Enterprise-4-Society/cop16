import ee
import os
import csv
from tqdm import tqdm
import time
import numpy as np
import pandas as pd
from shapely.geometry import mapping

###### DEFINE PARAMETERS ###########################################################################################
# Define ee parameters
crs = 'EPSG:4326'
ee_project = "ee-e4s-biodiversity-loss"
satellite = "LANDSAT/LE07/C02/T1"
span = 2
scale = 30
summer = True
post_date = 2023

# Authenticate to the Earth Engine account
ee.Authenticate()
ee.Initialize(project=ee_project)

###### GET POLYGONS #################################################################################################
def create_ee_geometry(shapely_geometry):
    """
    Returns the geometry of a GeoDataFrame that can be used in Earth Engine
    :geometry: Shapely geometry to convert to ee.geometry
    return - ee.geometry.Geometry
    """

    try:
        # Check if the geometry is valid
        if not shapely_geometry.is_valid:
            raise ValueError("Invalid Shapely geometry")

        # Convert Shape geometry to GeoJSON (compatible with EE)
        geojson_geom = mapping(shapely_geometry)

        # Create and return Earth Engine Geometry
        return ee.Geometry(geojson_geom)

    except Exception as e:
        # Handle any exceptions
        print(f"Error converting geometry: {e}")
        return None


def get_polygons(gdf_to_convert, crs=crs):
        """
        Ensure that gdf geometries have the right characteristics
        :return: gdf
        """

        # Verify if the gdf applies the right crs
        if gdf_to_convert.crs.to_string() != crs:
            print(f"\n>>> Reprojecting shapefile to {crs}\n")
            gdf_converted = gdf_to_convert.to_crs(epsg=int(crs.split(":")[1]))

        else:
            gdf_converted = gdf_to_convert

        print(f"\n>>> Converting shapefile\n")
        # Convert gdf geometry to the Earth Engine system
        gdf_converted["geometry_ee"] = gdf_converted["geometry"].apply(create_ee_geometry)

        return gdf_converted # try by implementing .simplify(tolerance=0.01)

###### GET UNPROCESSED DATE-SITE PAIRS ##############################################################################

def add_post_event(gdf, date_label, post_date=post_date):
    # Copy existing dataframe and change date with post date
    gdf_post = gdf.copy()
    gdf_post[date_label] = post_date

    # Concate gdfs
    gdf_with_post = pd.concat([gdf, gdf_post])

    return gdf_with_post


def remove_processed_sites(df, processed_sites_path, site_code_label, date_label):
    # Define sites already processed
    processed_sites = pd.read_csv(processed_sites_path)[[site_code_label, "End year"]].rename({"End year": date_label}, axis="columns")

    # Define sites that need to be processed
    to_process_sites = df[[site_code_label, date_label]]

    # Merge with indicator to track row existence
    merged = to_process_sites.merge(processed_sites, on=list(to_process_sites.columns), how='outer', indicator=True)

    # Get date-site pairs that have not been processed yet
    sites = merged[merged["_merge"] != "both"][[site_code_label, date_label]]

    # Return error if all sites have been processed
    if sites.empty:
        print("All sites have been already processed.")

    df_filtered = df.merge(sites, on=[site_code_label, date_label], how="inner")

    return df_filtered

def remove_error_sites(df, site_code_label):
    sites_with_errors = ["NL9801023", "ES1110016", "FI0100005", "SE0810436"]
    df_treated = df[~df[site_code_label].isin(sites_with_errors)]
    return df_treated

###### RETRIEVE SATELLITE IMAGES ####################################################################################

def retrieve_image(area, start_year, end_year, summer=summer, satellite=satellite):
    """
    Return the satellite image clipped to the associated area. This image is based on an image collection from "start_date" to "end_date"
    using the satellite "satellite" and the base coordinate system "crs". The function also reduces and corrects for atmospheric events. The number of
    images present in the image collection is also returned.
    :param area: ee.geometry.Geometry - geometry of the area of interest
    :param start_year: int - starting year of the image collection, the span between end_year and start_year is 3
    :param end_year: int - ending year of the image collection
    :return: ee.image.Image, float
    """

    # Retrieve the image collection
    if summer is True:
        start_month = "-05-01"
        end_month = "-09-30"
        list_year_start = [str(x) + start_month for x in range(start_year, end_year + 1)]
        list_year_end = [str(x) + end_month for x in range(start_year, end_year + 1)]

        landsat_img = ee.ImageCollection(satellite) \
            .filter(ee.Filter.Or(
            ee.Filter.date(list_year_start[0], list_year_end[0]),
            ee.Filter.date(list_year_start[1], list_year_end[1]),
            ee.Filter.date(list_year_start[2], list_year_end[2]))) \
            .filterBounds(area)  # Crop the images only for the interested area

        limited_images = landsat_img.sort('CLOUD_COVER').limit(70)

    else:
        start_month = "-01-01"
        end_month = "-12-31"
        start_date = str(start_year) + start_month
        end_date = str(end_year) + end_month
        landsat_img = ee.ImageCollection(satellite) \
            .filterDate(start_date, end_date) \
            .filterBounds(area)  # Crop the images only for the interested area

        limited_images = landsat_img.sort('CLOUD_COVER').limit(70)

    # Count how many images are in this image collection
    image_count_number = limited_images.size().getInfo()

    # Correct the images for atmospheric events (doing an average between a collection)
    if satellite == 'LANDSAT/LE07/C02/T1':
        composite = ee.Algorithms.Landsat.simpleComposite(**{
            'collection': limited_images, 'asFloat': True})
        print("The code is available")

    else:
        print("The code for other satellites is not yet available.")

    # Return the image cropped for the selected area
    return composite.clip(area), image_count_number

###### COMPUTE EVI #################################################################################################

def compute_evi(image, area, crs, scale):
        """
        Returns the Enhanced Vegetation Index by NASA of an image  "image" and associated area "area" using the base coordinate system "crs". "Scale"
        represents the pixel resolution at which the image will be reduced.
        :param image: ee.image.Image - satellite image for which the EVI will be calculated
        :param area: ee.geometry.Geometry - geometry of the area of interest
        :return: float
        """

        # Select the bands
        nir = image.select('B4')
        red = image.select('B3')
        blue = image.select('B1')

        # Get EVI for all pixels
        EVI = image.expression(
            '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
                'NIR': nir,
                'RED': red,
                'BLUE': blue
            }).rename("EVI")

        # Take mean of EVI
        dict_evi = EVI.reduceRegion(**{
            "reducer": ee.Reducer.mean(),
            "geometry": area,
            "scale": scale,
            "maxPixels": 10000000,  # Default number of pixels to reduce
            "bestEffort": True,
            # Boolean which if True compute and use a larger scale which would allow the operation to succeed in case the polygon would contain too many pixels at the given scale
            "crs": crs})

        # Get a number (not an object)
        return dict_evi.getInfo()["EVI"]

###### PROCESS EVI #################################################################################################

def get_satellite(gdf, site_code_label, date_label, output_folder, crs=crs, span=span, scale=scale, summer=summer):
    """
    output_csv (str): folder in which the evi file should be created
    """

    # Get ee geometry
    geometries = get_polygons(gdf, crs)

    # Add the post date
    geometries = add_post_event(geometries, date_label, post_date)

    # Check if the file already exists
    final_path = output_folder + f'evi_{post_date}_{span}y_{"summer" if summer is True else ""}.csv'
    file_exists = os.path.isfile(final_path)

    # Remove date-site pairs already processed in the output file
    if file_exists:
        geometries = remove_processed_sites(geometries, final_path, site_code_label, date_label)
        if geometries.empty:
            print("Opening existing file")
            return pd.read_csv(final_path)

    # Remove sites with errors
    geometries = remove_error_sites(geometries, site_code_label)

    print("\nSite geometries processed\n")

        # Open the CSV file for writing
    with open(final_path, mode="a" if file_exists else "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [site_code_label, "End year", "Start year", "Summer only", "EVI score", "Number of images"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write header only if the file is new
        if not file_exists:
            writer.writeheader()

        total_sites = len(geometries[site_code_label])
        for i, sitecode in enumerate(tqdm(geometries[site_code_label], total=total_sites)):
            # Define the area of interest
            print("\nSite ", i + 1, "/", total_sites, " : ", sitecode, "\n")
            start_time = time.time()

            area = geometries["geometry_ee"].iloc[i]
            end_year = geometries[date_label].iloc[i]  # Retrieve the associated year for the site
            start_year = end_year - span

            print(end_year, "\n")

            # Retrieve the images
            image, image_count = retrieve_image(area, start_year, end_year)

            # Skip sites with no available images
            if image_count == 0:
                print(f"Skipping site {sitecode} due to no available images for the specified period.")
                evi_score = np.nan

            else:
                print(f"Time after retrieving images: {time.time() - start_time} seconds")
                # Get the EVI score of the location on a given timespan
                evi_score = compute_evi(image, area, crs, scale)

            print(f"Time after computing EVI: {time.time() - start_time} seconds")
            # Append the new data to the main csv
            new_row = {
                site_code_label: sitecode,
                'End year': int(end_year),
                'Start year': int(start_year),
                'Summer only': summer,
                'EVI score': evi_score,
                'Number of images': image_count
            }

            writer.writerow(new_row)
            i += 1

        return pd.read_csv(final_path)