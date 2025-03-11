# Packages
import numpy as np
import pandas as pd
import os
import geopandas as gpd
import ee
import requests
from shapely.geometry import mapping
import geemap.foliumap as geemap
from tqdm import tqdm
import time
import math


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


class EVI:
    def __init__(self, site_file, shp_path, gdf_path, image_path, evi_path, param_evi, summer, image) -> object:
        """
        :param site_file: dataframe - dataset including the (unique) identifiers of the sites of interest
        :param shp_path: str - path leading to the original shape files
        :param gdf_path: str - path leading to shape files that have already been opened and converted in the
        appropriate crs. The related file has two variables: the site identifier and the geometry.
        :param image_path: str - # path leading to the folder where generated images are stored
        :param evi_path: str - # path leading to the csv file listing the EVI scores of the site of interests and
        previously downloaded sites
        :param param_evi: dictionary - parameters including:
                - "project": the client project ID or number to be used for making API calls on Google Earth Engine
                - "crs": the authority code (e.g. 'EPSG:4326') for the base coordinate system of the satellite images,
                - "satellite": the Google Earth Engine code for the satellite chosen (e.g. 'LANDSAT/LE07/C02/T1'),
                - "scale": resolution of the satellite image - which generally depends on the satellite chosen e.g.
                Landsat has a resolution of 30m so scale=30,
                - "start": (float or dict) year before the treatment (e.g. before the mine started or before the sites
                started being protected) and
                - "end": year after the treatment (e.g. today or most recent images from the satellite database)
        :param summer: boolean - if True, the algorithm only takes images from May 1 to Sept 30, else all year considered
        :param image: boolean - if True, the algorithm downloads all images generated
        """

        # Paths
        self.shp_path = shp_path
        self.gdf_path = gdf_path
        self.image_path = image_path
        self.evi_path = evi_path

        # Datasets
        self.site_file = site_file

        # Parameters
        self.project = param_evi.get("project")
        self.crs = param_evi.get("crs")
        self.satellite = param_evi.get("satellite")
        self.scale = param_evi.get("scale")
        self.start = param_evi.get("start")
        self.end = param_evi.get("end")

        self.summer = summer
        self.get_image = image

        # Variable names
        self.site_code_variable = param_evi.get("variable")

        # Observations of interest
        self.site_codes = pd.Series(self.site_file[self.site_code_variable], name=self.site_code_variable) \
            .drop_duplicates().reset_index(drop=True)

    def get_polygons(self, site_codes):
        """
        Returns a Geodataframe including the geometries of all sites from the site file and that have been used at least once.
        Shapely and EE geometries are both available.
        :param site_codes: data series - the codes of the sites of interest
        :return gdf
        """
        print("\nGetting geometries\n")
        # Read or create a gdf from the shape files
        if os.path.exists(self.gdf_path):  # TODO: replace by df_geometry_wo_EVI
            print("\n> Directory exists:\n")
            # Read the existing GeoDataFrame
            gdf = gpd.read_file(self.gdf_path) # TODO: replace by df_geometry_wo_EVI

            # Check if all the sites of interest are in the GeoDataFrame
            boolean_series = site_codes.isin(gdf[self.site_code_variable])

            # If not, update the file with the sites of interest
            if not boolean_series.all():
                print("\n>> Updating gdf with new sites\n")
                # Get the site codes that are not in the existing GeoDataFrame
                missing_site_codes = site_codes[~boolean_series]
                gdf_all = gpd.read_file(self.shp_path)  # TODO: replace by df_geometry_wo_EVI

                if gdf_all.crs.to_string() != self.crs:
                    print("\n>>> Reprojecting shapefile to EPSG:4326\n")
                    gdf_all = gdf_all.to_crs(epsg=int(self.crs.split(":")[1]))

                # Extract the missing sites of interest
                gdf_new_sites = gdf_all[[self.site_code_variable, "geometry"]].loc[(gdf_all[self.site_code_variable]
                                                                                    .isin(missing_site_codes))]

                # Append the new sites to the existing GeoDataFrame
                gdf = gpd.GeoDataFrame(pd.concat([gdf, gdf_new_sites], ignore_index=True))

                # Export the updated GeoDataFrame
                print("\n>> Gdf exported\n")
                gdf.to_file(self.gdf_path) # TODO: replace by df_geometry_wo_EVI

            # Convert gdf geometry to the Earth Engine system
            gdf["geometry_ee"] = gdf["geometry"].apply(create_ee_geometry)

        else:
            print("\n> Directory does not exist:\n")
            # Create the directory if it does not exist
            os.makedirs(os.path.dirname(self.gdf_path), exist_ok=True)   # TODO: replace by df_geometry_wo_EVI

            print("\n> Creating gdf with sites of interest\n")
            # Read the shape files, extract the sites of interest and export
            gdf_all = gpd.read_file(self.shp_path) # TODO: replace by df_geometry_wo_EVI

            if gdf_all.crs.to_string() != self.crs:
                print("\n>>> Reprojecting shapefile to EPSG:4326\n")
                gdf_all = gdf_all.to_crs(epsg=int(self.crs.split(":")[1]))

            gdf = gdf_all[[self.site_code_variable, "geometry"]].loc[
                (gdf_all[self.site_code_variable].isin(self.site_codes))]

            # Export the GeoDataFrame
            print("\n> Gdf exported\n")
            gdf.to_file(self.gdf_path)  # TODO: replace by df_geometry_wo_EVI

            # Convert gdf geometry to the Earth Engine system
            gdf["geometry_ee"] = gdf["geometry"].apply(create_ee_geometry)

        return gdf[gdf[self.site_code_variable].isin(site_codes)]

    def compute_evi(self, image, area):
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
            "scale": self.scale,
            "maxPixels": 10000000,  # Default number of pixels to reduce
            "bestEffort": True,
            # Boolean which if True compute and use a larger scale which would allow the operation to succeed in case the polygon would contain too many pixels at the given scale
            "crs": self.crs})

        # Get a number (not an object)
        return dict_evi.getInfo()["EVI"]

    def retrieve_image(self, area, start_year, end_year):
        """
        Return the satellite image clipped to the associated area. This image is based on an image collection from "start_date" to "end_date"
        using the satellite "satellite" and the base coordinate system "crs". The function also reduces and corrects for atmospheric events. The number of
        images present in the image collection is also returned.
        :param area: ee.geometry.Geometry - geometry of the area of interest
        :param start_year: int - starting year of the image collection, the span between end_year and start_year is 3
        :param end_year: int - ending year of the image collection
        :return: ee.image.Image, float
        """

        # Check if difference in years is 3 years
        if (end_year - start_year != 2) and (self.summer is True):
            raise ValueError("Invalid timespan")

        # Retrieve the image collection
        if self.summer is True:
            start_month = "-05-01"
            end_month = "-09-30"
            list_year_start = [str(x) + start_month for x in range(start_year, end_year + 1)]
            list_year_end = [str(x) + end_month for x in range(start_year, end_year + 1)]

            landsat_img = ee.ImageCollection(self.satellite) \
                .filter(ee.Filter.Or(
                ee.Filter.date(list_year_start[0], list_year_end[0]),
                ee.Filter.date(list_year_start[1], list_year_end[1]),
                ee.Filter.date(list_year_start[2], list_year_end[2]))) \
                .filterBounds(area)  # Crop the images only for the interested area

        else:
            start_month = "-01-01"
            end_month = "-12-31"
            start_date = str(start_year) + start_month
            end_date = str(end_year) + end_month
            landsat_img = ee.ImageCollection(self.satellite) \
                .filterDate(start_date, end_date) \
                .filterBounds(area)  # Crop the images only for the interested area

        # Count how many images are in this image collection
        image_count = landsat_img.size()

        # Get the image count as a number
        image_count_number = image_count.getInfo()

        # Correct the images for atmospheric events (doing an average between a collection)
        if self.satellite == 'LANDSAT/LE07/C02/T1':
            composite = ee.Algorithms.Landsat.simpleComposite(**{
                'collection': landsat_img, 'asFloat': True})

        else:
            print("The code for other satellites is not yet available.")

        # Return the image cropped for the selected area
        return composite.clip(area), image_count_number

    def get_satellite(self, site_codes):
        """
        Extracts the satellite images for the sites of interest over the period of interest. Returns the EVI score
        and number of images used the computation of the EVI.
        :param site_codes: data series - series of the codes of the sites of interest
        :return: dataframe
        """
        print("\nGetting satellite images:\n")
        # Authenticate to the Earth Engine account
        ee.Authenticate()
        if self.project is None:
            ee.Initialize()
        else:
            ee.Initialize(project=self.project)

        # Extract the geometries of interest by first, updating the gdf_file with potentially missing sites, and second,
        # only selecting the geometries of the sites of interest
        geometries = self.get_polygons(site_codes)

        # Create dataframe for saving vegetation indicator
        evi_df = pd.DataFrame(
            columns=[self.site_code_variable, 'Start year', 'End year', 'EVI score', 'Number of images'])

        # Compute the EVI for each site and for all timespans of interest
        i = 0
        total_sites = len(geometries[self.site_code_variable])
        for i, sitecode in enumerate(tqdm(geometries[self.site_code_variable], total=total_sites)):
            # Define the area of interest
            print("\nSite ", i + 1, "/", len(geometries))
            start_site = time.time()

            area = geometries["geometry_ee"].iloc[i]
            i += 1

            if isinstance(self.start, dict):
                start = self.start.get(sitecode)
            else:
                start = self.start

            for year in [start, self.end]:

                # Skip the iteration if the year is NaN
                if year is None or (isinstance(year, float) and math.isnan(year)):
                    continue

                print(year, "\n")
                print(sitecode)
                start_time = time.time()

                # Set the start and end months
                start_year = year - 2
                end_year = year

                # Retrieve the images
                image, image_count = self.retrieve_image(area, start_year, end_year)

                # Skip sites with no available images
                if image_count == 0:
                    print(f"Skipping site {sitecode} due to no available images for the specified period.")
                    evi_score = np.nan

                else:
                    print(f"Time after retrieving images: {time.time() - start_time} seconds")
                    # Get the EVI score of the location on a given timespan
                    evi_score = self.compute_evi(image, area)

                print(f"Time after computing EVI: {time.time() - start_time} seconds")

                # Append the new data to the main dataframe
                new_row = {
                    self.site_code_variable: sitecode,
                    'Start year': int(start_year),
                    'End year': int(end_year),
                    'Summer only': self.summer,
                    'EVI score': evi_score,
                    'Number of images': image_count
                }
                new_row_df = pd.DataFrame([new_row])
                evi_df = pd.concat([evi_df, new_row_df], ignore_index=True)
                print(f"Total time for year {year}: {time.time() - start_time} seconds")

                # Download image for the sites of interest
                if self.get_image:
                    path_image_year = str(sitecode) + "_" + str(year)

                    self.download_images(image, area, path_image_year, "rgb", tif=True)
                    self.download_images(image, area, path_image_year, "inf", tif=True)
                    print(f"Total time for downloading images {year}: {time.time() - start_time} seconds")

            print(f"Total time for site {sitecode}: {time.time() - start_site} seconds\n")

        return evi_df

    def get_evi(self):
        """
        Returns the EVI for the start and end dates of interest, for the sites of interest.
        :return: dataframe
        """
        print("\nGetting EVI values:\n")
        if os.path.exists(self.evi_path):
            print("\n > Directory exists:\n")
            # Note: Run time for generating the evi score is very high, so we export and automatize opening if it has
            # already been generated
            evi_df = pd.read_csv(self.evi_path, index_col=0).reset_index()

            # Check if the current version of the EVI file includes all sites of interest
            evi_pilot_sites = set(evi_df[(evi_df["Start year"] == self.start) & (evi_df["End year"] == self.end)][self.site_code_variable].unique())
            site_code_pilot_forest_set = set(self.site_codes)

            # Find missing sites
            missing_sites = site_code_pilot_forest_set.difference(evi_pilot_sites)

            if missing_sites:
                print(">> Missing sites and associated site codes: ", missing_sites)

                # Convert site_codes to a pandas Series
                missing_sites_series = pd.Series(list(missing_sites))

                # Get the EVI for the missing sites
                evi_df_new = self.get_satellite(missing_sites_series)

                # Append the EVI to the main EVI file
                evi_df_updated = pd.concat([evi_df, evi_df_new], ignore_index=True)  # TODO maybe add _final
                print("\n>> Updating EVI file\n")
                evi_df_updated.to_csv(self.evi_path, index=False)
                evi_df = evi_df_updated

            evi_df = evi_df.loc[evi_df[self.site_code_variable].isin(self.site_codes)]

        else:
            print("\n> Directory does not exist:\n")
            # Check if the directory exists and if not create it
            directory = os.path.dirname(self.evi_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

            # Get the EVI for all the sites in the sample
            evi_df = self.get_satellite(self.site_codes)
            print("\n >> Saving EVI file\n")
            evi_df.to_csv(self.evi_path, index=False)

        return evi_df

    def get_images(self, sitecode, year, param_image, tif=False):
        """
        Returns a single image for a specific site at a specific year
        # TODO remove the downloading part in the get satellite images and separate the download of images in another
        # function: call for the GEE data with Initialize et
        :param sitecode: str - code for the site of interest
        :param year: int - year of interest
        :param param_image: str - file name for the image
        :param tif: boolean - if True, downloads also .tif satellite image
        :return: .jpg file and if tif is True, .tif file
        """
        # Authenticate to the Earth Engine account
        ee.Authenticate()
        ee.Initialize(project=self.project)

        # Get area
        geometries = self.get_polygons(sitecode)
        area = geometries.loc[geometries[self.site_code_variable] == sitecode, "geometry_ee"].iloc[0]

        # Set the start and end months
        start_year = year - 2
        end_year = year

        # Get satellite images
        image, image_count = self.retrieve_image(area, start_year, end_year)

        # Download image
        path = str(sitecode) + "_" + str(year)
        return self.download_images(image, area, path, param_image, tif=tif)

    def download_images(self, image, area, path_image_year, param_image, tif=False):
        """
        Downloads the .jpg and .tif (optional) satellite images of an area of interest using RGB or infrared bands
        :param path_image_year:str - file name for the image
        :param image: ee.image.Image - satellite image to be downloaded according to the parameters specification (time and location set already
        :param area: ee.geometry.Geometry - geometry of the area of interest
        :param param_image: str - "rgb" for downloading RGB images or "inf" for downloading infrared images
        :param tif: boolean - if True, downloads also .tif satellite image
        :return: .jpg and optionally .tif images of the area of interest
        """
        if param_image == "rgb":
            # Set RGB satellites bands parameters
            # TODO Understand the options for the band definition
            rgb = {'bands': ["B3", "B2", "B1"], 'dimensions': "500x500", 'min': 0, 'max': 0.3}

            # Download in .jpg
            print("Downloading RGB in JPG")
            url = image.getThumbURL(rgb)
            img_data = requests.get(url).content
            name = self.image_path + "_" + path_image_year + "_" + param_image
            with open(name + '.jpg', "wb") as handler:
                handler.write(img_data)

            # Download in .tif
            if tif:
                print("Downloading RGB in TIF")
                geemap.download_ee_image(image, region=area, crs=self.crs, scale=self.scale,
                                         filename=self.image_path + "_" + path_image_year + "_" + param_image + ".tif")

        if param_image == "inf":
            # Set INF, G, B satellites bands parameters
            # TODO Understand the options for the band definition
            inf = {'bands': ["B4", "B2", "B1"], 'dimensions': "500x500", 'min': 0, 'max': [0.5, 0.3, 0.3]}

            print("Downloading INF in JPG")
            # Download in .jpg
            url = image.getThumbURL(inf)
            img_data = requests.get(url).content
            name = self.image_path + "_" + path_image_year + "_" + param_image
            with open(name + '.jpg', "wb") as handler:
                handler.write(img_data)

            # Download in .tif
            if tif:
                print("Downloading INF in TIF")
                geemap.download_ee_image(image, region=area, crs=self.crs, scale=self.scale,
                                         filename=self.image_path + "_" + path_image_year + "_" + param_image + ".tif")

        return None
