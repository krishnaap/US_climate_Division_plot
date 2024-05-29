#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 18:30:38 2024

@author: krishna
"""


import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from shapely.geometry import Polygon, MultiPolygon, Point
from shapely.ops import unary_union
import numpy as np
import cartopy.feature as cfeature
import sys



# Load the datasets
era5_data = xr.open_dataset('./data/US_airtemp_era5.nc')
climdiv_polygons = xr.open_dataset('./data/climdiv_polygons.nc')

# =============================================================================
# Plot different Climate divisions based on the State and region
# https://psl.noaa.gov/data/usclimdivs/descript.html
# =============================================================================

def get_polygon_from_division(division):
    latitudes = division.attrs['lat']
    longitudes = division.attrs['lon']
    return Polygon(zip(longitudes, latitudes))

fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.LambertConformal()})
ax.set_extent([-122, -68, 25, 50], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE, linewidth = 1)
ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth = 1)
ax.add_feature(cfeature.STATES, linestyle=':', linewidth = 0.3)
ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

NGP = ['MT', 'WY', 'ND', 'SD', 'NE']                                            #North Great Plains
SGP = ['CO', 'NM', 'KS', 'TX', 'OK']                                            #South Great Plains
MW = ['MN', 'IA', 'MO', 'IL', 'WI', 'MI', 'IN', 'OH', 'KY']                     #MidWest
SE = ['AR', 'LA', 'AL', 'FL', 'GA', 'MS', 'NC', 'SC', 'TN']                     #Southeast
NE = ['VA', 'WV', 'PA', 'MD', 'DE', 'NJ', 'NY', 'CT', 'RI', 'MA', 'NH', 'VT', 'ME'] #Northeast

regions = [NGP, SGP, MW, SE, NE]

for region in regions:
    clim_divisions = [f"{state}_CD{num}" for state in region for num in range(1, 11)]
    polygons = [get_polygon_from_division(climdiv_polygons[div]) for div in clim_divisions if div in climdiv_polygons.data_vars]
    combined_polygon = unary_union(polygons)

    ax.add_geometries([combined_polygon], ccrs.PlateCarree(), facecolor='none', edgecolor='blue', linewidth=2)

plt.show()
#%%

# Define the function to create polygons from division data
def get_polygon_from_division(division):
    latitudes = division.attrs['lat']
    longitudes = division.attrs['lon']
    return Polygon(zip(longitudes, latitudes))

# List of relevant division codes for North Great Plains from specified states
states = ['MT', 'WY', 'ND', 'SD', 'NE']  # States in the North Great Plains
ngp_divisions = [f"{state}_CD{num}" for state in states for num in range(1, 11)]

# Extract polygons and combine into a single MultiPolygon
polygons = [get_polygon_from_division(climdiv_polygons[div]) for div in ngp_divisions if div in climdiv_polygons.data_vars]
combined_polygon = unary_union(polygons)

# Create a mask for each ERA5 data point
latitudes = era5_data.latitude.values
longitudes = era5_data.longitude.values
mask = np.fromiter((combined_polygon.contains(Point(lon, lat)) for lat in latitudes for lon in longitudes), dtype=bool)
mask = mask.reshape((len(latitudes), len(longitudes)))

# Apply mask to the t2m data for a specific time (e.g., first time point)
t2m_masked = era5_data.t2m.isel(time=0).where(mask)

# Plotting
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_extent([-125, -65, 25, 50], crs=ccrs.PlateCarree())  # Continental USA extent

# Add map features for context
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
# ax.add_feature(cfeature.STATES, linestyle='--')

# Plot masked t2m data
t2m_masked.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='coolwarm', add_colorbar=True)

# Outline the North Great Plains
ax.add_geometries([combined_polygon], ccrs.PlateCarree(), facecolor='none', edgecolor='blue', linewidth=2)

ax.set_title('ERA5 2m Temperature over North Great Plains')
plt.show()



#%%# Function to create a spatial mask as an xarray.DataArray
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
def create_spatial_mask(era5_data, combined_polygon):
    # Create a mesh grid of longitude and latitude coordinates
    lon, lat = np.meshgrid(era5_data.longitude, era5_data.latitude)
    # Check if each point is within the combined polygon
    mask = np.array([[combined_polygon.contains(Point(lon_val, lat_val))
                      for lon_val in lon_row] for lat_val, lon_row in zip(lat, lon)])

    return xr.DataArray(mask, coords=[era5_data.latitude, era5_data.longitude], dims=["latitude", "longitude"])

# Create the mask
mask = create_spatial_mask(era5_data, combined_polygon)

# Apply the mask to the data
masked_data = era5_data.t2m.where(mask, drop=True)

# Compute the spatial mean over time
spatial_mean = masked_data.mean(dim=["latitude", "longitude"])

# Plot the time series of the spatial mean
plt.figure(figsize=(10, 6))
spatial_mean.plot()
plt.title('Time Series of Spatial Mean Temperature over North Great Plains')
plt.xlabel('Time')
plt.ylabel('Temperature (K)')
plt.grid(True)
plt.show()
#%%
t2m_values = spatial_mean.t2m.values  # This extracts the numpy array directly

# Extract time values for plotting
time_values = spatial_mean.time.values

# Create the plot
plt.figure(figsize=(12, 6))
plt.plot(time_values, t2m_values, label='Spatial Mean Temperature', color='blue')
plt.title('Time Series of Spatial Mean Temperature over North Great Plains')
plt.xlabel('Time')
plt.ylabel('Temperature (K)')
plt.grid(True)
plt.legend()
plt.show()


#%% Deprecated
