import numpy as np
import xarray as xr
import shapely.geometry as sgeom
from shapely.prepared import prep
from shapely.geometry import Point
import geopandas as gpd

def polygon_to_mask(polygons, lon, lat):
    """
    Generate a mask array of points falling into one or multiple polygons.
    
    Parameters:
        polygons (list or shapely.geometry.Geometry): List of polygons or a single polygon.
        lon (numpy.ndarray): 1D or 2D array of longitude values.
        lat (numpy.ndarray): 1D or 2D array of latitude values.
        
    Returns:
        numpy.ndarray: A 2D Boolean mask array where True indicates the point is within the polygons.
    """
    
    # Validate inputs
    if lon.ndim not in {1, 2} or lat.ndim not in {1, 2}:
        raise ValueError("Longitude and latitude must be either 1D or 2D arrays.")
    
    if lon.ndim == 1 and lat.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)
    elif lon.shape != lat.shape:
        raise ValueError("Longitude and latitude arrays must have the same shape when 2D.")
    
    # Ensure polygons is a list
    if isinstance(polygons, sgeom.base.BaseGeometry):  # Check if it's a single polygon
        polygons = [polygons]

    # Prepare polygons for efficient lookup
    prepared_polygons = [prep(poly) for poly in polygons]
    
    # Flatten lon/lat for vectorized point generation
    flat_lon, flat_lat = lon.ravel(), lat.ravel()
    points = np.column_stack((flat_lon, flat_lat))
    
    # Create mask for all points
    mask_flat = np.zeros(points.shape[0], dtype=bool)
    for prepared in prepared_polygons:
        contains = np.array([prepared.contains(Point(p)) for p in points])
        mask_flat |= contains  # Combine masks for all polygons
    
    # Reshape mask back to the original grid shape
    mask = mask_flat.reshape(lon.shape)
    return mask

def calculate_mean_within_mask(dataset, variable_name, mask,
                                  dims=['latitude', 'longitude']):
    """
    Calculate the mean value of a variable within a polygonal region from a NetCDF file using a shapefile.
    
    Parameters:
        ncfile (str): Path to the NetCDF file.
        shapefile (str): Path to the shapefile.
        variable_name (str): Name of the variable in the NetCDF file.
        lon_var (str): Name of the longitude variable in the dataset (default: "longitude").
        lat_var (str): Name of the latitude variable in the dataset (default: "latitude").
        dims (list of str): Dimensions to calculate the mean over (default: all).
        
    Returns:
        float: Mean value of the variable within the polygonal region.
    """
    # Validate presence of variable
    if variable_name not in dataset:
        raise ValueError(f"Variable '{variable_name}' not found in the dataset.")
    
    # Extract data
    data = dataset[variable_name]
    
    if dims is None:
        dims = data.dims  # Default to all dimensions
    mask_da = xr.DataArray(mask, dims=dims)
    
    # Apply mask to data
    masked_data = data.where(mask_da)
    
    # Calculate mean
    mean_value = masked_data.mean(dim=dims, skipna=True).item()
    
    return mean_value

if __name__ == '__main__':
    # File paths and variable
    ncfile = 'data.nc'
    shapefile = 'shapefile.shp'
    variable_name = 'data'
    lon_var = 'longitude'  # Customize if necessary
    lat_var = 'latitude'   # Customize if necessary
    
    # Read NetCDF data
    dataset = xr.open_dataset(ncfile)
    lon = dataset[lon_var].values
    lat = dataset[lon_var].values
    
    # Read shapefile and combine polygons
    shp = gpd.read_file(shapefile)
    # polygons = shp.geometry[0]
    polygons = shp.geometry.unary_union  # Merge all polygons into a single geometry
    mask = polygon_to_mask(polygons, lon, lat)
    
    # Calculate mean within polygon
    mean_value = calculate_mean_within_mask(dataset, variable_name, mask, dims=[lon_var, lon_var])
    
    # Output result
    print(f"Mean value within polygon: {mean_value}")
    
'''
# mask for CMAQ data

def polygon_to_mask(polygon, x, y):

    # Generate a mask array of points falling into the polygon

    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    mask = np.zeros(x.shape, dtype=bool)

    # if each point falls into a polygon, without boundaries
    prepared = prep(polygon)
    for index in np.ndindex(x.shape):
        point = sgeom.Point(x[index], y[index])
        if prepared.contains(point):
            mask[index] = True

    return mask

if __name__ == '__main__':
    import geopandas as gpd
    
    # read data
    ncfile = xr.open_dataset('data.nc')
    lon = ncfile.longitude
    lat = ncfile.latitude
    data = ncfile.data

    # read shapefile
    shp = gpd.read_file('shapefile.shp')

    # generate mask
    mask = polygon_to_mask(shp.geometry[0],lon,lat)
    mask_da = xr.DataArray(mask,dims=('y','x'))

    # apply mask
    masked_data = data.where(mask_da)
    data_mean   = masked_data.mean(dim=['y','x'],skipna=True)
'''