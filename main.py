import geopandas as gpd
import matplotlib.pyplot as plt
import warnings
import pandas as pd
from shapely.geometry import Polygon, LineString
from tkinter import Tk, filedialog, simpledialog
from tkinter.filedialog import askopenfilename
from tkinter.simpledialog import askinteger

# Function to convert a Polygon into LineString
def polygon_to_linestring(polygon):
    exterior_coords = list(polygon.exterior.coords)
    return LineString(exterior_coords)

# Function to convert all polygons in a GeoDataFrame into LineStrings
def convert_polygons_to_linestrings(gdf):
    linestrings = []
    for _, row in gdf.iterrows():
        linestrings.append(polygon_to_linestring(row['geometry']))
    return gpd.GeoDataFrame(geometry=linestrings, crs=gdf.crs)


# Function to ask for file selection with a custom message
def custom_askopenfilename(title, message):
    root = Tk()
    root.withdraw()

    # Custom dialog with message
    file_path = filedialog.askopenfilename(title=title, filetypes=[("Shapefile", ".shp")], initialdir="./")
    if not file_path:
        # If the user cancels the dialog or doesn't select a file, return None
        return None

    # Show the custom message
    simpledialog.messagebox.showinfo(title, message)

    root.destroy()
    return file_path

# Function to remove a specific part from a polygon based on pkuid
def remove_polygon_part(polygon, pkuid_to_remove):
    # Check if the part with pkuid_to_remove exists
    part_to_remove = find_part_with_pkuid(polygon, pkuid_to_remove)
    if part_to_remove and pkuid_to_remove != 100:  # Prevent removal of part with pkuid equal to 100
        # Create a copy of the original exterior
        exterior = list(polygon.exterior.coords)
        
        # Remove the specific part from the interiors
        interiors = [list(interior.coords) for interior in polygon.interiors if interior != part_to_remove]

        # Reconstruct the polygon with the remaining parts
        if exterior and interiors:
            return Polygon(exterior, interiors)
        elif exterior:
            return Polygon(exterior)
    return polygon  # Return the original polygon if part_to_remove is not found or pkuid_to_remove is 100

# Function to find the part with pkuid equal to 100
def find_part_with_pkuid(polygon, pkuid_to_find):
    for idx, interior in enumerate(polygon.interiors):
        if idx + 1 == pkuid_to_find:
            return interior
    return None

# Function to remove subdivisions based on criteria
def remove_subdivision(subdivisions, pkuid_to_remove):
    remaining_subdivisions = []
    polygon_removed = False
    for subdivision in subdivisions:
        if should_keep_subdivision(subdivision, pkuid_to_remove) or polygon_removed:
            remaining_subdivisions.append(subdivision)
        elif not polygon_removed:
            polygon_removed = True
            print("\nRemoved a polygon.")
    return remaining_subdivisions

# Function to determine whether to keep or remove a subdivision
def should_keep_subdivision(subdivision, pkuid_to_remove):
    # Implement your criteria for deciding whether to keep or remove a subdivision
    # For example, remove a specific part of a polygon based on pkuid
    return subdivision['pkuid'] != pkuid_to_remove

# Function to merge two polygons based on their pkuid values
def merge_polygons(gdf, pkuid1, pkuid2, new_pkuid):
    polygon1 = gdf[gdf['pkuid'] == pkuid1].unary_union
    polygon2 = gdf[gdf['pkuid'] == pkuid2].unary_union

    # Check if polygon1 and polygon2 are not None before performing the union
    if polygon1 is None or polygon2 is None:
        print(f"Polygon with pkuid {pkuid1} or {pkuid2} not found.")
        return gdf

    merged_polygon = polygon1.union(polygon2)

    # Create a new GeoDataFrame with the merged polygon and the new pkuid
    merged_gdf = gpd.GeoDataFrame({'geometry': [merged_polygon], 'pkuid': [new_pkuid]}, crs=gdf.crs)

    # Copy attributes from pkuid1 to merged_gdf
    for col in gdf.columns:
        if col != 'geometry' and col != 'pkuid':
            merged_gdf[col] = gdf[gdf['pkuid'] == pkuid1][col].iloc[0]

    # Drop the rows with pkuid1 and pkuid2 from the original GeoDataFrame
    gdf_remaining = gdf[gdf['pkuid'] != pkuid1]
    gdf_remaining = gdf_remaining[gdf_remaining['pkuid'] != pkuid2]

    # Concatenate the original GeoDataFrame with the merged GeoDataFrame
    gdf_remaining = gpd.GeoDataFrame(pd.concat([gdf_remaining, merged_gdf]))

    return gdf_remaining


# Function to add area_before_removal column
def add_area_before_removal(gdf):
    gdf['area_before_removal'] = gdf.apply(lambda row: row['geometry'].area, axis=1)
    return gdf

# Open file explorer to select shapefile
Tk().withdraw()
shapefile_path = askopenfilename(filetypes=[("Shapefile", ".shp")])

if shapefile_path:
    # Read the shapefile
    gdf = gpd.read_file(shapefile_path)

    # Get the CRS of the shapefile
    crs = gdf.crs

    # Get the units from the CRS and append "square" prefix
    units = "square " + crs.axis_info[0].unit_name

    # Conversion factor from square meters to acres
    conversion_factor = 0.000247105

    # Initialize gdf_remaining with a copy of the original GeoDataFrame
    gdf_remaining = gdf.copy()

    # Add the 'area_after_removal' column to the DataFrame with default values
    gdf['area_after_removal'] = 0

    #Path file access
    Tk().withdraw()
    shapefile_path_temp = custom_askopenfilename("Open the path", "Select the path shapefile (*.shp)") # Replace with the path to your shapefile

    # Ask the user for the operation to perform (1: Area Calculation, 2: Removal of Polygon, 3: Merge of Polygon)
    operation = askinteger("Select Operation", "Choose the operation:\n1: Area Calculation\n2: Removal of Polygon\n3: Merge of Polygon")
    
    # Initialize pkuid_to_remove to None
    pkuid_to_remove = None

    if operation == 1:
        # Area Calculation

        # Calculate the areas of individual polygons and store them in a new column
        gdf = add_area_before_removal(gdf)

        # Calculate the total usable and unusable land areas
        total_usable_acres = gdf[gdf['pkuid'] != 100]['area_before_removal'].sum() * conversion_factor
        total_unusable_acres = (gdf['area_before_removal'].sum() - gdf[gdf['pkuid'] != 100]['area_before_removal'].sum()) * conversion_factor

        # Print the total areas
        print("\nTotal Usable Land Area:")
        print(f"In Acres: {total_usable_acres:.2f} acres")

        print("\nTotal Unusable Land Area:")
        print(f"In Acres: {total_unusable_acres:.2f} acres\n")


        #path Area Calculation
        
        # Function to calculate the area of the rectangle given two parallel LineStrings and width
        def calculate_rectangle_area(shapefile_path_temp, width1):
            # Read the shapefile into a GeoDataFrame
            gdf1 = gpd.read_file(shapefile_path_temp)

            # Buffer each LineString to create a polygon (rectangle)
            buffered_geometry = gdf1['geometry'].buffer(width1)

            # Create a new GeoDataFrame with the buffered polygons
            gdf_buffered1 = gpd.GeoDataFrame(geometry=buffered_geometry, crs=gdf1.crs)

            # Calculate the area of the polygon (rectangle) in square units of the CRS
            area_square_units1 = gdf_buffered1.unary_union.area

            # Convert the area to acres
            area_acres1 = area_square_units1 * conversion_factor

            return area_acres1
        
        #Each part Calculation using pkuid number
        for index, row in gdf.iterrows():
            pkuid = row['pkuid']
            area_square_units_before = row['area_before_removal']
            area_acres_before = area_square_units_before * conversion_factor  # Convert area to acres
            print(f"PKUID {pkuid} area : {area_acres_before:.15f} acres")


        #path Area Calculation
        width = 2  # Width in meters
        if shapefile_path_temp.endswith('.shp'):
            area_acres1 = calculate_rectangle_area(shapefile_path_temp, width)
            print(f"\nArea of the Path: {area_acres1:.2f} acres\n")
        else:
            print("\nInvalid shapefile format. Please provide a valid shapefile with two parallel LineStrings.")

    elif operation == 2:
        # Removal of Polygon

        # Add area_before_removal column to the DataFrame
        gdf = add_area_before_removal(gdf)

         # Remove subdivisions based on criteria
        subdivisions_to_remove = gdf.to_dict('records')
        remaining_subdivisions = remove_subdivision(subdivisions_to_remove, pkuid_to_remove)

        # Ask the user for the pkuid to remove (excluding 100)
        while True:
            pkuid_to_remove = askinteger("Remove Polygon", "Enter the pkuid of the part to remove (excluding 100): ")
            if pkuid_to_remove is None:
                print("Invalid value. Please enter a valid pkuid.")
            elif pkuid_to_remove == 100:
                print("Invalid value. Cannot remove the part with pkuid = 100. Please enter a different pkuid.")
            elif pkuid_to_remove not in gdf['pkuid'].values:
                print(f"Part with pkuid {pkuid_to_remove} is already removed. Please enter a different pkuid.")
            else:
                break

        # Calculate the area before removal and store it in a new column
        gdf['area_before_removal'] = gdf.apply(lambda row: row['geometry'].area, axis=1)

        # Check if the specified pkuid is already removed
        if pkuid_to_remove not in gdf['pkuid'].values:
            print(f"Part with pkuid {pkuid_to_remove} is already removed.")
            exit()

        # Display the area before removal in the console with appropriate units
        for index, row in gdf.iterrows():
            area_square_units_before = row['area_before_removal']
            area_acres_before = area_square_units_before * conversion_factor  # Convert area to acres
            print(f"Polygon {index + 1} area before removal: {area_acres_before:.15f} acres")

        # Remove subdivisions based on criteria
        subdivisions_to_remove = gdf.to_dict('records')
        remaining_subdivisions = remove_subdivision(subdivisions_to_remove, pkuid_to_remove)

        # Create a new GeoDataFrame with the remaining subdivisions
        gdf_remaining = gpd.GeoDataFrame.from_records(remaining_subdivisions)

        # Set the CRS for the GeoDataFrame
        gdf_remaining.crs = crs

        # Perform geometric difference to remove the specific parts from the original shapefile polygons
        modified_polygons = []
        areas_after_removal = []  # Store areas after removal in a list
        for index, row in gdf_remaining.iterrows():
            original_polygon = row['geometry']
            new_polygon = remove_polygon_part(original_polygon, pkuid_to_remove)
            if new_polygon:
                modified_polygons.append(new_polygon)
                areas_after_removal.append(new_polygon.area)  # Calculate and store area after removal

        # Create a new GeoDataFrame with the modified polygons after cutting
        gdf_cut = gpd.GeoDataFrame(geometry=modified_polygons)

        # Set the CRS for the GeoDataFrame
        gdf_cut.crs = crs

        # Add the 'area_after_removal' column to the gdf_remaining DataFrame
        column_name = 'area_after_removal'
        i = 1
        while column_name in gdf_remaining.columns:
            i += 1
            column_name = f"area_after_removal_{i}"

        gdf_remaining[column_name] = areas_after_removal

        # Plot the original and modified polygons
        fig, ax = plt.subplots(1, 2, figsize=(12, 4))
        gdf.plot(ax=ax[0], facecolor='none', edgecolor='red', linewidth=1)
        gdf_cut.plot(ax=ax[1], facecolor='none', edgecolor='blue', linewidth=1)
        ax[0].set_title('Original Polygons')
        ax[1].set_title('Modified Polygons')
        plt.tight_layout()
        plt.show()

        # Print the remaining subdivisions in the console, along with areas before and after removal
        print("\nRemaining Subdivisions:\n")
        for index, row in gdf_remaining.iterrows():
            area_square_units_before = row['area_before_removal']
            area_acres_before = area_square_units_before * conversion_factor  # Convert area to acres

            if column_name in row:
                area_square_units_after = row[column_name]
                area_acres_after = area_square_units_after * conversion_factor  # Convert area to acres
                print(f"Polygon {index + 1} area: {area_acres_after:.15f} acres")
            else:
                print(f"Polygon {index + 1} area: {area_acres_before:.15f} acres")
                print("Polygon area after removal: Not available")

        # Update the 'area_after_removal' column after performing the removal operation
        gdf_remaining['area_after_removal'] = gdf_remaining.apply(lambda row: row['geometry'].area, axis=1)


    elif operation == 3:
        # Merge of Polygon
        while True:
            # Ask the user for two pkuid values to merge
            pkuid1 = askinteger("Merge Polygons", "Enter the first pkuid to merge:")
            pkuid2 = askinteger("Merge Polygons", "Enter the second pkuid to merge:")
            
            if pkuid1 is None or pkuid2 is None:
                print("Invalid value. Please enter valid pkuid values.")
            else:
                break
        
        # Perform the merge operation
        new_pkuid = askinteger("Enter the New PKUID","Enter the resultant pkuid : ")
        gdf_remaining = merge_polygons(gdf, pkuid1, pkuid2, new_pkuid)

         # Update the geometry column with the merged polygons
        gdf_remaining['geometry'] = gdf_remaining.apply(lambda row: row['geometry'].buffer(0), axis=1)

        # Update the 'area_after_removal' column with the calculated areas after merge
        gdf_remaining['area_after_removal'] = gdf_remaining.apply(lambda row: row['geometry'].area, axis=1)
        
        # Calculate the area after the merge and store it in a new column
        gdf_remaining = add_area_before_removal(gdf_remaining)
        
        # Plot the original and merged polygons
        fig, ax = plt.subplots(1, 2, figsize=(12, 4))
        gdf.plot(ax=ax[0], facecolor='none', edgecolor='red', linewidth=1)
        gdf_remaining.plot(ax=ax[1], facecolor='none', edgecolor='blue', linewidth=1)
        ax[0].set_title('Original Polygons')
        ax[1].set_title('Merged Polygons')
        plt.tight_layout()
        plt.show()

    # Specify the output file location
    output_shapefile_path = shapefile_path

    if output_shapefile_path:
            # Suppress the warning about column names longer than 10 characters being truncated
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            # Save the GeoDataFrame with the removed parts to a new shapefile
            gdf_remaining.to_file(output_shapefile_path)
        
        print("\nShapefile with appropriate operation saved successfully.")
    else:
        print("Save canceled.")

    #Modification of path file
    if output_shapefile_path:
        # Read the shapefile
        gdf = gpd.read_file(output_shapefile_path)

        # Check if the GeoDataFrame contains multiple polygons
        if len(gdf) < 1 or not all(isinstance(geom, Polygon) for geom in gdf['geometry']):
            print("Invalid shapefile format. Please provide a shapefile containing at least one polygon.")
        else:
            # Convert the polygons to LineStrings
            linestrings_gdf = convert_polygons_to_linestrings(gdf)

            # Visualize the original polygons and the LineStrings
            fig, ax = plt.subplots(figsize=(8, 8))
            gdf.boundary.plot(ax=ax, color='red', linewidth=2, label='Original Polygons')
            linestrings_gdf.plot(ax=ax, color='blue', linewidth=2, linestyle='--', label='LineStrings')
            ax.set_title('Polygons and LineStrings')
            ax.legend()
            plt.show()

            # Prompt the user to choose the location and filename to save the LineStrings shapefile
            save_path = shapefile_path_temp

            if save_path:
                # Save the LineStrings shapefile
                linestrings_gdf.to_file(save_path)
                print("LineStrings shapefile saved successfully.\n")
            else:
                print("Save canceled.")
    else:
        print("No shapefile selected.")

else:
    print("No shapefile selected.")