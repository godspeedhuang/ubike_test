from qgis.core import (
    QgsApplication,
    QgsVectorLayer,
    QgsRasterLayer,
)
from osgeo import gdal
import glob
import processing
from processing.core.Processing import Processing

from QgsCustomAlgorithms.customProvider import CustomProvider

# Initialize QGIS
qgs = QgsApplication([], False)
qgs.initQgis()

# Initialize processing plugin
Processing.initialize()  # TODO: 會自動生成processing folder，再處理一下

# Add the custom provider and custom algorithms
provider = CustomProvider()
QgsApplication.processingRegistry().addProvider(provider)


# Model input
# dtm = QgsRasterLayer(r"dataset\DTMerge.tif", "DTMerge", "gdal")
grid_layer = QgsVectorLayer(r"dataset\FET_2023_grid_97.geojson", "grid", "ogr")

# read dtm tiles
DTM_input_files = [file for file in glob.glob("dataset/分幅_臺北市20MDEM/*.grd")]
processing_options = gdal.WarpOptions(
    format="VRT",  # Gdal virtual format
    dstSRS="EPSG:3826",
    dstNodata=0,
)

# merge dtm tiles
in_memory_vrt_path = "/vsimem/merged.vrt"  # set memory path
g = gdal.Warp(in_memory_vrt_path, DTM_input_files, options=processing_options)
# Release the GDAL dataset object to free up memory and resources

DTM_layer = QgsRasterLayer(in_memory_vrt_path, "Merged DTM", "gdal")

params = {
    "contourinterval": 10,
    "dtm": DTM_layer,
    "grid": grid_layer,
    "gridinterval": 250,
    "Output": "TEMPORARY_OUTPUT",
}

result = processing.run("customProvider:calculate_slope", params)
print("RESULT:", result)
# print("output:", output)
gdal.Unlink(
    in_memory_vrt_path
)  # Release the GDAL dataset object to free up memory and resources
# qgs.exitQgis()
