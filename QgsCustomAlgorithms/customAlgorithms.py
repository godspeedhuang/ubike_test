from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingMultiStepFeedback,
    QgsProcessingParameterNumber,
    QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterFeatureSink,
    QgsCoordinateReferenceSystem,
    QgsRasterLayer,
    QgsRasterBandStats,
)
import processing


class CalculateSlope(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        # Setting the contour interval
        self.addParameter(
            QgsProcessingParameterNumber(
                "contourinterval",    # name
                "contour_interval",    # description
                type=QgsProcessingParameterNumber.Double,
                defaultValue=10,
            ))

        # Get the input DTM raster file
        self.addParameter(
            QgsProcessingParameterRasterLayer("dtm", "DTM", defaultValue=None))

        # Get the grid
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                "grid",
                "grid",
                types=[QgsProcessing.TypeVectorPolygon],
                defaultValue=None,
            ))

        # Get the grid interval
        # TODO: Update this setting
        self.addParameter(
            QgsProcessingParameterNumber(
                "gridinterval",
                "grid_interval",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=250,
            ))

        # Setting the output of this algorithm
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                "Output",
                "Output",
                type=QgsProcessing.TypeVectorAnyGeometry,
                createByDefault=True,
                supportsAppend=True,
                defaultValue=None,
            ))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(12, model_feedback)
        results = {}
        outputs = {}

        # Contour_n
        alg_params = {
            "BAND": 1,
            "CREATE_3D": False,
            "EXTRA": "",
            "FIELD_NAME": "ELEV",
            "IGNORE_NODATA": False,
            "INPUT": parameters["dtm"],
            "INTERVAL": parameters["contourinterval"],
            "NODATA": None,
            "OFFSET": 0,
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["Contour_n"] = processing.run(
            "gdal:contour",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("Contour_n.....Done")

        feedback.setCurrentStep(1)
        if feedback.isCanceled():
            return {}

        # Polygons to lines(grid lines)
        alg_params = {
            "INPUT": parameters["grid"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["PolygonsToLinesgridLines"] = processing.run(
            "native:polygonstolines",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("PolygonsToLinesgridLines.....Done")

        feedback.setCurrentStep(2)
        if feedback.isCanceled():
            return {}

        # Reprojected contour(3826)
        alg_params = {
            "INPUT": outputs["Contour_n"]["OUTPUT"],
            "OPERATION": "",
            "TARGET_CRS": QgsCoordinateReferenceSystem("EPSG:3826"),
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["ReprojectedContour3826"] = processing.run(
            "native:reprojectlayer",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("ReprojectedContour3826.....Done")

        feedback.setCurrentStep(3)
        if feedback.isCanceled():
            return {}

        # Split with lines(overlapped single lines)
        alg_params = {
            "INPUT": outputs["PolygonsToLinesgridLines"]["OUTPUT"],
            "LINES": outputs["PolygonsToLinesgridLines"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["SplitWithLinesoverlappedSingleLines"] = processing.run(
            "native:splitwithlines",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("SplitWithLinesoverlappedSingleLines.....Done")

        feedback.setCurrentStep(4)
        if feedback.isCanceled():
            return {}

        # Delete duplicate geometries(single lines)
        alg_params = {
            "INPUT": outputs["SplitWithLinesoverlappedSingleLines"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["DeleteDuplicateGeometriessingleLines"] = processing.run(
            "native:deleteduplicategeometries",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("DeleteDuplicateGeometriessingleLines.....Done")

        feedback.setCurrentStep(5)
        if feedback.isCanceled():
            return {}

        # Line intersections(grid and contour)
        alg_params = {
            "INPUT":
                outputs["ReprojectedContour3826"]["OUTPUT"],
            "INPUT_FIELDS": [""],
            "INTERSECT":
                outputs["DeleteDuplicateGeometriessingleLines"]["OUTPUT"],
            "INTERSECT_FIELDS": [""],
            "INTERSECT_FIELDS_PREFIX":
                "",
            "OUTPUT":
                QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["LineIntersectionsgridAndContour"] = processing.run(
            "native:lineintersections",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("LineIntersectionsgridAndContour.....Done")

        feedback.setCurrentStep(6)
        if feedback.isCanceled():
            return {}

        # Create spatial index(line intersection)
        alg_params = {
            "INPUT": outputs["LineIntersectionsgridAndContour"]["OUTPUT"]
        }
        outputs["CreateSpatialIndexlineIntersection"] = processing.run(
            "native:createspatialindex",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("CreateSpatialIndexlineIntersection.....Done")

        feedback.setCurrentStep(7)
        if feedback.isCanceled():
            return {}

        # Join attributes by location (points summary)
        alg_params = {
            "DISCARD_NONMATCHING": False,
            "INPUT": parameters["grid"],
            "JOIN": outputs["CreateSpatialIndexlineIntersection"]["OUTPUT"],
            "JOIN_FIELDS": [""],
            "PREDICATE": [0],    # intersects
            "SUMMARIES": [0],    # count
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["JoinAttributesByLocationPointsSummary"] = processing.run(
            "qgis:joinbylocationsummary",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("JoinAttributesByLocationPointsSummary.....Done")

        feedback.setCurrentStep(8)
        if feedback.isCanceled():
            return {}

        # Calculate formula

        contourIntervalNumber = self.parameterAsDouble(parameters,
                                                       "contourinterval",
                                                       context)
        gridIntervelNumber = self.parameterAsDouble(parameters, "gridinterval",
                                                    context)

        # Field calculator(calculate slope)
        alg_params = {
            "FIELD_LENGTH":
                10,
            "FIELD_NAME":
                "Slope",
            "FIELD_PRECISION":
                3,
            "FIELD_TYPE":
                0,    # Float
        # "FORMULA": calculate_slope_expr,
            "FORMULA":
                f'if ("ELEV_count" is null, 0, \
                ("ELEV_count" * 3.14 *  {contourIntervalNumber} )/ (8* {gridIntervelNumber})*100)',
            "INPUT":
                outputs["JoinAttributesByLocationPointsSummary"]["OUTPUT"],
            "OUTPUT":
                QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["FieldCalculatorcalculateSlope"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("FieldCalculatorcalculateSlope.....Done")

        feedback.setCurrentStep(9)
        if feedback.isCanceled():
            return {}

        # Field calculator(calculate slope class)
        alg_params = {
            "FIELD_LENGTH":
                1,
            "FIELD_NAME":
                "class",
            "FIELD_PRECISION":
                0,
            "FIELD_TYPE":
                1,    # Integer
            "FORMULA":
                'CASE WHEN "slope" <= 5 THEN 1\r\n \
                WHEN "slope" <= 15 THEN 2\r\n \
                WHEN "slope" <= 30 THEN 3\r\n \
                WHEN "slope" <= 40 THEN 4\r\n \
                WHEN "slope" <= 55 THEN 5\r\n \
                WHEN "slope" <= 100 THEN 6\r\n \
                WHEN "slope" > 100 THEN 7\r\n \
                    END',
            "INPUT":
                outputs["FieldCalculatorcalculateSlope"]["OUTPUT"],
            "OUTPUT":
                parameters["Output"],
        }
        outputs["FieldCalculatorcalculateSlopeClass"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("FieldCalculatorcalculateSlopeClass.....Done")

        feedback.setCurrentStep(10)
        if feedback.isCanceled():
            return {}

        # Zonal statistics(Mean ELEV)
        alg_params = {
            "COLUMN_PREFIX": "ELEV_",
            "INPUT": parameters["grid"],
            "INPUT_RASTER": parameters["dtm"],
            "RASTER_BAND": 1,
            "STATISTICS": [2],    # Mean
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["ZonalStatisticsmeanElev"] = processing.run(
            "native:zonalstatisticsfb",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        print("ZonalStatisticsmeanElev.....done")

        feedback.setCurrentStep(11)
        if feedback.isCanceled():
            return {}

        # Join attributes by field value(Combination)
        alg_params = {
            "DISCARD_NONMATCHING": False,
            "FIELD": "gridid",
            "FIELDS_TO_COPY": ["ELEV_mean"],
            "FIELD_2": "gridid",
            "INPUT": outputs["FieldCalculatorcalculateSlopeClass"]["OUTPUT"],
            "INPUT_2": outputs["ZonalStatisticsmeanElev"]["OUTPUT"],
            "METHOD":
                1,    # Take attributes of the first matching feature only (one-to-one)
            "PREFIX": "",
            "OUTPUT": parameters["Output"],
        }
        outputs["JoinAttributesByFieldValuecombination"] = processing.run(
            "native:joinattributestable",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        results["Output"] = outputs["JoinAttributesByFieldValuecombination"][
            "OUTPUT"]
        return results

    def name(self):
        return "calculate_slope"

    def displayName(self):
        return "calculate_slope"

    def group(self):
        return "terrain"

    def groupId(self):
        return "terrain"

    def createInstance(self):
        return CalculateSlope()


class VectorizeNDVI(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                "grid",
                "grid",
                types=[QgsProcessing.TypeVectorPolygon],
                defaultValue=None,
            ))
        self.addParameter(
            QgsProcessingParameterRasterLayer("rasterlayer",
                                              "RasterLayer(NDVI_classify)",
                                              defaultValue=None))
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                "Output",
                "Output",
                type=QgsProcessing.TypeVectorAnyGeometry,
                createByDefault=True,
                supportsAppend=True,
                defaultValue=None,
            ))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(9, model_feedback)
        results = {}
        outputs = {}

        # Create spatial index(grid)
        alg_params = {"INPUT": parameters["grid"]}
        outputs["CreateSpatialIndexgrid"] = processing.run(
            "native:createspatialindex",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(1)
        if feedback.isCanceled():
            return {}

        # Warp (reproject)
        alg_params = {
            "DATA_TYPE": 0,    # Use Input Layer Data Type
            "EXTRA": "",
            "INPUT": parameters["rasterlayer"],
            "MULTITHREADING": False,
            "NODATA": None,
            "OPTIONS": "",
            "RESAMPLING": 0,    # Nearest Neighbour
            "SOURCE_CRS": QgsCoordinateReferenceSystem("EPSG:4326"),
            "TARGET_CRS": QgsCoordinateReferenceSystem("EPSG:3826"),
            "TARGET_EXTENT": None,
            "TARGET_EXTENT_CRS": None,
            "TARGET_RESOLUTION": None,
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["WarpReproject"] = processing.run(
            "gdal:warpreproject",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(2)
        if feedback.isCanceled():
            return {}

        # Polygonize (raster to vector)
        alg_params = {
            "BAND": 1,
            "EIGHT_CONNECTEDNESS": False,
            "EXTRA": "",
            "FIELD": "DN",
            "INPUT": outputs["WarpReproject"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["PolygonizeRasterToVector"] = processing.run(
            "gdal:polygonize",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(3)
        if feedback.isCanceled():
            return {}

        # Create spatial index(polygonize)
        alg_params = {"INPUT": outputs["PolygonizeRasterToVector"]["OUTPUT"]}
        outputs["CreateSpatialIndexpolygonize"] = processing.run(
            "native:createspatialindex",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(4)
        if feedback.isCanceled():
            return {}

        # Extract by expression (select tree)
        alg_params = {
            "EXPRESSION": "DN = 2",
            "INPUT": outputs["CreateSpatialIndexpolygonize"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["ExtractByExpressionSelectTree"] = processing.run(
            "native:extractbyexpression",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(5)
        if feedback.isCanceled():
            return {}

        # Intersection(gridNtree)
        alg_params = {
            "INPUT": outputs["CreateSpatialIndexgrid"]["OUTPUT"],
            "INPUT_FIELDS": ["gridid"],
            "OVERLAY": outputs["ExtractByExpressionSelectTree"]["OUTPUT"],
            "OVERLAY_FIELDS": ["DN"],
            "OVERLAY_FIELDS_PREFIX": "",
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["Intersectiongridntree"] = processing.run(
            "native:intersection",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(6)
        if feedback.isCanceled():
            return {}

        # Field calculator(tree_area)
        alg_params = {
            "FIELD_LENGTH": 10,
            "FIELD_NAME": "tree_area",
            "FIELD_PRECISION": 3,
            "FIELD_TYPE": 0,    # Float
            "FORMULA": "$area",
            "INPUT": outputs["Intersectiongridntree"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["FieldCalculatortree_area"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(7)
        if feedback.isCanceled():
            return {}

        # Field calculator(grid_area)
        alg_params = {
            "FIELD_LENGTH": 10,
            "FIELD_NAME": "grid_area",
            "FIELD_PRECISION": 3,
            "FIELD_TYPE": 0,    # Float
            "FORMULA": "62500",
            "INPUT": outputs["FieldCalculatortree_area"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["FieldCalculatorgrid_area"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(8)
        if feedback.isCanceled():
            return {}

        # Field calculator(Coverage)
        alg_params = {
            "FIELD_LENGTH": 10,
            "FIELD_NAME": "coverage",
            "FIELD_PRECISION": 3,
            "FIELD_TYPE": 0,    # Float
            "FORMULA": '"tree_area" / "grid_area"',
            "INPUT": outputs["FieldCalculatorgrid_area"]["OUTPUT"],
            "OUTPUT": parameters["Output"],
        }
        outputs["FieldCalculatorcoverage"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )
        results["Output"] = outputs["FieldCalculatorcoverage"]["OUTPUT"]
        return results

    def name(self):
        return "vectorize_NDVI"

    def displayName(self):
        return "vectorize_NDVI"

    def group(self):
        return "terrain"

    def groupId(self):
        return "terrain"

    def createInstance(self):
        return VectorizeNDVI()


class AggregatePoiToGrid(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterNumber(
                "filternumber",
                "Filter_number",
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=5,
            ))
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                "grid",
                "grid",
                types=[QgsProcessing.TypeVectorPolygon],
                defaultValue=None,
            ))
        self.addParameter(
            QgsProcessingParameterMultipleLayers(
                "pois",
                "POIs",
                layerType=QgsProcessing.TypeVectorPoint,
                defaultValue=None,
            ))
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                "Output",
                "Output",
                type=QgsProcessing.TypeVectorAnyGeometry,
                createByDefault=True,
                supportsAppend=True,
                defaultValue=None,
            ))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(12, model_feedback)
        results = {}
        outputs = {}

        # Create spatial index_grid
        alg_params = {"INPUT": parameters["grid"]}
        outputs["CreateSpatialIndex_grid"] = processing.run(
            "native:createspatialindex",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(1)
        if feedback.isCanceled():
            return {}

        # Merge vector layers
        alg_params = {
            "CRS": QgsCoordinateReferenceSystem("EPSG:4326"),
            "LAYERS": parameters["pois"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["MergeVectorLayers"] = processing.run(
            "native:mergevectorlayers",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(2)
        if feedback.isCanceled():
            return {}

        # Get filter number
        filterNumber = self.parameterAsInt(parameters, "filternumber", context)

        # Extract by expression
        alg_params = {
            "EXPRESSION": f'"rating_num"> {filterNumber} ',
            "INPUT": outputs["MergeVectorLayers"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["ExtractByExpression"] = processing.run(
            "native:extractbyexpression",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(3)
        if feedback.isCanceled():
            return {}

        # Reproject layer
        alg_params = {
            "INPUT": outputs["ExtractByExpression"]["OUTPUT"],
            "OPERATION": "",
            "TARGET_CRS": QgsCoordinateReferenceSystem("EPSG:3826"),
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["ReprojectLayer"] = processing.run(
            "native:reprojectlayer",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(4)
        if feedback.isCanceled():
            return {}

        # Create spatial index_poi
        alg_params = {"INPUT": outputs["ReprojectLayer"]["OUTPUT"]}
        outputs["CreateSpatialIndex_poi"] = processing.run(
            "native:createspatialindex",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(5)
        if feedback.isCanceled():
            return {}

        # Join attributes by location (summary)_rating_star
        alg_params = {
            "DISCARD_NONMATCHING": False,
            "INPUT": outputs["CreateSpatialIndex_grid"]["OUTPUT"],
            "JOIN": outputs["CreateSpatialIndex_poi"]["OUTPUT"],
            "JOIN_FIELDS": ["rating"],
            "PREDICATE": [0],    # intersects
            "SUMMARIES": [5, 6],    # sum,mean
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["JoinAttributesByLocationSummary_rating_star"] = processing.run(
            "qgis:joinbylocationsummary",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(6)
        if feedback.isCanceled():
            return {}

        # Join attributes by location (summary)_rating_num
        alg_params = {
            "DISCARD_NONMATCHING": False,
            "INPUT": outputs["CreateSpatialIndex_grid"]["OUTPUT"],
            "JOIN": outputs["CreateSpatialIndex_poi"]["OUTPUT"],
            "JOIN_FIELDS": ["rating_num"],
            "PREDICATE": [0],    # intersects
            "SUMMARIES": [0, 5],    # count,sum
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["JoinAttributesByLocationSummary_rating_num"] = processing.run(
            "qgis:joinbylocationsummary",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(7)
        if feedback.isCanceled():
            return {}

        # Join attributes by field value
        alg_params = {
            "DISCARD_NONMATCHING":
                False,
            "FIELD":
                "gridid",
            "FIELDS_TO_COPY": ["rating_sum"],
            "FIELD_2":
                "gridid",
            "INPUT":
                outputs["JoinAttributesByLocationSummary_rating_num"]["OUTPUT"],
            "INPUT_2":
                outputs["JoinAttributesByLocationSummary_rating_star"]
                ["OUTPUT"],
            "METHOD":
                1,    # Take attributes of the first matching feature only (one-to-one)
            "PREFIX":
                "",
            "OUTPUT":
                QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["JoinAttributesByFieldValue"] = processing.run(
            "native:joinattributestable",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(8)
        if feedback.isCanceled():
            return {}

        # Field calculator(rating_num_count)
        alg_params = {
            "FIELD_LENGTH": 10,
            "FIELD_NAME": "rating_num_count",
            "FIELD_PRECISION": 0,
            "FIELD_TYPE": 1,    # Integer
            "FORMULA": 'if ("rating_num_count" is null, 0, "rating_num_count")',
            "INPUT": outputs["JoinAttributesByFieldValue"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["FieldCalculatorrating_num_count"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(9)
        if feedback.isCanceled():
            return {}

        # Field calculator(rating_num_sum)
        alg_params = {
            "FIELD_LENGTH": 10,
            "FIELD_NAME": "rating_num_sum",
            "FIELD_PRECISION": 0,
            "FIELD_TYPE": 1,    # Integer
            "FORMULA": 'if ("rating_num_sum" is null, 0, "rating_num_sum")',
            "INPUT": outputs["FieldCalculatorrating_num_count"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["FieldCalculatorrating_num_sum"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(10)
        if feedback.isCanceled():
            return {}

        # Field calculator(rating_sum)
        alg_params = {
            "FIELD_LENGTH": 10,
            "FIELD_NAME": "rating_sum",
            "FIELD_PRECISION": 1,
            "FIELD_TYPE": 0,    # Float
            "FORMULA": 'if ("rating_sum" is null, 0, "rating_sum")',
            "INPUT": outputs["FieldCalculatorrating_num_sum"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["FieldCalculatorrating_sum"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(11)
        if feedback.isCanceled():
            return {}

        # Refactor fields
        alg_params = {
            "FIELDS_MAPPING": [
                {
                    "expression": '"gridid"',
                    "length": 6,
                    "name": "gridid",
                    "precision": 0,
                    "type": 2,
                },
                {
                    "expression": '"rating_num_count"',
                    "length": 10,
                    "name": "rating_num_count",
                    "precision": 0,
                    "type": 2,
                },
                {
                    "expression": '"rating_num_sum"',
                    "length": 10,
                    "name": "rating_num_sum",
                    "precision": 0,
                    "type": 2,
                },
                {
                    "expression": '"rating_sum"',
                    "length": 10,
                    "name": "rating_star_sum",
                    "precision": 1,
                    "type": 6,
                },
            ],
            "INPUT": outputs["FieldCalculatorrating_sum"]["OUTPUT"],
            "OUTPUT": parameters["Output"],
        }
        outputs["RefactorFields"] = processing.run(
            "native:refactorfields",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )
        results["Output"] = outputs["RefactorFields"]["OUTPUT"]
        return results

    def name(self):
        return "aggregate_poi_to_grid"

    def displayName(self):
        return "aggregate_poi_to_grid"

    def group(self):
        return ""

    def groupId(self):
        return ""

    def createInstance(self):
        return AggregatePoiToGrid()


class AggregateStationToGrid(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                "busstation",
                "Bus_station",
                types=[QgsProcessing.TypeVectorPoint],
                defaultValue=None,
            ))
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                "grid",
                "grid",
                types=[QgsProcessing.TypeVectorPolygon],
                defaultValue=None,
            ))
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                "Output",
                "Output",
                type=QgsProcessing.TypeVectorAnyGeometry,
                createByDefault=True,
                supportsAppend=True,
                defaultValue="",
            ))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(7, model_feedback)
        results = {}
        outputs = {}

        # Create spatial index_grid
        alg_params = {"INPUT": parameters["grid"]}
        outputs["CreateSpatialIndex_grid"] = processing.run(
            "native:createspatialindex",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(1)
        if feedback.isCanceled():
            return {}

        # Reproject layer
        alg_params = {
            "INPUT": parameters["busstation"],
            "OPERATION": "",
            "TARGET_CRS": QgsCoordinateReferenceSystem("EPSG:3826"),
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["ReprojectLayer"] = processing.run(
            "native:reprojectlayer",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(2)
        if feedback.isCanceled():
            return {}

        # Create spatial index_station
        alg_params = {"INPUT": outputs["ReprojectLayer"]["OUTPUT"]}
        outputs["CreateSpatialIndex_station"] = processing.run(
            "native:createspatialindex",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(3)
        if feedback.isCanceled():
            return {}

        # Join attributes by location (summary)_RouteNum
        alg_params = {
            "DISCARD_NONMATCHING": False,
            "INPUT": outputs["CreateSpatialIndex_grid"]["OUTPUT"],
            "JOIN": outputs["CreateSpatialIndex_station"]["OUTPUT"],
            "JOIN_FIELDS": ["RouteNum"],
            "PREDICATE": [0],    # intersects
            "SUMMARIES": [0, 5],    # count,sum
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["JoinAttributesByLocationSummary_routenum"] = processing.run(
            "qgis:joinbylocationsummary",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(4)
        if feedback.isCanceled():
            return {}

        # Field calculator(RouteNum_count)
        alg_params = {
            "FIELD_LENGTH":
                10,
            "FIELD_NAME":
                "RouteNum_count",
            "FIELD_PRECISION":
                0,
            "FIELD_TYPE":
                1,    # Integer
            "FORMULA":
                'if ("RouteNum_count" is null, 0, "RouteNum_count")',
            "INPUT":
                outputs["JoinAttributesByLocationSummary_routenum"]["OUTPUT"],
            "OUTPUT":
                QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["FieldCalculatorroutenum_count"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(5)
        if feedback.isCanceled():
            return {}

        # Field calculator(RouteNum_sum)
        alg_params = {
            "FIELD_LENGTH": 10,
            "FIELD_NAME": "RouteNum_sum",
            "FIELD_PRECISION": 0,
            "FIELD_TYPE": 1,    # Integer
            "FORMULA": 'if ("RouteNum_sum" is null, 0, "RouteNum_sum")',
            "INPUT": outputs["FieldCalculatorroutenum_count"]["OUTPUT"],
            "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
        }
        outputs["FieldCalculatorroutenum_sum"] = processing.run(
            "native:fieldcalculator",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.setCurrentStep(6)
        if feedback.isCanceled():
            return {}

        # Refactor fields
        alg_params = {
            "FIELDS_MAPPING": [
                {
                    "expression": '"gridid"',
                    "length": 6,
                    "name": "gridid",
                    "precision": 0,
                    "type": 2,
                },
                {
                    "expression": '"RouteNum_count"',
                    "length": 10,
                    "name": "RouteNum_count",
                    "precision": 0,
                    "type": 2,
                },
                {
                    "expression": '"RouteNum_sum"',
                    "length": 10,
                    "name": "RouteNum_sum",
                    "precision": 0,
                    "type": 2,
                },
            ],
            "INPUT": outputs["FieldCalculatorroutenum_sum"]["OUTPUT"],
            "OUTPUT": parameters["Output"],
        }
        outputs["RefactorFields"] = processing.run(
            "native:refactorfields",
            alg_params,
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )
        results["Output"] = outputs["RefactorFields"]["OUTPUT"]
        return results

    def name(self):
        return "aggregate_station_to_grid"

    def displayName(self):
        return "aggregate_station_to_grid"

    def group(self):
        return "Bus"

    def groupId(self):
        return "Bus"

    def createInstance(self):
        return AggregateStationToGrid()
