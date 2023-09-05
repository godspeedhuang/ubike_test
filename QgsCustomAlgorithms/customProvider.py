from qgis.core import QgsProcessingProvider
from .customAlgorithms import (
    CalculateSlope,
    VectorizeNDVI,
    AggregatePoiToGrid,
    AggregateStationToGrid,
)


class CustomProvider(QgsProcessingProvider):

    def loadAlgorithms(self, *args, **kwargs):
        self.addAlgorithm(CalculateSlope())
        self.addAlgorithm(VectorizeNDVI())
        self.addAlgorithm(AggregatePoiToGrid())
        self.addAlgorithm(AggregateStationToGrid())

    def id(self, *args, **kwargs):
        """Used for identifying the provider.

        This string should be a unique, short, character only string,
        eg "qgis" or "gdal". This string should not be localised.
        """
        return "customProvider"

    def name(self, *args, **kwargs):
        """The human friendly name of your plugin in Processing.

        This string should be as short as possible (e.g. "Lastools", not
        "Lastools version 1.0.1 64-bit") and localised.
        """
        return self.tr("customProvider")
