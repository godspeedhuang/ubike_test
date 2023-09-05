import pandas as pd
from PyQt5.QtCore import QVariant
from qgis.core import (
    QgsVectorLayer,
    QgsFields,
    QgsField,
    QgsCoordinateReferenceSystem,
    QgsPointXY,
    QgsGeometry,
    QgsFeature,
    QgsWkbTypes,
    QgsMemoryProviderUtils,
)


def QgsVectorLayer_to_dataframe(vector_layer: QgsVectorLayer,
                                target_fields: list) -> pd.DataFrame:
    """
    Extract designated columns from QgsVector Layer
    and convert datatype into dataframe
    """
    data_container = []
    for feature in vector_layer.getFeatures():
        values = [feature[target_field] for target_field in target_fields]
        data_container.append(values)

    df = pd.DataFrame(data_container, columns=target_fields)
    return df


# Read csv
def csv_to_QgsVectorLayer(
    fp: str,
    layer_name: str = "csv_to_QgsVectorLayer",
    lat: float = "lat",
    lng: float = "lng",
    **kwargs: dict[str, str],
) -> QgsVectorLayer:
    """
    Conver csv(Point) to QgsVectorLayer in order to fit in custom algorithms (As a valid parameter)
    """
    # Read csv
    df = pd.read_csv(fp).fillna(0)

    # Define Fields
    fields = QgsFields()
    for field_name, field_type in kwargs.items():
        if field_type == "Double":
            field = QgsField(field_name, QVariant.Double)
        elif field_type == "Int":
            field = QgsField(field_name, QVariant.Int)
        else:
            raise ValueError(f"Unsupported field type: {field_type}")
        fields.append(field)

    # Create Empty Memory Layer
    vector_layer = QgsMemoryProviderUtils.createMemoryLayer(
        layer_name,
        fields,
        QgsWkbTypes.Point,
        QgsCoordinateReferenceSystem("EPSG:4326"),
    )

    def add_features(row):
        """Write features into Memory Layer"""
        # Add attribute to the feature(row)
        feature = QgsFeature(fields)
        for key in kwargs.keys():
            feature.setAttribute(key, row[key])

        # Set geometry info
        geom = QgsGeometry.fromPointXY(QgsPointXY(row[lng], row[lat]))
        feature.setGeometry(geom)

        # Append one completed feature(row) into the layer
        vector_layer.dataProvider().addFeature(feature)

    df.apply(add_features, axis=1)

    return vector_layer
