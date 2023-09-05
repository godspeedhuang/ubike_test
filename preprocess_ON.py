import datetime
import glob
import itertools
# import xmltodict
# import pykml.parser
import json
import os
import pickle
import re
# from tqdm.notebook import tqdm
import warnings

import constants
# import twd97
import geopandas as gpd
import numpy as np
import pandas as pd
import processing
import rasterio as rio
import rasterio.features
import rasterstats
# Geoprocessing
from osgeo import gdal
from processing.core.Processing import Processing
from qgis.core import QgsApplication, QgsRasterLayer, QgsVectorLayer
# Import custom algorithms
from QgsCustomAlgorithms.customProvider import CustomProvider
from QgsCustomAlgorithms.generalFunction import QgsVectorLayer_to_dataframe, csv_to_QgsVectorLayer
from rasterio.warp import calculate_default_transform
from shapely.geometry import LineString, Point, Polygon, shape
from shapely.geometry.multipolygon import MultiPolygon

pd.set_option("display.max_columns", None)
warnings.filterwarnings("ignore")


class Preprocess:

    def __init__(self, DateList):
        # date
        self.DateList = DateList
        self.HyphenDateList = [
            x[0:4] + "-" + x[4:6] + "-" + x[6:8] for x in self.DateList
        ]

        # file path
        # self.TransactionFilePath = "input/data/transaction"
        self.PopulationFilePath = "input/data/population"
        self.TrafficFilePath = "input/data/traffic"
        self.RoadFilePath = "input/data/road"
        self.TerrainFilePath = "input/data/terrain"
        self.DevelopmentFilePath = "input/data/development"
        self.LandFilePath = "input/data/land"
        self.PoiFilePath = "input/data/POI"
        self.DoneFilePath = "output_ON"

        self.gridFilePath = f"{self.PopulationFilePath}/FET_2023_grid_97.geojson"
        # grid information
        self.gridDF = gpd.read_file(f"{self.gridFilePath}").set_crs("epsg:3826")

    def Point_WGS84toTWD97(self, DF, colName, reverse):
        # get all points' wgs84 coordinate
        # e.g. "POINT (121.56654 25.08235)", "POINT (121.50449 25.10045)", ...

        wgs84List = DF[colName].astype(str).str[6:].str.replace(" ",
                                                                ", ").tolist()
        wgs84List = [eval(x) for x in wgs84List]

        # WGS84: lng, lat
        if reverse:
            twd97List = [twd97.fromwgs84(x[1], x[0]) for x in wgs84List]

        # WGS84: lat, lng
        else:
            twd97List = [twd97.fromwgs84(x[0], x[1]) for x in wgs84List]

        # output
        DF["twd97_lat"] = [x[0] for x in twd97List]
        DF["twd97_lng"] = [x[1] for x in twd97List]

        return DF

    # get the intersection area within two Polygon series
    def OverlayWithinPPL(self, df1, df2, key1, key2, method):
        # df1: "polygon df" or "point df" or "line df"
        # df2: "grid df"
        # method: "Polygon" or "Point" or "Line"

        overlayDF = gpd.overlay(df1, df2, how="union").explode().reset_index()
        overlayDF = overlayDF[[key1, key2, "geometry"]]
        overlayDF = overlayDF[overlayDF[key1].notnull()]
        overlayDF = overlayDF[overlayDF[key2].notnull()]
        overlayDF[key2] = overlayDF[key2].astype(int)
        overlayDF.index = range(len((overlayDF)))

        if method == "Polygon":
            overlayDF["area"] = overlayDF.geometry.area
        elif method == "Line":
            overlayDF["length"] = overlayDF.geometry.length
        elif method == "Point":
            pass

        return overlayDF

    # transaction
    def Transaction(self):
        print("Preprocessing for Transaction Data...")

        # read stop data
        StopDF = pd.read_csv(
            f"{self.TransactionFilePath}/86ec099baa2d36c22ab3a87350b718de_export.csv"
        )
        StopDF = StopDF[["sno", "lat", "lng"]]
        StopDF["sno"] = StopDF["sno"].astype(str)
        StopDF["sno"] = "U" + StopDF["sno"].str[3:]
        StopDF["geometry"] = [
            Point(twd97.fromwgs84(x, y))
            for x, y in zip(StopDF.lat, StopDF.lng)
        ]

        # read transaction data
        DF = pd.DataFrame()
        for date in self.DateList:
            tem = open(
                f"{self.TransactionFilePath}/202303_txn_identified_transfer/{date}.pkl",
                "rb",
            )
            temFile = pickle.load(tem)
            DF = pd.concat([DF, temFile])
            tem.close()
        DF.index = range(0, len(DF))

        # convert on_time(datetime), off_time(datetime) to on_hour(int), off_hour(int)
        # hardcore the location of date and hour
        DF["on_date"] = DF["on_time"].astype(str).str[:10].str.split(
            "-").str.join("")
        DF["off_date"] = DF["off_time"].astype(str).str[:10].str.split(
            "-").str.join("")
        DF["on_hour"] = DF["on_time"].astype(str).str[11:13].astype(int)
        DF["off_hour"] = DF["off_time"].astype(str).str[11:13].astype(int)

        # complete mapping table
        ON_STOP_ID = list(set(DF["on_stop_id"]))
        ONDATE, OFFDATE = list(set(DF["on_date"])), list(set(DF["off_date"]))
        ONHOUR, OFFHOUR = list(set(DF["on_hour"])), list(set(DF["off_hour"]))
        onColumns, offColumns = ["on_stop_id", "on_date", "on_hour"], [
            "off_stop_id",
            "off_date",
            "off_hour",
        ]
        productElements = list(itertools.product(ON_STOP_ID, ONDATE, ONHOUR))
        OnMappingDF = pd.DataFrame(productElements, columns=onColumns)
        OffMappingDF = pd.DataFrame(productElements, columns=offColumns)

        # on (借車)
        onDF = DF.groupby(onColumns).size().reset_index(name="counts")
        onDF = onDF.merge(OnMappingDF, how="right", on=onColumns)
        onDF.fillna(0, inplace=True)

        # off (還車)
        offDF = DF.groupby(offColumns).size().reset_index(name="counts")
        offDF = offDF.merge(OffMappingDF, how="right", on=offColumns)
        offDF.fillna(0, inplace=True)

        # stop_id's information include lat, lng
        onDF = onDF.merge(StopDF,
                          how="left",
                          left_on="on_stop_id",
                          right_on="sno").drop(["sno"], axis=1)
        offDF = offDF.merge(StopDF,
                            how="left",
                            left_on="off_stop_id",
                            right_on="sno").drop(["sno"], axis=1)

        # get "on" geometry, "off" geometry and Grid geometry
        On_point = gpd.GeoDataFrame(onDF[["on_stop_id", "geometry"]])
        Off_point = gpd.GeoDataFrame(offDF[["off_stop_id", "geometry"]])
        Grid_poly = self.gridDF[["gridid", "geometry"]]

        # drop duplicate
        On_point.drop_duplicates(inplace=True)
        Off_point.drop_duplicates(inplace=True)
        On_point.index = range(len(On_point))
        Off_point.index = range(len(Off_point))

        # overlay
        OverlayOn_DF = self.OverlayWithinPPL(On_point,
                                             Grid_poly,
                                             "on_stop_id",
                                             "gridid",
                                             method="Point")
        OverlayOff_DF = self.OverlayWithinPPL(Off_point,
                                              Grid_poly,
                                              "off_stop_id",
                                              "gridid",
                                              method="Point")

        # mapping
        onDF = onDF.merge(OverlayOn_DF.drop(["geometry"], axis=1),
                          how="left",
                          on="on_stop_id")
        offDF = offDF.merge(OverlayOff_DF.drop(["geometry"], axis=1),
                            how="left",
                            on="off_stop_id")

        # kick out NA, reset index and astype
        onDF = onDF[onDF.gridid.notnull()]
        offDF = offDF[offDF.gridid.notnull()]
        onDF.index = range(len(onDF))
        offDF.index = range(len(offDF))
        onDF[["counts", "gridid"]] = onDF[["counts", "gridid"]].astype(int)
        offDF[["counts", "gridid"]] = offDF[["counts", "gridid"]].astype(int)

        print("Done!")

        return onDF, offDF

    # population (人口信令)
    def Population(self):
        print("Preprocessing for Population Data...")

        # read data
        populationDF = pd.read_csv(
            f"{self.PopulationFilePath}/台北市停留人口_資料集_1.csv")
        workDF = pd.read_csv(
            f"{self.PopulationFilePath}/台北市停留人口_資料集_2_工作人口.csv")
        liveDF = pd.read_csv(
            f"{self.PopulationFilePath}/台北市停留人口_資料集_2_居住人口.csv")
        tourDF = pd.read_csv(
            f"{self.PopulationFilePath}/台北市停留人口_資料集_2_遊客人口.csv")
        populationDF["日期"] = populationDF["日期"].str.split("-").str.join("")
        workDF["日期"] = workDF["日期"].str.split("-").str.join("")
        liveDF["日期"] = liveDF["日期"].str.split("-").str.join("")
        tourDF["日期"] = tourDF["日期"].str.split("-").str.join("")

        # by age and total
        DF = (populationDF.groupby(["日期", "時間", "網格編號",
                                    "年齡別"]).sum("放大後人數").reset_index())
        DF = (pd.pivot_table(DF,
                             values="放大後人數",
                             index=["日期", "時間", "網格編號"],
                             columns="年齡別").reset_index().fillna(0))
        DF.columns = [
            "日期",
            "時間",
            "網格編號",
            "Age_15_17_Counts",
            "Age_18_21_Counts",
            "Age_22_29_Counts",
            "Age_30_39_Counts",
            "Age_40_49_Counts",
            "Age_50_59_Counts",
            "Age_60_64_Counts",
            "Age_Over65_Counts",
        ]
        DF["Age_Total_Counts"] = (
            DF["Age_15_17_Counts"] + DF["Age_18_21_Counts"] +
            DF["Age_22_29_Counts"] + DF["Age_30_39_Counts"] +
            DF["Age_40_49_Counts"] + DF["Age_50_59_Counts"] +
            DF["Age_60_64_Counts"] + DF["Age_Over65_Counts"])

        mergeCol = ["日期", "時間", "網格編號"]

        # by work
        DF = DF.merge(workDF, how="left", on=mergeCol)
        DF.rename(columns={"放大後人數": "WorkPopulationCounts"}, inplace=True)

        # by live
        DF = DF.merge(liveDF, how="left", on=mergeCol)
        DF.rename(columns={"放大後人數": "LivePopulationCounts"}, inplace=True)

        # by tour
        DF = DF.merge(tourDF, how="left", on=mergeCol)
        DF.rename(columns={"放大後人數": "TourPopulationCounts"}, inplace=True)

        print("Done!")

        return DF

    # traffic for MRT
    def TrafficMRT(self):
        print("Preprocessing for MRT Data...")

        # read data
        MRT_DF = pd.read_csv(f"{self.TrafficFilePath}/臺北捷運車站出入口座標.csv",
                             encoding="big5")
        population = pd.read_csv(
            f"{self.TrafficFilePath}/臺北捷運每日分時各站OD流量統計資料_202303.csv")

        ## information of MRT Station and Exit
        # generate twd97 coordinate
        MRT_DF["Exit_twd97"] = MRT_DF.apply(
            lambda x: twd97.fromwgs84(x.緯度, x.經度), axis=1)
        MRT_DF[["Exit_twd97_lat",
                "Exit_twd97_lng"]] = MRT_DF["Exit_twd97"].tolist()

        # get Station Name
        MRT_DF["Station"] = MRT_DF["出入口名稱"].str.split("站出口").str[0]

        # split MRT_DF into Exit_DF and Station_DF
        Exit_DF = MRT_DF.drop(columns=["項次", "Station"])
        Station_DF = (MRT_DF[["Exit_twd97_lat", "Exit_twd97_lng",
                              "Station"]].groupby("Station").agg(
                                  Station_twd97_lat=("Exit_twd97_lat", "mean"),
                                  Station_twd97_lng=("Exit_twd97_lng", "mean"),
                              ).reset_index())

        # create column 'geometry'
        Exit_DF["geometry"] = Exit_DF.Exit_twd97.apply(Point)
        Station_DF["geometry"] = [
            Point(x) for x in list(
                zip(Station_DF.Station_twd97_lat, Station_DF.Station_twd97_lng))
        ]

        # get exit geometry, station geometry and Grid geometry
        Exit_point = gpd.GeoDataFrame(Exit_DF[["出入口名稱", "geometry"]])
        Station_point = gpd.GeoDataFrame(Station_DF[["Station", "geometry"]])
        Grid_poly = self.gridDF[["gridid", "geometry"]]

        # overlay
        OverlayExit_DF = self.OverlayWithinPPL(Exit_point,
                                               Grid_poly,
                                               "出入口名稱",
                                               "gridid",
                                               method="Point")
        OverlayStation_DF = self.OverlayWithinPPL(Station_point,
                                                  Grid_poly,
                                                  "Station",
                                                  "gridid",
                                                  method="Point")

        ## population of MRT Stop (各站點進出人次)
        population = population[population["日期"].isin(self.HyphenDateList)]
        populationIN = population[["日期", "時段", "進站", "人次"]]
        populationOUT = population[["日期", "時段", "出站", "人次"]]
        populationIN.index = range(0, len(populationIN))
        populationOUT.index = range(0, len(populationOUT))

        # calculate number of people by day, hour, station
        populationIN = populationIN.groupby(["日期", "時段",
                                             "進站"]).sum("人次").reset_index()
        populationOUT = (populationOUT.groupby(["日期", "時段",
                                                "出站"]).sum("人次").reset_index())

        # get GRID ID
        populationIN = populationIN.merge(OverlayStation_DF,
                                          how="inner",
                                          left_on="進站",
                                          right_on="Station")
        populationOUT = populationOUT.merge(OverlayStation_DF,
                                            how="inner",
                                            left_on="出站",
                                            right_on="Station")

        print("Done!")

        return OverlayExit_DF, populationIN, populationOUT

    # traffic for BUS
    def TrafficBus(self):
        print("Preprocessing for Bus Data...")

        # read data
        DF = gpd.read_file(f"{self.TrafficFilePath}/busstop/busstop.shp")

        # get TWD97 coordinate
        DF = self.Point_WGS84toTWD97(DF, "geometry", reverse=True)
        DF.drop(columns=["geometry"], inplace=True)

        # create twd97 geometry series
        DF["LatLng"] = list(zip(DF.twd97_lat, DF.twd97_lng))
        DF["geometry"] = DF.LatLng.apply(Point)

        # get bus geometry and Grid geometry
        Bus_point = DF[["BSM_BUSSTO", "geometry"]]
        Grid_poly = self.gridDF[["gridid", "geometry"]]

        # overlay
        overlayDF = self.OverlayWithinPPL(Bus_point, Grid_poly, "BSM_BUSSTO",
                                          "gridid", "Point")

        print("Done!")

        return overlayDF

    # road for bus route
    def TrafficBusRoute(self):
        print("Preprocessing for BusRoute Data...")

        bus_df = pd.read_csv(f"{self.TrafficFilePath}/bus_station_detailed.csv")
        bus_gdf = (gpd.GeoDataFrame(
            bus_df[["RouteNum"]],
            geometry=gpd.points_from_xy(bus_df["lon"], bus_df["lat"]),
        ).set_crs("epsg:4326").to_crs("epsg:3826"))
        bus_agg_df = ((gpd.sjoin(
            self.gridDF[["gridid", "geometry"]],
            bus_gdf,
            predicate="intersects",
            how="inner",
        )).groupby(["gridid"]).agg({
            "RouteNum": "sum"
        }).rename(columns={
            "sum": "BusRouteCounts"
        }).reset_index())

        bus_merge_df = pd.merge(self.gridDF[["gridid"]],
                                bus_agg_df,
                                how="left",
                                on="gridid").fillna(0)
        # # read data
        # column_definition = {
        #     "RouteNum": "Int",
        # }

        # bus_station_points = csv_to_QgsVectorLayer(
        #     f"{self.TrafficFilePath}/bus_station_detailed.csv",  # csv filepath
        #     "bus_station_points",
        #     "lat",
        #     "lon",
        #     **column_definition,
        # )

        # params = {
        #     "busstation": bus_station_points,
        #     "grid": self.gridFilePath,
        #     "Output": "TEMPORARY_OUTPUT",
        # }

        # result = processing.run("customProvider:aggregate_station_to_grid", params)

        # BusRouteDF = QgsVectorLayer_to_dataframe(
        #     result["Output"], ["gridid", "RouteNum_sum"]
        # )
        # BusRouteDF.rename(columns={"RouteNum_sum": "BusRouteCounts"}, inplace=True)

        print("Done!")
        print(bus_merge_df)
        return bus_merge_df

    # road for side walk
    def RoadSideWalk(self):
        print("Preprocessing for SideWalk Data...")

        # read SideWalk data
        f = open(f"{self.RoadFilePath}/TP_SIDEWORK.json")
        sidewalkDF = json.load(f)
        sidewalkDF = pd.json_normalize(sidewalkDF["features"])
        sidewalkDF = gpd.GeoDataFrame(sidewalkDF)
        sidewalkDF["geometry"] = [
            Polygon(x[0][0]) for x in sidewalkDF["geometry.coordinates"]
        ]

        # get sidewalk geometry and Grid geometry
        SW_poly = sidewalkDF[["properties.ObjectID", "geometry"]]
        Grid_poly = self.gridDF[["gridid", "geometry"]]

        # overlay
        overlayDF = self.OverlayWithinPPL(SW_poly,
                                          Grid_poly,
                                          "properties.ObjectID",
                                          "gridid",
                                          method="Polygon")

        print("Done!")

        return overlayDF

    # road for Marked SideWalk data
    def RoadMarkedSideWalk(self):
        print("Preprocessing for Marked SideWalk Data...")

        # read Marked SideWalk data
        MarkedSideWalkDF = gpd.read_file(
            f"{self.RoadFilePath}/(交工處)標線型人行道圖資_202304171730/grapline_21_15.shp"
        )

        # get sidewalk geometry and Grid geometry
        MSW_poly = MarkedSideWalkDF[["KEYID", "geometry"]]
        Grid_poly = self.gridDF[["gridid", "geometry"]]

        # overlay
        overlayDF = self.OverlayWithinPPL(MSW_poly,
                                          Grid_poly,
                                          "KEYID",
                                          "gridid",
                                          method="Polygon")

        print("Done!")

        return overlayDF

    # road for Bike
    def RoadBikeRoute(self):
        print("Preprocessing for Bike Route Data...")

        # read data
        root = pykml.parser.fromstring(
            open(f"{self.RoadFilePath}/台北市_自行車道-市區自行車道1120505.kml", "r").read())
        # print(root.Document.Folder.Placemark.name)
        # print(root.Document.Folder.Placemark.LineString.coordinates)
        Name = [
            x.name.text for x in root.findall(
                ".//{http://www.opengis.net/kml/2.2}Placemark")
        ]
        Coor = [
            str(x.LineString.coordinates).split(",0 ")[0:-1] for x in
            root.findall(".//{http://www.opengis.net/kml/2.2}Placemark")
        ]

        # get all twd97 coordinates of road of bike
        twd97Coor = list()
        for x in Coor:
            tem = []
            for y in x:
                e = y.split(",")
                t = twd97.fromwgs84(e[1], e[0])
                tem.append(t)
            twd97Coor.append(tem)

        DF = gpd.GeoDataFrame(data={"Name": Name, "twd97Coor": twd97Coor})

        # remove some road of bike, which's len(coordinate) is 1
        removeIndex = list()
        for x in range(len(DF)):
            if len(DF.loc[x, "twd97Coor"]) == 1:
                removeIndex.append(x)
        DF.drop(index=removeIndex, inplace=True)

        DF["geometry"] = DF.twd97Coor.apply(LineString)

        # get road of bike DataFrame and Grid DataFrame
        gdf_lines = DF[["Name", "geometry"]]
        gdf_poly = self.gridDF[["gridid", "geometry"]]

        # get intersection LineString between all road of bike and all Grid
        OverlayDF = (gpd.overlay(gdf_lines, gdf_poly,
                                 how="union").explode().reset_index(drop=True))

        # get Taipei DF
        TpeOverlayDF = OverlayDF[OverlayDF["gridid"].notnull()]
        TpeOverlayDF.index = range(len(TpeOverlayDF))
        TpeOverlayDF["gridid"] = TpeOverlayDF["gridid"].astype(int)

        # calaculate Length of road of bike
        TpeOverlayDF["length"] = TpeOverlayDF.geometry.length

        print("Done!")

        return TpeOverlayDF

    # road for tree
    def RoadTree(self):
        print("Preprocessing for Tree Data...")

        # read data
        f = open(f"{self.RoadFilePath}/TaipeiTree.json")
        TreeDF = json.load(f)
        TreeDF = pd.json_normalize(TreeDF)
        TreeDF["LatLng"] = list(zip(TreeDF.X, TreeDF.Y))
        TreeDF = gpd.GeoDataFrame(TreeDF)
        TreeDF["geometry"] = TreeDF.LatLng.apply(Point)

        # get Tree geometry and Grid geometry
        Tree_point = TreeDF[["TreeID", "geometry"]]
        Grid_poly = self.gridDF[["gridid", "geometry"]]

        # overlay
        overlayDF = self.OverlayWithinPPL(Tree_point,
                                          Grid_poly,
                                          "TreeID",
                                          "gridid",
                                          method="Point")

        print("Done!")

        return overlayDF

    # road for light
    def RoadLight(self):
        print("Preprocessing for Light Data...")

        # read data
        f = open(f"{self.RoadFilePath}/TaipeiLight.json")
        LightDF = json.load(f)
        LightDF = pd.json_normalize(LightDF)
        LightDF = LightDF[LightDF.X != ""]
        LightDF = LightDF[LightDF.Y != ""]
        LightDF.index = range(len(LightDF))
        LightDF["LatLng"] = list(zip(LightDF.X, LightDF.Y))
        LightDF = gpd.GeoDataFrame(LightDF)
        LightDF["geometry"] = LightDF.LatLng.apply(Point)

        # get Tree geometry and Grid geometry
        Light_point = LightDF[["LIGHTID", "geometry"]]
        Grid_poly = self.gridDF[["gridid", "geometry"]]

        # overlay
        overlayDF = self.OverlayWithinPPL(Light_point,
                                          Grid_poly,
                                          "LIGHTID",
                                          "gridid",
                                          method="Point")

        print("Done!")

        return overlayDF

    # road for road net
    def RoadRoadNet(self):
        print("Preprocessing for RoadNet Data...")

        # read data
        road_polygon = gpd.read_file(os.path.join(self.RoadFilePath,
                                                  "Road.shp"))

        # Overlay grid and road_polygon
        OverlayRoadArea_DF = self.OverlayWithinPPL(
            road_polygon, self.gridDF[["gridid", "geometry"]], "ID", "gridid",
            "Polygon")

        # Sum the road area of each grid
        grid_road_area = OverlayRoadArea_DF.groupby("gridid")["area"].sum()

        RoadNetDF = (self.gridDF[["gridid",
                                  "Area"]].merge(grid_road_area,
                                                 how="left",
                                                 on="gridid").fillna(0))
        RoadNetDF["ratio"] = RoadNetDF["area"] / RoadNetDF["Area"]
        RoadNetDF = RoadNetDF.round(3)
        RoadNetDF.rename(columns={
            "area": "RoadArea",
            "ratio": "RoadAreaRatio"
        },
                         inplace=True)
        print("Done!")

    # road for road length
    def RoadRoadLength(self):
        print("Preprocessing for Road Length Data...")

        # read data
        road_line = gpd.read_file(os.path.join(self.RoadFilePath,
                                               "MidRoad.shp"))

        road_lengths = {}
        for category, road_type_list in constants.ROAD_CATEGORY.items():
            selected_roads = road_line[road_line["ROADTYPE"].isin(
                road_type_list)]
            road_lengths[category] = (
                self.OverlayWithinPPL(
                    selected_roads[["ROADTYPE", "geometry"]],  # selected roads
                    self.gridDF[["gridid", "geometry"]],
                    "ROADTYPE",
                    "gridid",
                    "Line",  # data type
                ).groupby("gridid")["length"].sum()
            )  # return road length by each grid

        RoadLengthDF = self.gridDF[["gridid"]].copy()

        for idx, length_by_grid in enumerate(road_lengths.values()):
            RoadLengthDF = RoadLengthDF.merge(length_by_grid,
                                              how="left",
                                              on="gridid",
                                              suffixes=(None,
                                                        f"_{idx}")).fillna(0)

        RoadLengthDF = RoadLengthDF.rename(
            columns={
                RoadLengthDF.columns[1]:
                    list(constants.ROAD_CATEGORY.keys())[0],
                RoadLengthDF.columns[2]:
                    list(constants.ROAD_CATEGORY.keys())[1],
                RoadLengthDF.columns[3]:
                    list(constants.ROAD_CATEGORY.keys())[2],
                RoadLengthDF.columns[4]:
                    list(constants.ROAD_CATEGORY.keys())[3],
            })

        RoadLengthDF["TotalRoadLength"] = RoadLengthDF.iloc[:, 1:].sum(axis=1)
        RoadLengthDF = RoadLengthDF.round(3)
        print("Done!")

        return RoadLengthDF

    # Terrain
    def Terrain(self):
        print("Preprocessing for Terrain Data...")

        # read dtm tiles
        DTM_input_files = [
            file
            for file in glob.glob(f"{self.TerrainFilePath}/分幅_臺北市20MDEM/*.grd")
        ]
        processing_options = gdal.WarpOptions(
            format="GTiff",
            dstSRS="EPSG:3826",
            dstNodata=0,
        )

        # merge dtm tiles
        outputPath = os.path.join(self.TerrainFilePath,
                                  "DTMerge.tif")  # set output filepath
        g = gdal.Warp(outputPath, DTM_input_files, options=processing_options)
        g = None  # Release the GDAL dataset object to free up memory and resources

        # Initialize input data
        DTM_layer = QgsRasterLayer(outputPath, "Merged DTM", "gdal")
        grid_layer = QgsVectorLayer(
            self.gridFilePath,
            "grid",
            "ogr",
        )
        params = {
            "contourinterval": 10,
            "dtm": DTM_layer,
            "grid": grid_layer,
            "gridinterval": 250,
            "Output": "TEMPORARY_OUTPUT",
        }
        dtm_result = processing.run("customProvider:calculate_slope", params)
        target_fields = ["gridid", "Slope", "class", "ELEV_mean"]
        dtm_df = QgsVectorLayer_to_dataframe(dtm_result["Output"], target_fields)

        dtm_df = dtm_df.round(3)

        # Import satellite data (Download from https://apps.sentinel-hub.com/eo-browser/)
        ndvi_b04_fp = os.path.join(
            self.TerrainFilePath,
            r"EO_Browser_images\2023-03-06-00_00_2023-03-06-23_59_Sentinel-2_L2A_B04_(Raw).tiff",
        )
        ndvi_b08_fp = os.path.join(
            self.TerrainFilePath,
            r"EO_Browser_images\2023-03-06-00_00_2023-03-06-23_59_Sentinel-2_L2A_B08_(Raw).tiff",
        )
        with rio.open(ndvi_b04_fp) as src:
            ndvi_b04 = src.read(1)
            bounds = src.bounds
            meta = src.meta.copy()
        with rio.open(ndvi_b08_fp) as src:
            ndvi_b08 = src.read(1)

        # Calculate NDVI
        ndvi_values = (ndvi_b08 - ndvi_b04) / (ndvi_b08 + ndvi_b04)


        def calculate_ndvi_mean(self,meta,bounds,ndvi_values)->pd.DataFrame:

            # Convert projection from epsg:4326 to epsg:3826
            transform, _, _ = calculate_default_transform(
                meta["crs"], {"init": "EPSG:3826"}, meta["width"], meta["height"],
                *bounds)

            # zonal statistics for getting mean NDVI of each grid
            ndvi_zonal = rasterstats.zonal_stats(
                self.gridDF[["geometry"]],  # geometry object/ epsg:3826
                ndvi_values,  # ndarray
                affine=transform,  # the transformed affine (epsg:3826)
                stats="mean",  # Get the mean NDVI by grid
                nodata=-999,
            )

            # Assign gridid for the output of zonal statistics
            ndvi_zonal_df = pd.concat(
                [
                    self.gridDF[["gridid"]],
                    pd.DataFrame(ndvi_zonal, columns=["mean"])
                ],
                axis="columns",
            )
            return ndvi_zonal_df

        def calculate_ndvi_coverage(self, ndvi_values)->pd.DataFrame:

            # Defined tree and non-tree pixels code
            tree_value = 2
            non_tree_value = 1

            # Classify tree and non-tree by NDVI_THRESHOLD
            classified_data = np.where(ndvi_values >= constants.NDVI_THRESHOLD,
                                    tree_value, non_tree_value)

            # Polygonize the raster data with tree value
            shapes = list(
                rio.features.shapes(classified_data, transform=meta["transform"]))
            polygons = [
                shape(geom) for geom, value in shapes if value == tree_value
            ]

            # Create a new geodataframe for vectorized polygon
            tree_vector = (gpd.GeoDataFrame(
                {
                    "fid": [tree_value]
                }, geometry=[MultiPolygon(polygons)
                            ]).set_crs("epsg:4326").to_crs("epsg:3826"))

            # Intersect between grid and polygons
            grid_sindex = self.gridDF.sindex
            ndvi_coverage_df = pd.DataFrame(columns=["gridid", "area"])
            ndvi_df_list = []

            tree_idx = grid_sindex.query(tree_vector.geometry,
                                        predicate="intersects")

            for inputIdx, treeIdx in zip(
                    tree_idx[0],
                    tree_idx[1]):  # idx[0]:input geometry, idx[1]:tree geometry
                selected_input = tree_vector.iloc[inputIdx]
                selected_tree = self.gridDF.iloc[treeIdx]
                intersects = selected_input.geometry.intersection(
                    selected_tree.geometry)
                gridid = selected_tree["gridid"]
                area = intersects.area
                row = {"gridid": gridid, "area": area}
                ndvi_df_list.append(row)

            ndvi_coverage_df = pd.concat([ndvi_coverage_df, pd.DataFrame(ndvi_df_list)],
                            ignore_index=True)

            ndvi_coverage_df = (self.gridDF[["gridid", "Area"]].merge(ndvi_coverage_df,
                                                            how="left",
                                                            on="gridid").fillna(0))
            # Adding new columns
            ndvi_coverage_df["coverage"] = ndvi_coverage_df["area"] / ndvi_coverage_df["Area"]
            ndvi_coverage_df.drop(["area", "Area"], axis=1, inplace=True)

            return ndvi_coverage_df

        ndvi_zonal_df = self.calculate_ndvi_mean(meta, bounds, ndvi_values)
        ndvi_coverage_df = self.calculate_ndvi_coverage(ndvi_values)

        # Concat the coverage and zonal together
        ndvi_df = pd.merge(ndvi_zonal_df, ndvi_coverage_df, how="inner", on="gridid")
        ndvi_df = ndvi_df.round(3)

        print(ndvi_df)
        terrain_df = dtm_df.merge(ndvi_df, how="inner", on="gridid")

        # with fiona.open(
        #     "output.shp", "w", "ESRI Shapefile", shp_schema, "epsg:4326"
        # ) as shp:
        #     pixel_value = 2
        #     polygons = [shape(geom) for geom, value in shapes if value == pixel_value]
        #     multipolygon = MultiPolygon(polygons)
        #     shp.write(
        #         {
        #             "geometry": mapping(multipolygon),
        #             "properties": {"pixelvalue": int(pixel_value)},
        #         }
        #     )

        # with rio.open(NDVI_output_fp, "w", **meta) as dst:
        #     dst.write(classified_data, 1)

        # src_ds = gdal.Open(NDVI_output_fp)
        # srcband = src_ds.GetRasterBand(1)
        # dst_layername = "polygonized"
        # drv = ogr.GetDriverByName("ESRI Shapefile")
        # dst_ds = drv.CreateDataSource(dst_layername + ".shp")
        # dst_layer = dst_ds.CreateLayer(dst_layername, srs=None)
        # gdal.Polygonize(srcband, None, dst_layer, -1, [], callback=None)

        # NDVI = QgsRasterLayer(NDVI_output_fp, "NDVI", "gdal")

        # Vectorize NDVI
        # params = {
        #     "grid": grid_layer,
        #     "rasterlayer": NDVI,
        #     "Output": "TEMPORARY_OUTPUT",
        # }
        # NDVI_result = processing.run("customProvider:vectorize_NDVI", params)
        # target_fields = ["gridid", "coverage"]
        # NDVI_df = QgsVectorLayer_to_dataframe(NDVI_result["Output"], target_fields)
        # print(NDVI_df)

        # with rio.open(NDVI_result["Output"]) as src:
        #     print(src.read(1))
        # target_field_name = ["NDVI_mean", "coverage"]
        # data_container = []
        # for feature in NDVI_result["Output"].getFeatures():
        #     values = [feature[field_name] for field_name in target_field_name]
        #     data_container.append(values)
        # df_NDVI = pd.DataFrame(data_container, columns=target_field_name)
        # print(df_NDVI)

        # NDVI = pd.read_csv(f"{self.TerrainFilePath}/FInal_Output3.csv")
        # terrain_df = DTM.merge(NDVI, how="left", on="gridid")
        # terrain_df.rename(
        #     columns={"S": "Slope", "CLASS": "SlopeClass", " Coverage": "Coverage"},
        #     inplace=True,
        # )

        print("Done!")

        return terrain_df

    # development(容積)
    def Development(self, land_df=None):
        print("Preprocessing for Development Data...")

        # read data
        building = gpd.read_file(os.path.join(self.DevelopmentFilePath,
                                              "building.gpkg"),
                                 layer="building")
        # Fix building geometry
        building.geometry = building.geometry.buffer(0)

        if not land_df:
            land_df = self.__LandCalculate()
        land_df = gpd.GeoDataFrame(land_df,
                                   geometry="geometry").set_crs("epsg:3826")

        print("land_df Get")

        input_columns = ["floor"]
        tree_columns = ["gridid", "code"]
        development_df = self.sindex_intersection(building,
                                                  land_df,
                                                  input_columns,
                                                  tree_columns,
                                                  area=True)

        # development_df = pd.DataFrame(
        #     columns=["gridid", "code", "area", "floor", "geometry"]
        # )
        # development_df_list = []

        # land_df_sindex = land_df.sindex

        # idx = land_df_sindex.query(building.geometry, predicate="intersects")

        # for index, (inputIdx, treeIdx) in enumerate(
        #     zip(idx[0], idx[1])
        # ):  # idx[0]:input geometry, idx[1]:tree geometry
        #     selected_input = building.iloc[inputIdx]
        #     selected_tree = land_df.iloc[treeIdx]
        #     intersects = selected_input.geometry.intersection(selected_tree.geometry)

        #     # eliminate invalid types
        #     if intersects.geom_type in ["Polygon", "MultiPolygon"]:
        #         floor = selected_input["1_floor"]
        #         if floor:
        #             gridid = selected_tree["gridid"]
        #             code = selected_tree["code"]
        #             area = intersects.area
        #             geometry = intersects
        #             row = {
        #                 "gridid": gridid,
        #                 "code": code,
        #                 "area": area,
        #                 "floor": floor,
        #                 "geometry": geometry,
        #             }
        #             development_df_list.append(row)

        # development_df = pd.concat(
        #     [development_df, pd.DataFrame(development_df_list)], ignore_index=True
        # )
        development_df[
            "floor_area"] = development_df["floor"] * development_df["area"]

        print(development_df)

        development_df_pt = development_df.pivot_table(index="gridid",
                                                       columns="code",
                                                       values="floor_area",
                                                       aggfunc="sum").fillna(0)

        # Sum of total building area
        development_df_pt["sum"] = development_df_pt.sum(axis=1)

        # Filter selected columns
        selected_category = [key for key in constants.BUILDING_CATEGORY.keys()]
        building_category = selected_category.copy()
        building_category.append("sum")

        # Filtered table
        development_df_pt = development_df_pt[building_category]
        development_df_pt["others"] = development_df_pt[
            "sum"] - development_df_pt.loc[:, selected_category].sum(
                axis=1)  # Others = Total - selected categories

        # Rename columns
        development_df_pt.rename(columns=constants.BUILDING_CATEGORY,
                                 inplace=True)

        # Merge back to the origin grid
        development_df_pt = (self.gridDF[["gridid",
                                          "Area"]].merge(development_df_pt,
                                                         how="left",
                                                         on="gridid").fillna(0))

        # Add ratio columns
        for value in constants.BUILDING_CATEGORY.values():
            development_df_pt[f"{value}_ratio"] = (development_df_pt[value] /
                                                   development_df_pt["Area"])
        development_df_pt["others_ratio"] = (development_df_pt["others"] /
                                             development_df_pt["Area"])
        development_df_pt["all_ratio"] = (development_df_pt["all"] /
                                          development_df_pt["Area"])

        development_df_pt = development_df_pt.round(3)
        development_df_pt = development_df_pt.drop(["Area"], axis=1)

        print(development_df_pt)

        print("Done!")
        return development_df

    def sindex_intersection(
        self,
        input_data: gpd.GeoDataFrame,
        tree_data: gpd.GeoDataFrame,
        input_columns: list = [],
        tree_columns: list = [],
        area: bool = True,
        geometry: bool = True,
    ) -> pd.DataFrame:
        """
        Using spatail index to intersect between two polygon data,
        tree_data means the data created spatial index
        input_data means the input of the function sindex.query()

        Args:
            input_data:
            tree_data:
            input_columns: the column you wanted to keep from the input data
            tree_columns: the column you wanted to keep from the tree data
            area: If true, return the area of each intersected polygon as a new column in the dataframe.
            geometry: If true, return the geometry info as a new column in the dataframe.
        Returns:

        Raises:

        """
        # the columns for the optput dataframe
        columns = input_columns + tree_columns
        if area:
            columns.append("area")
        if geometry:
            columns.append("geometry")
        outputDF = pd.DataFrame(columns=columns)

        # create the spatial index for the tree data
        sindex = tree_data.sindex
        # Using spatial index to intersects
        idx = sindex.query(input_data["geometry"], predicate="intersects")

        row_list = []
        for inputIdx, treeIdx in zip(idx[0], idx[1]):
            selected_input = input_data.iloc[inputIdx]
            selected_tree = tree_data.iloc[treeIdx]
            intersects = selected_tree.geometry.intersection(
                selected_input.geometry)

            # Eliminate the invalid geometry types and insert row data information
            row = {}
            if intersects.geom_type in ["Polygon", "MultiPolygon"]:
                if input_columns:
                    for col in input_columns:
                        row[col] = selected_input[col]
                if tree_columns:
                    for col in tree_columns:
                        row[col] = selected_tree[col]
                if area:
                    row["area"] = intersects.area
                if geometry:
                    row["geometry"] = intersects

                row_list.append(row)

        outputDF = pd.concat([outputDF, pd.DataFrame(row_list)],
                             axis=0,
                             ignore_index=True)
        return outputDF

    # 計算土地使用by grid
    def __LandCalculate(self):
        """Calculate land_df"""

        # read data
        landuse = gpd.read_file(f"{self.LandFilePath}/landuse_108.gpkg",
                                layer="landuse_108")

        # mapping landuse code to landuse category code
        landuse["mapping_category"] = landuse["code"].map(
            lambda code: constants.LANDUSE_CATEGORY.get(str(code), "unknown"))

        input_columns = ["code", "category"]
        tree_columns = ["gridid"]

        land_df = self.sindex_intersection(landuse,
                                           self.gridDF,
                                           input_columns,
                                           tree_columns,
                                           area=True,
                                           geometry=False)

        return land_df

    # Land(土地使用現況)
    def Land(self):
        print("Preprocessing for Land Data...")

        land_df = self.__LandCalculate()

        land_output_df = land_df.pivot_table(index="gridid",
                                           columns="category",
                                           values="area",
                                           aggfunc="sum").fillna(0)

        land_output_df = self.gridDF[["gridid", "Area"]].merge(land_output_df,
                                                             how="left",
                                                             on="gridid")

        for code, name in constants.LANDUSE_CATEGORY_CODE.items():
            land_output_df[f"{name}Ratio"] = (land_output_df[str(code)] /
                                            land_output_df["Area"])
            land_output_df.drop(str(code), axis=1, inplace=True)

        land_output_df = land_output_df.round(3)

        print("Done!")
        return land_output_df

    # POI
    def POI(self):
        print("Preprocessing for POI Data...")

        # read data
        poi_filepath_list = [
            file for file in glob.glob(f"{self.PoiFilePath}/*.csv")
        ]

        def poi_preprocessing(poi_filepath):
            # Get the POI name by the original file name
            layer_name = "".join(
                re.search(r"poi_(.*?)\.csv",
                          poi_filepath).group(1).capitalize().split("_"))
            poi_df = pd.read_csv(poi_filepath)

            # Convert csv to GeoDataFrame
            poi_gdf = (
                gpd.GeoDataFrame(
                    poi_df[["rating_num", "rating"]],
                    geometry=gpd.points_from_xy(poi_df.lng.astype(float),
                                                poi_df.lat.astype(float)),
                ).set_crs("epsg:4326")  # original crs
                .to_crs("epsg:3826")  # destinate crs
            )

            # POIs intersect with grid and compute the aggregation
            poi_agg_df = (
                gpd.sjoin(
                    self.gridDF[["gridid", "geometry"]],
                    poi_gdf,
                    predicate="intersects",
                    how="inner",
                ).loc[
                    lambda x: x["rating_num"] >= constants.
                    POI_FILTER_NUM]  # filter the POI rating less than the threshold POI_FILTER_NUM
                .groupby(["gridid"]).agg({
                    "rating_num": ["count", "sum"],  # create multiIndex columns
                    "rating": "sum",
                }).droplevel(
                    0, axis=1
                )  # flatten the multiIndex columns built by agg() function
                .reset_index(
                )  # let "gridid" as a new column instead of the index
            )

            # Update the columns name
            poi_agg_df.columns = [
                "gridid",
                f"{layer_name}POICounts",
                f"{layer_name}POIRatingCountsSum",
                f"{layer_name}POIRatingStarSum",
            ]

            # Merge with original grid, if no data, fill 0.
            poi_merge_df = pd.merge(
                self.gridDF[["gridid"]],
                poi_agg_df,
                on="gridid",
                how="left",  # ramain all gridid
            ).fillna(0)

            return poi_merge_df

        poi_df = pd.DataFrame()
        for filepath in poi_filepath_list:
            poi_merge_df = poi_preprocessing(filepath)
            if poi_df.empty:
                poi_df = poi_merge_df
            else:
                poi_df = pd.merge(poi_df, poi_merge_df, on="gridid")

        # # Setting
        # PoiDF = pd.DataFrame()
        # column_definition = {
        #     "rating_num": "Int",
        #     "rating": "Double",
        # }

        # def _poi_preprocessing(poi_filepath):
        #     layer_name = "".join(
        #         re.search(r"poi_(.*?)\.csv", poi_filepath)
        #         .group(1)
        #         .capitalize()
        #         .split("_")
        #     )
        #     vlyr = csv_to_QgsVectorLayer(
        #         poi_filepath,  # file path
        #         layer_name,
        #         "lat",  # The latitude column name in the original dataset
        #         "lng",  # The longitude column name in the original dataset
        #         **column_definition,  # The column wanted to be extracted {name: QVarint type}
        #     )
        #     params = {
        #         "filternumber": 5,
        #         "grid": self.gridFilePath,
        #         "pois": [vlyr],
        #         "Output": "TEMPORARY_OUTPUT",
        #     }
        #     result = processing.run("customProvider:aggregate_poi_to_grid", params)

        #     target_fields = [
        #         "gridid",
        #         "rating_num_count",
        #         "rating_num_sum",
        #         "rating_star_sum",
        #     ]
        #     df = QgsVectorLayer_to_dataframe(result["Output"], target_fields)

        #     df = df.rename(
        #         columns={
        #             "rating_num_count": f"{layer_name}POICounts",
        #             "rating_num_sum": f"{layer_name}POIRatingCountsSum",
        #             "rating_star_sum": f"{layer_name}POIRatingStarSum",
        #         }
        #     )
        #     return df

        # for filepath in poi_filepath_list:
        #     temDF = __poi_preprocessing(filepath)
        #     try:
        #         PoiDF = PoiDF.merge(temDF, how="inner", on="gridid")
        #     except:
        #         PoiDF = temDF

        print("Done!")
        print(poi_df)

        return poi_df

    # save preprocessed data
    def saveDF(self, DF):
        DF.to_csv(f"{self.DoneFilePath}/DF.csv", index=False)
        print("Save the dataframe after preprocessing!")

    # pipeline of preprocess
    def run(self):
        # create Base DF with complate gridID, Date, Hour
        GRIDID = list(set(self.gridDF["gridid"]))
        HOUR = [x for x in range(24)]
        productElements = list(itertools.product(GRIDID, self.DateList, HOUR))
        KEYS = ["GridID", "Date", "Hour"]
        DF = pd.DataFrame(productElements, columns=KEYS)
        DF["Weekday"] = pd.to_datetime(DF["Date"], format="%Y%m%d").dt.weekday
        DF["IsWeekend"] = DF["Weekday"].apply(lambda x: 1 if x > 4 else 0)

        # ## transaction
        # TranOnDF, TranOffDF = self.Transaction()
        # TranOnKeys, TranOffKeys = ["gridid", "on_date", "on_hour"], [
        #     "gridid",
        #     "off_date",
        #     "off_hour",
        # ]
        # TranOnDF = (
        #     TranOnDF[TranOnKeys + ["counts"]].groupby(TranOnKeys).sum().reset_index()
        # )
        # TranOffDF = (
        #     TranOffDF[TranOffKeys + ["counts"]].groupby(TranOffKeys).sum().reset_index()
        # )

        # # OnCounts
        # DF = DF.merge(TranOnDF, how="left", left_on=KEYS, right_on=TranOnKeys)
        # DF = DF[["GridID", "Date", "Hour", "IsWeekend", "counts"]]
        # DF.rename(columns={"counts": "OnCounts"}, inplace=True)
        # temDFCols = list(DF.columns)

        # # OffCounts
        # DF = DF.merge(TranOffDF, how="left", left_on=KEYS, right_on=TranOffKeys)
        # DF = DF[temDFCols + ["counts"]]
        # DF.rename(columns={"counts": "OffCounts"}, inplace=True)

        # # NetCounts
        # DF.fillna(0, inplace=True)
        # DF["NetCounts"] = DF["OffCounts"] - DF["OnCounts"]

        # ## Population
        # # Age_15_17_Counts, Age_18_21_Counts, ...,
        # # Age_Over65_Counts, Age_Total_Counts,
        # # WorkPopulationCounts, LivePopulationCounts, TourPopulationCounts
        # PopulationDF = self.Population()
        # PopuOffKeys = ["網格編號", "日期", "時間"]
        # DF = DF.merge(PopulationDF, how="left", left_on=KEYS, right_on=PopuOffKeys)
        # DF.drop(PopuOffKeys, axis=1, inplace=True)
        # DF.fillna(0, inplace=True)

        # ## Traffic MRT
        # OverlayExit_DF, populationIN, populationOUT = self.TrafficMRT()
        # MRTExitCountsDF = (
        #     OverlayExit_DF.drop("geometry", axis=1)
        #     .groupby("gridid")
        #     .size()
        #     .reset_index(name="counts")
        # )

        # populationIN["日期"] = populationIN["日期"].str.split("-").str.join("")
        # populationOUT["日期"] = populationOUT["日期"].str.split("-").str.join("")

        # MRTPopulationCols = ["gridid", "日期", "時段"]
        # MRTPopulationInCountsDF = (
        #     populationIN[MRTPopulationCols + ["人次"]]
        #     .groupby(MRTPopulationCols)
        #     .sum()
        #     .reset_index()
        # )
        # MRTPopulationOutCountsDF = (
        #     populationOUT[MRTPopulationCols + ["人次"]]
        #     .groupby(MRTPopulationCols)
        #     .sum()
        #     .reset_index()
        # )

        # # MRTExitCounts
        # DF = DF.merge(MRTExitCountsDF, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)
        # DF.rename(columns={"counts": "MRTExitCounts"}, inplace=True)
        # DF.fillna(0, inplace=True)
        # DF["MRTExitCounts"] = DF["MRTExitCounts"].astype(int)

        # # MRTPopulationInCounts
        # DF = DF.merge(
        #     MRTPopulationInCountsDF,
        #     how="left",
        #     left_on=KEYS,
        #     right_on=MRTPopulationCols,
        # )
        # DF.drop(MRTPopulationCols, axis=1, inplace=True)
        # DF.rename(columns={"人次": "MRTPopulationInCounts"}, inplace=True)
        # DF.fillna(0, inplace=True)
        # DF["MRTPopulationInCounts"] = DF["MRTPopulationInCounts"].astype(int)

        # # MRTPopulationOutCounts
        # DF = DF.merge(
        #     MRTPopulationOutCountsDF,
        #     how="left",
        #     left_on=KEYS,
        #     right_on=MRTPopulationCols,
        # )
        # DF.drop(MRTPopulationCols, axis=1, inplace=True)
        # DF.rename(columns={"人次": "MRTPopulationOutCounts"}, inplace=True)
        # DF.fillna(0, inplace=True)
        # DF["MRTPopulationOutCounts"] = DF["MRTPopulationOutCounts"].astype(int)

        # ## Traffic Bus
        # OverlayBusStop_DF = self.TrafficBus()
        # BusStopCountsDF = (
        #     OverlayBusStop_DF.groupby("gridid").size().reset_index(name="counts")
        # )

        # # BusStopCounts
        # DF = DF.merge(BusStopCountsDF, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)
        # DF.rename(columns={"counts": "BusStopCounts"}, inplace=True)
        # DF.fillna(0, inplace=True)
        # DF["BusStopCounts"] = DF["BusStopCounts"].astype(int)

        # ## Bus Route
        # BusRouteDF = self.TrafficBusRoute()

        # # BusRouteCounts
        # DF = DF.merge(BusRouteDF, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)

        # ## Road SideWalk
        # OverlaySideWalk_DF = self.RoadSideWalk()
        # SideWalkAreaDF = (
        #     OverlaySideWalk_DF[["gridid", "area"]].groupby("gridid").sum().reset_index()
        # )

        # # SideWalkArea
        # DF = DF.merge(SideWalkAreaDF, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)
        # DF.rename(columns={"area": "SideWalkArea"}, inplace=True)
        # DF.fillna(0, inplace=True)

        # ## Road MarkedSideWalk
        # OverlayMarkedSideWalk_DF = self.RoadMarkedSideWalk()
        # MarkedSideWalkAreaDF = (
        #     OverlayMarkedSideWalk_DF[["gridid", "area"]]
        #     .groupby("gridid")
        #     .sum()
        #     .reset_index()
        # )

        # # MarkedSideWalkArea
        # DF = DF.merge(
        #     MarkedSideWalkAreaDF, how="left", left_on="GridID", right_on="gridid"
        # )
        # DF.drop("gridid", axis=1, inplace=True)
        # DF.rename(columns={"area": "MarkedSideWalkArea"}, inplace=True)
        # DF.fillna(0, inplace=True)

        # ## Road BikeRoute
        # OverlayBikeRoute_DF = self.RoadBikeRoute()
        # BikeRouteLengthDF = (
        #     OverlayBikeRoute_DF[["gridid", "length"]]
        #     .groupby("gridid")
        #     .sum()
        #     .reset_index()
        # )

        # # BikeRouteLength
        # DF = DF.merge(
        #     BikeRouteLengthDF, how="left", left_on="GridID", right_on="gridid"
        # )
        # DF.drop("gridid", axis=1, inplace=True)
        # DF.rename(columns={"length": "BikeRouteLength"}, inplace=True)
        # DF.fillna(0, inplace=True)

        # ## Road Tree
        # OverlayTree_DF = self.RoadTree()
        # TreeCountsDF = (
        #     OverlayTree_DF.groupby("gridid").size().reset_index(name="counts")
        # )

        # # TreeCounts
        # DF = DF.merge(TreeCountsDF, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)
        # DF.rename(columns={"counts": "TreeCounts"}, inplace=True)
        # DF.fillna(0, inplace=True)
        # DF["TreeCounts"] = DF["TreeCounts"].astype(int)

        # ## Road Light
        # OverlayLight_DF = self.RoadLight()
        # LightCountsDF = (
        #     OverlayLight_DF.groupby("gridid").size().reset_index(name="counts")
        # )

        # # LightCounts
        # DF = DF.merge(LightCountsDF, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)
        # DF.rename(columns={"counts": "LightCounts"}, inplace=True)
        # DF.fillna(0, inplace=True)
        # DF["LightCounts"] = DF["LightCounts"].astype(int)

        # ## Road Net
        # RoadNetDF = self.RoadRoadNet()

        # # RoadArea, RoadAreaRatio
        # DF = DF.merge(RoadNetDF, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)

        # ## Road Length
        # RoadLengthDF = self.RoadRoadLength()

        # # ExpressWayLength, ProvincialHighwayLength, UrbanRoad_RoadStreetLength, UrbanRoad_LaneAlleyLength, TotalRoadLength
        # DF = DF.merge(RoadLengthDF, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)

        ## Terrain
        terrain_df = self.Terrain()

        # # Slope, SlopeClass, ELEV_mean, NDVImean, Coverage
        # DF = DF.merge(terrain_df, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)

        ## Development
        # development_df = self.Development()

        # # commercial_floor_area, commercial_ratio, mixed_floor_area, mixed_ratio, ...
        # DF = DF.merge(development_df, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)

        # ## Land
        # land_df = self.Land()

        # # NatureRatio, CommerceRatio, ResidenceRatio, ...
        # DF = DF.merge(land_df, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)

        ## POI
        # PoiDF = self.POI()

        # # TouristAttractionRatingNumCount, TouristAttractionRatingNumSum, TouristAttractionRatingStarSum, ...
        # DF = DF.merge(PoiDF, how="left", left_on="GridID", right_on="gridid")
        # DF.drop("gridid", axis=1, inplace=True)

        ## save
        # self.saveDF(DF)


#         return DF

if __name__ == "__main__":
    # Initialize QGIS
    qgs = QgsApplication([], False)
    qgs.initQgis()
    # Initialize QGIS processing framework
    Processing.initialize()
    # Add the custom provider and custom algorithms
    provider = CustomProvider()
    QgsApplication.processingRegistry().addProvider(provider)
    # QgsApplication.processingRegistry().addProvider(Grass7AlgorithmProvider())

    p = Preprocess(DateList=["20230305", "20230311", "20230317", "20230322"])
    #     TranOnDF, TranOffDF = p.Transaction()
    #     PopuDF = p.Population()
    #     OverlayExit_DF, populationIN, populationOUT = p.TrafficMRT()
    #     OverlayBusStop_DF = p.TrafficBus()
    #     # p.RoadRoadNet()
    #     OverlaySideWalk_DF = p.RoadSideWalk()
    #     OverlayMarkedSideWalk_DF = p.RoadMarkedSideWalk()
    #     OverlayBikeRoute_DF = p.RoadBikeRoute()
    #     OverlayTree_DF = p.RoadTree()
    #     OverlayLight_DF = p.RoadLight()

    StartTime = str(datetime.datetime.now())
    print("Data Preprocessing Start at {StartTime}!\n".format(
        StartTime=StartTime))

    p.run()

    EndTime = str(datetime.datetime.now())
    print("\nData Preprocessing Finished at {EndTime}!".format(EndTime=EndTime))

    qgs.exitQgis()
