import pandas as pd
import os
import pickle
import itertools
import twd97
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from tqdm.notebook import tqdm
import warnings
import xmltodict
import pykml.parser
import json
import datetime

pd.set_option('display.max_columns', None)
warnings.filterwarnings("ignore")


class Preprocess:

    def __init__(self, DateList):
        # date
        self.DateList = DateList
        self.HyphenDateList = [
            x[0:4] + '-' + x[4:6] + '-' + x[6:8] for x in self.DateList
        ]

        # file path
        self.TransactionFilePath = 'input/data/transaction'
        self.PopulationFilePath = 'input/data/population'
        self.TrafficFilePath = 'input/data/traffic'
        self.RoadFilePath = 'input/data/road'
        self.TerrainFilePath = 'input/data/terrain'
        self.DevelopmentFilePath = 'input/data/development'
        self.LandFilePath = 'input/data/land'
        self.PoiFilePath = 'input/data/POI'
        self.DoneFilePath = 'output_OFF'

        # grid information
        self.gridDF = gpd.read_file(
            f"{self.PopulationFilePath}/FET_2023_grid_97.geojson")

    def Point_WGS84toTWD97(self, DF, colName, reverse):
        # get all points' wgs84 coordinate
        # e.g. "POINT (121.56654 25.08235)", "POINT (121.50449 25.10045)", ...

        wgs84List = DF[colName].astype(str).str[6:].str.replace(' ',
                                                                ', ').tolist()
        wgs84List = [eval(x) for x in wgs84List]

        # WGS84: lng, lat
        if reverse:
            twd97List = [twd97.fromwgs84(x[1], x[0]) for x in wgs84List]

        # WGS84: lat, lng
        else:
            twd97List = [twd97.fromwgs84(x[0], x[1]) for x in wgs84List]

        # output
        DF['twd97_lat'] = [x[0] for x in twd97List]
        DF['twd97_lng'] = [x[1] for x in twd97List]

        return DF

    # get the intersection area within two Polygon series
    def OverlayWithinPPL(self, df1, df2, key1, key2, method):

        # df1: "polygon df" or "point df" or "line df"
        # df2: "grid df"
        # method: "Polygon" or "Point" or "Line"

        overlayDF = gpd.overlay(df1, df2, how='union').explode().reset_index()
        overlayDF = overlayDF[[key1, key2, 'geometry']]
        overlayDF = overlayDF[overlayDF[key1].notnull()]
        overlayDF = overlayDF[overlayDF[key2].notnull()]
        overlayDF[key2] = overlayDF[key2].astype(int)
        overlayDF.index = range(len((overlayDF)))

        if method == 'Polygon':
            overlayDF['area'] = overlayDF.geometry.area
        elif method == 'Line':
            overlayDF['length'] = overlayDF.geometry.length
        elif method == 'Point':
            pass

        return overlayDF

    # transaction
    def Transaction(self):

        print("Preprocessing for Transaction Data...")

        # read stop data
        StopDF = pd.read_csv(
            f"{self.TransactionFilePath}/86ec099baa2d36c22ab3a87350b718de_export.csv"
        )
        StopDF = StopDF[['sno', 'lat', 'lng']]
        StopDF['sno'] = StopDF['sno'].astype(str)
        StopDF['sno'] = 'U' + StopDF['sno'].str[3:]
        StopDF['geometry'] = [
            Point(twd97.fromwgs84(x, y))
            for x, y in zip(StopDF.lat, StopDF.lng)
        ]

        # read transaction data
        DF = pd.DataFrame()
        for date in self.DateList:
            tem = open(
                f"{self.TransactionFilePath}/202303_txn_identified_transfer/{date}.pkl",
                'rb')
            temFile = pickle.load(tem)
            DF = pd.concat([DF, temFile])
            tem.close()
        DF.index = range(0, len(DF))

        # convert on_time(datetime), off_time(datetime) to on_hour(int), off_hour(int)
        # hardcore the location of date and hour
        DF['on_date'] = DF['on_time'].astype(str).str[:10].str.split(
            '-').str.join('')
        DF['off_date'] = DF['off_time'].astype(str).str[:10].str.split(
            '-').str.join('')
        DF['on_hour'] = DF['on_time'].astype(str).str[11:13].astype(int)
        DF['off_hour'] = DF['off_time'].astype(str).str[11:13].astype(int)

        # complete mapping table
        ON_STOP_ID = list(set(DF['on_stop_id']))
        ONDATE, OFFDATE = list(set(DF['on_date'])), list(set(DF['off_date']))
        ONHOUR, OFFHOUR = list(set(DF['on_hour'])), list(set(DF['off_hour']))
        onColumns, offColumns = ['on_stop_id', 'on_date', 'on_hour'
                                ], ['off_stop_id', 'off_date', 'off_hour']
        productElements = list(itertools.product(ON_STOP_ID, ONDATE, ONHOUR))
        OnMappingDF = pd.DataFrame(productElements, columns=onColumns)
        OffMappingDF = pd.DataFrame(productElements, columns=offColumns)

        # on (借車)
        onDF = DF.groupby(onColumns).size().reset_index(name='counts')
        onDF = onDF.merge(OnMappingDF, how='right', on=onColumns)
        onDF.fillna(0, inplace=True)

        # off (還車)
        offDF = DF.groupby(offColumns).size().reset_index(name='counts')
        offDF = offDF.merge(OffMappingDF, how='right', on=offColumns)
        offDF.fillna(0, inplace=True)

        # stop_id's information include lat, lng
        onDF = onDF.merge(StopDF,
                          how='left',
                          left_on='on_stop_id',
                          right_on='sno').drop(['sno'], axis=1)
        offDF = offDF.merge(StopDF,
                            how='left',
                            left_on='off_stop_id',
                            right_on='sno').drop(['sno'], axis=1)

        # get "on" geometry, "off" geometry and Grid geometry
        On_point = gpd.GeoDataFrame(onDF[['on_stop_id', 'geometry']])
        Off_point = gpd.GeoDataFrame(offDF[['off_stop_id', 'geometry']])
        Grid_poly = self.gridDF[['gridid', 'geometry']]

        # drop duplicate
        On_point.drop_duplicates(inplace=True)
        Off_point.drop_duplicates(inplace=True)
        On_point.index = range(len(On_point))
        Off_point.index = range(len(Off_point))

        # overlay
        OverlayOn_DF = self.OverlayWithinPPL(On_point,
                                             Grid_poly,
                                             'on_stop_id',
                                             'gridid',
                                             method='Point')
        OverlayOff_DF = self.OverlayWithinPPL(Off_point,
                                              Grid_poly,
                                              'off_stop_id',
                                              'gridid',
                                              method='Point')

        # mapping
        onDF = onDF.merge(OverlayOn_DF.drop(['geometry'], axis=1),
                          how='left',
                          on='on_stop_id')
        offDF = offDF.merge(OverlayOff_DF.drop(['geometry'], axis=1),
                            how='left',
                            on='off_stop_id')

        # kick out NA, reset index and astype
        onDF = onDF[onDF.gridid.notnull()]
        offDF = offDF[offDF.gridid.notnull()]
        onDF.index = range(len(onDF))
        offDF.index = range(len(offDF))
        onDF[['counts', 'gridid']] = onDF[['counts', 'gridid']].astype(int)
        offDF[['counts', 'gridid']] = offDF[['counts', 'gridid']].astype(int)

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
        populationDF['日期'] = populationDF['日期'].str.split('-').str.join('')
        workDF['日期'] = workDF['日期'].str.split('-').str.join('')
        liveDF['日期'] = liveDF['日期'].str.split('-').str.join('')
        tourDF['日期'] = tourDF['日期'].str.split('-').str.join('')

        # by age and total
        DF = populationDF.groupby(['日期', '時間', '網格編號',
                                   '年齡別']).sum('放大後人數').reset_index()
        DF = pd.pivot_table(DF,
                            values='放大後人數',
                            index=['日期', '時間', '網格編號'],
                            columns='年齡別').reset_index().fillna(0)
        DF.columns = [
            '日期', '時間', '網格編號', 'Age_15_17_Counts', 'Age_18_21_Counts',
            'Age_22_29_Counts', 'Age_30_39_Counts', 'Age_40_49_Counts',
            'Age_50_59_Counts', 'Age_60_64_Counts', 'Age_Over65_Counts'
        ]
        DF['Age_Total_Counts'] = DF['Age_15_17_Counts'] + DF[
            'Age_18_21_Counts'] + DF['Age_22_29_Counts'] + DF[
                'Age_30_39_Counts'] + DF['Age_40_49_Counts'] + DF[
                    'Age_50_59_Counts'] + DF['Age_60_64_Counts'] + DF[
                        'Age_Over65_Counts']

        mergeCol = ['日期', '時間', '網格編號']

        # by work
        DF = DF.merge(workDF, how='left', on=mergeCol)
        DF.rename(columns={'放大後人數': 'WorkPopulationCounts'}, inplace=True)

        # by live
        DF = DF.merge(liveDF, how='left', on=mergeCol)
        DF.rename(columns={'放大後人數': 'LivePopulationCounts'}, inplace=True)

        # by tour
        DF = DF.merge(tourDF, how='left', on=mergeCol)
        DF.rename(columns={'放大後人數': 'TourPopulationCounts'}, inplace=True)

        print("Done!")

        return DF

    # traffic for MRT
    def TrafficMRT(self):

        print("Preprocessing for MRT Data...")

        # read data
        MRT_DF = pd.read_csv(f"{self.TrafficFilePath}/臺北捷運車站出入口座標.csv",
                             encoding='big5')
        population = pd.read_csv(
            f"{self.TrafficFilePath}/臺北捷運每日分時各站OD流量統計資料_202303.csv")

        ## information of MRT Station and Exit
        # generate twd97 coordinate
        MRT_DF['Exit_twd97'] = MRT_DF.apply(
            lambda x: twd97.fromwgs84(x.緯度, x.經度), axis=1)
        MRT_DF[['Exit_twd97_lat',
                'Exit_twd97_lng']] = MRT_DF['Exit_twd97'].tolist()

        # get Station Name
        MRT_DF['Station'] = MRT_DF['出入口名稱'].str.split('站出口').str[0]

        # split MRT_DF into Exit_DF and Station_DF
        Exit_DF = MRT_DF.drop(columns=['項次', 'Station'])
        Station_DF = MRT_DF[['Exit_twd97_lat', 'Exit_twd97_lng', 'Station']].\
                             groupby('Station').\
                             agg(Station_twd97_lat=('Exit_twd97_lat', 'mean'),
                                 Station_twd97_lng=('Exit_twd97_lng', 'mean')).reset_index()

        # create column 'geometry'
        Exit_DF['geometry'] = Exit_DF.Exit_twd97.apply(Point)
        Station_DF['geometry'] = [Point(x) for x in \
                                  list(zip(Station_DF.Station_twd97_lat, Station_DF.Station_twd97_lng))]

        # get exit geometry, station geometry and Grid geometry
        Exit_point = gpd.GeoDataFrame(Exit_DF[['出入口名稱', 'geometry']])
        Station_point = gpd.GeoDataFrame(Station_DF[['Station', 'geometry']])
        Grid_poly = self.gridDF[['gridid', 'geometry']]

        # overlay
        OverlayExit_DF = self.OverlayWithinPPL(Exit_point,
                                               Grid_poly,
                                               '出入口名稱',
                                               'gridid',
                                               method='Point')
        OverlayStation_DF = self.OverlayWithinPPL(Station_point,
                                                  Grid_poly,
                                                  'Station',
                                                  'gridid',
                                                  method='Point')

        ## population of MRT Stop (各站點進出人次)
        population = population[population['日期'].isin(self.HyphenDateList)]
        populationIN = population[['日期', '時段', '進站', '人次']]
        populationOUT = population[['日期', '時段', '出站', '人次']]
        populationIN.index = range(0, len(populationIN))
        populationOUT.index = range(0, len(populationOUT))

        # calculate number of people by day, hour, station
        populationIN = populationIN.groupby(['日期', '時段',
                                             '進站']).sum('人次').reset_index()
        populationOUT = populationOUT.groupby(['日期', '時段',
                                               '出站']).sum('人次').reset_index()

        # get GRID ID
        populationIN = populationIN.merge(OverlayStation_DF,
                                          how='inner',
                                          left_on='進站',
                                          right_on='Station')
        populationOUT = populationOUT.merge(OverlayStation_DF,
                                            how='inner',
                                            left_on='出站',
                                            right_on='Station')

        print("Done!")

        return OverlayExit_DF, populationIN, populationOUT

    # traffic for BUS
    def TrafficBus(self):

        print("Preprocessing for Bus Data...")

        # read data
        DF = gpd.read_file(f"{self.TrafficFilePath}/busstop/busstop.shp")

        # get TWD97 coordinate
        DF = self.Point_WGS84toTWD97(DF, 'geometry', reverse=True)
        DF.drop(columns=['geometry'], inplace=True)

        # create twd97 geometry series
        DF['LatLng'] = list(zip(DF.twd97_lat, DF.twd97_lng))
        DF['geometry'] = DF.LatLng.apply(Point)

        # get bus geometry and Grid geometry
        Bus_point = DF[['BSM_BUSSTO', 'geometry']]
        Grid_poly = self.gridDF[['gridid', 'geometry']]

        # overlay
        overlayDF = self.OverlayWithinPPL(Bus_point, Grid_poly, 'BSM_BUSSTO',
                                          'gridid', 'Point')

        print("Done!")

        return overlayDF

    # road for bus route
    # by 軒柏
    def TrafficBusRoute(self):
        print("Preprocessing for BusRoute Data...")

        # read data
        BusRouteDF = pd.read_csv(f"{self.TrafficFilePath}/bus_final.csv")
        BusRouteDF.drop(columns=["RouteNum_count"], axis=1, inplace=True)
        BusRouteDF.rename(columns={"RouteNum_sum": "BusRouteCounts"},
                          inplace=True)

        print("Done!")

        return BusRouteDF

    # road for side walk
    def RoadSideWalk(self):

        print("Preprocessing for SideWalk Data...")

        # read SideWalk data
        f = open(f"{self.RoadFilePath}/TP_SIDEWORK.json")
        sidewalkDF = json.load(f)
        sidewalkDF = pd.json_normalize(sidewalkDF['features'])
        sidewalkDF = gpd.GeoDataFrame(sidewalkDF)
        sidewalkDF['geometry'] = [
            Polygon(x[0][0]) for x in sidewalkDF['geometry.coordinates']
        ]

        # get sidewalk geometry and Grid geometry
        SW_poly = sidewalkDF[['properties.ObjectID', 'geometry']]
        Grid_poly = self.gridDF[['gridid', 'geometry']]

        # overlay
        overlayDF = self.OverlayWithinPPL(SW_poly,
                                          Grid_poly,
                                          'properties.ObjectID',
                                          'gridid',
                                          method='Polygon')

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
        MSW_poly = MarkedSideWalkDF[['KEYID', 'geometry']]
        Grid_poly = self.gridDF[['gridid', 'geometry']]

        # overlay
        overlayDF = self.OverlayWithinPPL(MSW_poly,
                                          Grid_poly,
                                          'KEYID',
                                          'gridid',
                                          method='Polygon')

        print("Done!")

        return overlayDF

    # road for Bike
    def RoadBikeRoute(self):

        print("Preprocessing for Bike Route Data...")

        # read data
        root = pykml.parser.fromstring(
            open(f"{self.RoadFilePath}/台北市_自行車道-市區自行車道1120505.kml", 'r').read())
        # print(root.Document.Folder.Placemark.name)
        # print(root.Document.Folder.Placemark.LineString.coordinates)
        Name = [
            x.name.text for x in root.findall(
                './/{http://www.opengis.net/kml/2.2}Placemark')
        ]
        Coor = [
            str(x.LineString.coordinates).split(',0 ')[0:-1] for x in
            root.findall('.//{http://www.opengis.net/kml/2.2}Placemark')
        ]

        # get all twd97 coordinates of road of bike
        twd97Coor = list()
        for x in Coor:
            tem = []
            for y in x:
                e = y.split(',')
                t = twd97.fromwgs84(e[1], e[0])
                tem.append(t)
            twd97Coor.append(tem)

        DF = gpd.GeoDataFrame(data={'Name': Name, "twd97Coor": twd97Coor})

        # remove some road of bike, which's len(coordinate) is 1
        removeIndex = list()
        for x in range(len(DF)):
            if len(DF.loc[x, 'twd97Coor']) == 1:
                removeIndex.append(x)
        DF.drop(index=removeIndex, inplace=True)

        DF['geometry'] = DF.twd97Coor.apply(LineString)

        # get road of bike DataFrame and Grid DataFrame
        gdf_lines = DF[['Name', 'geometry']]
        gdf_poly = self.gridDF[['gridid', 'geometry']]

        # get intersection LineString between all road of bike and all Grid
        OverlayDF = gpd.overlay(gdf_lines, gdf_poly,
                                how='union').explode().reset_index(drop=True)

        # get Taipei DF
        TpeOverlayDF = OverlayDF[OverlayDF['gridid'].notnull()]
        TpeOverlayDF.index = range(len(TpeOverlayDF))
        TpeOverlayDF['gridid'] = TpeOverlayDF['gridid'].astype(int)

        # calaculate Length of road of bike
        TpeOverlayDF['length'] = TpeOverlayDF.geometry.length

        print("Done!")

        return TpeOverlayDF

    # road for tree
    def RoadTree(self):

        print("Preprocessing for Tree Data...")

        # read data
        f = open(f"{self.RoadFilePath}/TaipeiTree.json")
        TreeDF = json.load(f)
        TreeDF = pd.json_normalize(TreeDF)
        TreeDF['LatLng'] = list(zip(TreeDF.X, TreeDF.Y))
        TreeDF = gpd.GeoDataFrame(TreeDF)
        TreeDF['geometry'] = TreeDF.LatLng.apply(Point)

        # get Tree geometry and Grid geometry
        Tree_point = TreeDF[['TreeID', 'geometry']]
        Grid_poly = self.gridDF[['gridid', 'geometry']]

        # overlay
        overlayDF = self.OverlayWithinPPL(Tree_point,
                                          Grid_poly,
                                          'TreeID',
                                          'gridid',
                                          method='Point')

        print("Done!")

        return overlayDF

    # road for light
    def RoadLight(self):

        print("Preprocessing for Light Data...")

        # read data
        f = open(f"{self.RoadFilePath}/TaipeiLight.json")
        LightDF = json.load(f)
        LightDF = pd.json_normalize(LightDF)
        LightDF = LightDF[LightDF.X != '']
        LightDF = LightDF[LightDF.Y != '']
        LightDF.index = range(len(LightDF))
        LightDF['LatLng'] = list(zip(LightDF.X, LightDF.Y))
        LightDF = gpd.GeoDataFrame(LightDF)
        LightDF['geometry'] = LightDF.LatLng.apply(Point)

        # get Tree geometry and Grid geometry
        Light_point = LightDF[['LIGHTID', 'geometry']]
        Grid_poly = self.gridDF[['gridid', 'geometry']]

        # overlay
        overlayDF = self.OverlayWithinPPL(Light_point,
                                          Grid_poly,
                                          'LIGHTID',
                                          'gridid',
                                          method='Point')

        print("Done!")

        return overlayDF

    # road for road net
    # by 軒柏
    def RoadRoadNet(self):
        print("Preprocessing for RoadNet Data...")

        # read data
        RoadNetDF = pd.read_csv(f"{self.RoadFilePath}/grid_road_area.csv")
        RoadNetDF.rename(columns={
            "Area": "RoadArea",
            "Ratio": "RoadAreaRatio"
        },
                         inplace=True)

        print("Done!")

        return RoadNetDF

    # road for road length
    # by 軒柏
    def RoadRoadLength(self):
        print("Preprocessing for Road Length Data...")

        # read data
        RoadLengthDF = pd.read_csv(
            f"{self.RoadFilePath}/road_length_category.csv")
        RoadLengthDF.rename(
            columns={
                "class1": "ExpressWayLength",    # 一般快速道路(市區高架道路)
                "class2": "ProvincialHighwayLength",    # 省道
                "class3": "UrbanRoad_RoadStreetLength",    # 市區道路(路、街)
                "class4": "UrbanRoad_LaneAlleyLength",    # 市區道路(巷、弄)
                "Sum": "TotalRoadLength"
            }    # 總長度
            ,
            inplace=True)

        print("Done!")

        return RoadLengthDF

    # Terrain
    # by 軒柏
    def Terrain(self):

        print("Preprocessing for Terrain Data...")

        # read data
        DTM = pd.read_csv(f"{self.TerrainFilePath}/DTM_Slope_2.csv")
        NDVI = pd.read_csv(f"{self.TerrainFilePath}/FInal_Output3.csv")
        TerrainDF = DTM.merge(NDVI, how='left', on='gridid')
        TerrainDF.rename(columns={
            "S": "Slope",
            "CLASS": "SlopeClass",
            " Coverage": "Coverage"
        },
                         inplace=True)

        print("Done!")

        return TerrainDF

    # development(容積)
    # by 軒柏
    def Development(self):
        print("Preprocessing for Development Data...")

        # read data
        DevelopmentDF = pd.read_csv(
            f"{self.DevelopmentFilePath}/building_ratio_final.csv")
        DevelopmentDF.drop(columns=['Area'], axis=1, inplace=True)

        print("Done!")

        return DevelopmentDF

    # Land(土地使用現況)
    # by 軒柏
    def Land(self):
        print("Preprocessing for Land Data...")

        # read data
        LandDF = pd.read_csv(f"{self.LandFilePath}/landuse_category_area.csv")
        LandDF.rename(columns={
            "1": "NatureRatio",
            "2": "CommerceRatio",
            "3": "ResidenceRatio",
            "4": "MixedResidenceRatio",
            "5": "IndustryRatio",
            "6": "InfrastructureRatio",
            "7": "EducationRatio",
            "8": "LeisureRatio",
            "9": "OpenSpaceRatio",
            "10": "TrafficRatio",
            "11": "OtherRatio"
        },
                      inplace=True)

        print("Done!")

        return LandDF

    # POI
    # by 軒柏
    def POI(self):
        print("Preprocessing for POI Data...")

        # read data
        AllPoiDF = os.listdir(f"{self.PoiFilePath}")
        PoiDF = pd.DataFrame()

        for x in AllPoiDF:
            temName = ''.join(
                [t.capitalize() for t in x.split('.')[0].split('_')])
            temDF = pd.read_csv(f"{self.PoiFilePath}/{x}")
            if x == 'store.csv':
                temDF.rename(
                    columns={
                        "rating_num": f"{temName}POICounts",    # 點位數量
                        "rating_n_1":
                            f"{temName}POIRatingCountsSum",    # google評論數量加總
                        "rating_sta": f"{temName}POIRatingStarSum"
                    },    # google評論星星數加總
                    inplace=True)
            else:
                temDF.rename(
                    columns={
                        "rating_num_count": f"{temName}POICounts",    # 點位數量
                        "rating_num_sum":
                            f"{temName}POIRatingCountsSum",    # google評論數量加總
                        "rating_star_sum": f"{temName}POIRatingStarSum"
                    },    # google評論星星數加總
                    inplace=True)
            try:
                PoiDF = PoiDF.merge(temDF, how='inner', on='gridid')
            except:
                PoiDF = temDF

        print("Done!")

        return PoiDF

    # save preprocessed data
    def saveDF(self, DF):
        DF.to_csv(f"{self.DoneFilePath}/DF.csv", index=False)
        print("Save the dataframe after preprocessing!")

    # pipeline of preprocess
    def run(self):

        # create Base DF with complate gridID, Date, Hour
        GRIDID = list(set(self.gridDF['gridid']))
        HOUR = [x for x in range(24)]
        productElements = list(itertools.product(GRIDID, self.DateList, HOUR))
        KEYS = ['GridID', 'Date', 'Hour']
        DF = pd.DataFrame(productElements, columns=KEYS)
        DF['Weekday'] = pd.to_datetime(DF["Date"], format="%Y%m%d").dt.weekday
        DF['IsWeekend'] = DF['Weekday'].apply(lambda x: 1 if x > 4 else 0)

        ## transaction
        TranOnDF, TranOffDF = self.Transaction()
        TranOnKeys, TranOffKeys = ['gridid', 'on_date', 'on_hour'
                                  ], ['gridid', 'off_date', 'off_hour']
        TranOnDF = TranOnDF[TranOnKeys +
                            ['counts']].groupby(TranOnKeys).sum().reset_index()
        TranOffDF = TranOffDF[TranOffKeys + ['counts']].groupby(
            TranOffKeys).sum().reset_index()

        # OnCounts
        DF = DF.merge(TranOnDF, how='left', left_on=KEYS, right_on=TranOnKeys)
        DF = DF[['GridID', 'Date', 'Hour', 'IsWeekend', 'counts']]
        DF.rename(columns={'counts': 'OnCounts'}, inplace=True)
        temDFCols = list(DF.columns)

        # OffCounts
        DF = DF.merge(TranOffDF, how='left', left_on=KEYS, right_on=TranOffKeys)
        DF = DF[temDFCols + ['counts']]
        DF.rename(columns={'counts': 'OffCounts'}, inplace=True)

        # NetCounts
        DF.fillna(0, inplace=True)
        DF['NetCounts'] = DF['OffCounts'] - DF['OnCounts']

        ## Population
        # Age_15_17_Counts, Age_18_21_Counts, ...,
        # Age_Over65_Counts, Age_Total_Counts,
        # WorkPopulationCounts, LivePopulationCounts, TourPopulationCounts
        PopulationDF = self.Population()
        PopuOffKeys = ['網格編號', '日期', '時間']
        DF = DF.merge(PopulationDF,
                      how='left',
                      left_on=KEYS,
                      right_on=PopuOffKeys)
        DF.drop(PopuOffKeys, axis=1, inplace=True)
        DF.fillna(0, inplace=True)

        ## Traffic MRT
        OverlayExit_DF, populationIN, populationOUT = self.TrafficMRT()
        MRTExitCountsDF = OverlayExit_DF.drop('geometry', axis=1).\
                          groupby('gridid').size().reset_index(name='counts')

        populationIN['日期'] = populationIN['日期'].str.split('-').str.join('')
        populationOUT['日期'] = populationOUT['日期'].str.split('-').str.join('')

        MRTPopulationCols = ['gridid', '日期', '時段']
        MRTPopulationInCountsDF = populationIN[MRTPopulationCols+['人次']].\
                                  groupby(MRTPopulationCols).sum().reset_index()
        MRTPopulationOutCountsDF = populationOUT[MRTPopulationCols+['人次']].\
                                   groupby(MRTPopulationCols).sum().reset_index()

        # MRTExitCounts
        DF = DF.merge(MRTExitCountsDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)
        DF.rename(columns={'counts': 'MRTExitCounts'}, inplace=True)
        DF.fillna(0, inplace=True)
        DF['MRTExitCounts'] = DF['MRTExitCounts'].astype(int)

        # MRTPopulationInCounts
        DF = DF.merge(MRTPopulationInCountsDF,
                      how='left',
                      left_on=KEYS,
                      right_on=MRTPopulationCols)
        DF.drop(MRTPopulationCols, axis=1, inplace=True)
        DF.rename(columns={'人次': 'MRTPopulationInCounts'}, inplace=True)
        DF.fillna(0, inplace=True)
        DF['MRTPopulationInCounts'] = DF['MRTPopulationInCounts'].astype(int)

        # MRTPopulationOutCounts
        DF = DF.merge(MRTPopulationOutCountsDF,
                      how='left',
                      left_on=KEYS,
                      right_on=MRTPopulationCols)
        DF.drop(MRTPopulationCols, axis=1, inplace=True)
        DF.rename(columns={'人次': 'MRTPopulationOutCounts'}, inplace=True)
        DF.fillna(0, inplace=True)
        DF['MRTPopulationOutCounts'] = DF['MRTPopulationOutCounts'].astype(int)

        ## Traffic Bus
        OverlayBusStop_DF = self.TrafficBus()
        BusStopCountsDF = OverlayBusStop_DF.groupby(
            'gridid').size().reset_index(name='counts')

        # BusStopCounts
        DF = DF.merge(BusStopCountsDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)
        DF.rename(columns={'counts': 'BusStopCounts'}, inplace=True)
        DF.fillna(0, inplace=True)
        DF['BusStopCounts'] = DF['BusStopCounts'].astype(int)

        ## Bus Route
        BusRouteDF = self.TrafficBusRoute()

        # BusRouteCounts
        DF = DF.merge(BusRouteDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)

        ## Road SideWalk
        OverlaySideWalk_DF = self.RoadSideWalk()
        SideWalkAreaDF = OverlaySideWalk_DF[[
            'gridid', 'area'
        ]].groupby('gridid').sum().reset_index()

        # SideWalkArea
        DF = DF.merge(SideWalkAreaDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)
        DF.rename(columns={'area': 'SideWalkArea'}, inplace=True)
        DF.fillna(0, inplace=True)

        ## Road MarkedSideWalk
        OverlayMarkedSideWalk_DF = self.RoadMarkedSideWalk()
        MarkedSideWalkAreaDF = OverlayMarkedSideWalk_DF[[
            'gridid', 'area'
        ]].groupby('gridid').sum().reset_index()

        # MarkedSideWalkArea
        DF = DF.merge(MarkedSideWalkAreaDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)
        DF.rename(columns={'area': 'MarkedSideWalkArea'}, inplace=True)
        DF.fillna(0, inplace=True)

        ## Road BikeRoute
        OverlayBikeRoute_DF = self.RoadBikeRoute()
        BikeRouteLengthDF = OverlayBikeRoute_DF[[
            'gridid', 'length'
        ]].groupby('gridid').sum().reset_index()

        # BikeRouteLength
        DF = DF.merge(BikeRouteLengthDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)
        DF.rename(columns={'length': 'BikeRouteLength'}, inplace=True)
        DF.fillna(0, inplace=True)

        ## Road Tree
        OverlayTree_DF = self.RoadTree()
        TreeCountsDF = OverlayTree_DF.groupby('gridid').size().reset_index(
            name='counts')

        # TreeCounts
        DF = DF.merge(TreeCountsDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)
        DF.rename(columns={'counts': 'TreeCounts'}, inplace=True)
        DF.fillna(0, inplace=True)
        DF['TreeCounts'] = DF['TreeCounts'].astype(int)

        ## Road Light
        OverlayLight_DF = self.RoadLight()
        LightCountsDF = OverlayLight_DF.groupby('gridid').size().reset_index(
            name='counts')

        # LightCounts
        DF = DF.merge(LightCountsDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)
        DF.rename(columns={'counts': 'LightCounts'}, inplace=True)
        DF.fillna(0, inplace=True)
        DF['LightCounts'] = DF['LightCounts'].astype(int)

        ## Road Net
        RoadNetDF = self.RoadRoadNet()

        # RoadArea, RoadAreaRatio
        DF = DF.merge(RoadNetDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)

        ## Road Length
        RoadLengthDF = self.RoadRoadLength()

        # ExpressWayLength, ProvincialHighwayLength, UrbanRoad_RoadStreetLength, UrbanRoad_LaneAlleyLength, TotalRoadLength
        DF = DF.merge(RoadLengthDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)

        ## Terrain
        TerrainDF = self.Terrain()

        # Slope, SlopeClass, ELEV_mean, NDVImean, Coverage
        DF = DF.merge(TerrainDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)

        ## Development
        DevelopmentDF = self.Development()

        # commercial_floor_area, commercial_ratio, mixed_floor_area, mixed_ratio, ...
        DF = DF.merge(DevelopmentDF,
                      how='left',
                      left_on='GridID',
                      right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)

        ## Land
        LandDF = self.Land()

        # NatureRatio, CommerceRatio, ResidenceRatio, ...
        DF = DF.merge(LandDF, how='left', left_on='GridID', right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)

        ## POI
        PoiDF = self.POI()

        # TouristAttractionRatingNumCount, TouristAttractionRatingNumSum, TouristAttractionRatingStarSum, ...
        DF = DF.merge(PoiDF, how='left', left_on='GridID', right_on='gridid')
        DF.drop('gridid', axis=1, inplace=True)

        ## save
        self.saveDF(DF)


#         return DF

if __name__ == '__main__':

    p = Preprocess(DateList=['20230305', '20230311', '20230317', '20230322'])
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
