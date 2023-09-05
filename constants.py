LANDUSE_CATEGORY_CODE = {
    "1": "Nature",
    "2": "Commerce",
    "3": "Residence",
    "4": "MixedResidence",
    "5": "Industry",
    "6": "Infrastructure",
    "7": "Education",
    "8": "Leisure",
    "9": "OpenSpace",
    "10": "Traffic",
    "11": "Other",
}

# LANDUSE CODE: LANDUSE CATEGORY CODE
LANDUSE_CATEGORY = {
    "101": "1",  # 水田
    "102": "1",  # 旱田
    "103": "1",  # 果園
    "105": "1",  # 畜牧
    "106": "1",  # 農業相關設施
    "204": "1",  # 混淆林
    "206": "1",  # 其他森林利用土地
    "302": "10",  # 一般鐵路及相關設施
    "303": "10",  # 高速鐵路及相關設施
    "304": "10",  # 捷運及相關設施
    "309": "10",  # 道路相關設施
    "310": "10",  # 港口
    "401": "1",  # 河道
    "408": "1",  # 水利構造物
    "501": "2",  # 商業
    "503": "4",  # 混合使用住宅
    "508": "11",  # 其他建築用地
    "602": "7",  # 學校
    "605": "11",  # 公用設備
    "701": "8",  # 文化設施
    "703": "8",  # 休閒設施
    "903": "1",  # 裸露地
    "905": "1",  # 空置地
    "10400": "1",  # 水產善殖
    "20100": "1",  # 針葉林
    "20200": "1",  # 闊葉林
    "20300": "1",  # 竹林
    "30100": "10",  # 機場
    "30500": "10",  # 國道
    "30600": "10",  # 省道
    "30700": "10",  # 快速公路
    "30800": "10",  # 一般道路
    "40200": "1",  # 堤防
    "40300": "1",  # 溝渠
    "40500": "1",  # 湖泊
    "40600": "1",  # 蓄水池
    "40700": "1",  # 水道沙洲灘地
    "40900": "1",  # 防汛道路
    "50200": "3",  # 純住宅
    "50400": "5",  # 製造業
    "50500": "5",  # 倉儲
    "50600": "8",  # 宗教
    "50700": "11",  # 殯葬設旄
    "60100": "6",  # 政府機關
    "60300": "6",  # 醫療保健
    "60400": "6",  # 社會福利設施
    "60600": "11",  # 環保設施
    "70103": "8",  # 其他文化設施
    "70200": "9",  # 公園綠地廣場
    "80100": "1",  # 礦業及相關設施
    "80200": "1",  # 土石及相關設施
    "90100": "1",  # 濕地
    "90200": "1",  # 草生地
    "90400": "1",  # 營建剩餘土石收容處理相關設施
}

BUILDING_CATEGORY = {
    50200: "residential",  # 純住宅
    503: "mixed_residential",  # 混合使用住宅
    501: "commercial",  # 商業
}

ROAD_CATEGORY = {
    "ExpressWayLength": ["94212"],  # 一般快速道路(含市區高架道路)
    "ProvincialHighwayLength": ["94213", "94213a"],  # 省道  # 省道快速公路
    "UrbanRoad_RoadStreetLength": ["94214"],  # 市區道路(路、街)
    "UrbanRoad_LaneAlleyLength": [
        "94214b",  # 市區道路(巷、弄)
        "94214c",  # 區塊道路、專用道路
        "94216",  # 鄉(鎮)道路
        "94219",  # 產業道路
    ],
}


NDVI_THRESHOLD = 0.4
POI_FILTER_NUM = 5
