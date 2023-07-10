# -*- coding: utf-8 -*-
"""
許可なく転載を禁止します。また、商用利用を禁止します。
AMD_Tools3.py
    メッシュ気象データの利用に必要な関数(計算で使う道具)のコレクション。
    1. GetMetData：GetMetData_Tile(isTile=True)またはGetMetData_Area(isTile=Falseまたは無)の呼び出し関数。
    2. GetSceData：GetSceData_Tile(isTile=True)またはGetSceData_Area(isTile=Falseまたは無)の呼び出し関数。
    3. GetGeoData：GetGeoData_Tile(isTile=True)またはGetGeoData_Area(isTile=Falseまたは無)の呼び出し関数。
    4. GetCSV_Table：CSV形式の表を読み込み、文字列のリストにする関数。
    5. GetCSV_Map：CSV形式のテキストファイルを数値として配列変数に読み込む関数。
    6. PutCSV_TS：時系列のデータをCSV形式のファイルで出力する関数。
    7. PutCSV_Map：2次元の浮動小数点配列をCSV形式のファイルで出力する関数。
    8. PutCSV_MT:3次元の配列を、3次メッシュコードをキーとするテーブルの形式のCSVファイルで出力する関数。
    9. PutNC_Map：2次元(空間分布)の気象変量をnetCDF形式のファイルで出力する関数。
   10. PutNC_3D：3次元(空間分布×時間変化)の気象変量をnetCDF形式のファイルで出力する関数。
   11. PutGSI_Map：2次元(空間分布)の配列を地理院地図用HTMLで出力する関数。
   12. mesh2lalo:緯度・経度を3次メッシュコードに変換する関数
   13. lalo2mesh:3次メッシュコード(文字列)を緯度・経度に変換する関数
   14. timedom:日付文字列の配列[初日,最終日]から、この期間のdatetimeオブジェクトの配列を返す関数。
   15. lalodom:数値の配列[緯度,緯度, 経度,経度]から、この区間を含むメッシュ範囲の中心緯度の配列を返す関数。
   16. GetMetData_Tile：メッシュ農業気象データを取得する関数(1次メッシュ区切り対応版)。
   17. GetSceData_Tile：メッシュ温暖化シナリオデータを取得する関数(1次メッシュ区切り対応版)。
   18. GetGeoData_Tile：土地利用や都道府県域などの地理情報を取得する関数(1次メッシュ区切り対応版)。
   19. GetMetData_Area：メッシュ農業気象データを取得する関数(Area区切り対応版)。
   20. GetSceData_Area：メッシュ温暖化シナリオデータを取得する関数(Area区切り対応版)。
   21. GetGeoData_Area：土地利用や都道府県域などの地理情報を取得する関数(Area区切り対応版)。
   22. GetMetDataHourly:メッシュ農業気象データ（時別値）を取得する関数

    改変履歴：
    20201121 GetMetDataHourly 関数を追加 （時別値取得用関数）
    20190221 GetCSVのエンコード問題に対応、１次メッシュ区分に対応
    20180913 GetGeoDataのバグを修正
    20180628 プロキシーサーバーに対応
    20180405 IDパスワード認証に対応
    20171204 Matplotlib2に対応
    20171125 PutGSI_Map関数を追加
    20170603 lonrange関数他を追加
    20170518 GetCSV_List関数の追加
    20170502 コメント文の修正
    20170425 PutKMZ関数の改良
    20170208 Python3バージョン
    20140129 PurCSV_MTに機能を追加
    20131118 3次メッシュコードをキーとする表として出力する関数(PurCSV_MT)の追加
    20130314 地理情報を読み取る関数の追加
    20121205 安定版初版
    Copyright (C)  OHNO, Hiroyuki
"""
#_引用符の中を 通知された認証情報で書き換えてください。___________
USER="ogawa-hiro2023"
PASSWORDS=["kiPNitXR2023", ""] #二つのうち、どちらかが正しければデータは取得できます

#_ プロキシーサーバー経由で接続する方は下も設定してください。______
#　（使用しない場合はこのままにしてください）
PROXY_IP = ""  #プロキシーサーバーのIPアドレス(文字列で)
PROXY_PORT = ""  #使用するポート番号(文字列で)

#_______________________________________________________
#///////////// 以下には変更を加えないでください ////////////////////////
from sys import exit
from os import unlink #,fdopen
from os.path import join,exists,isdir,basename
from datetime import datetime as dt
from datetime import timedelta as td
from math import floor, ceil
import numpy as np
import numpy.ma as ma #PutKMZ_Mapで使ってる
import tempfile
from random import randint
from netCDF4 import date2num, num2date, Dataset
import codecs
import urllib
import urllib.request
import ssl
import copy
import xarray as xr

# リモートで動かすときに必要
#import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import DateFormatter,DayLocator
from threading import Thread

TIMEZERO = dt.strptime("1900-01-01","%Y-%m-%d")
ssl._create_default_https_context = ssl._create_unverified_context

MESHLIST = ["3622","3623","3624","3631","3724","3725",
        "3823","3824","3831","3926","3927",
        "3928","4027","4028","4128",
        "4129","4229","4230","4328","4329","4429",
        "4529","4530","4531","4629","4630",
        "4631","4728","4729","4730","4731",
        "4828","4829","4830","4831","4839","4928","4929",
        "4930","4931","4932","4933","4934","4939","5029",
        "5030","5031","5032","5033","5034","5035","5036",
        "5038","5039","5129","5130","5131","5132","5133",
        "5134","5135","5136","5137","5138","5139","5229",
        "5231","5232","5233","5234","5235","5236","5237",
        "5238","5239","5240","5332","5333","5334","5335",
        "5336","5337","5338","5339","5340","5432","5433",
        "5435","5436","5437","5438","5439","5440",
        "5536","5537","5538","5539","5540","5541","5636",
        "5637","5638","5639","5640","5641","5738","5739",
        "5740","5741","5839","5840","5841","5939","5940",
        "5941","5942","6039","6040","6041","6139","6140",
        "6141","6239","6240","6241","6243","6339","6340",
        "6341","6342","6343","6439","6440","6441","6442",
        "6443","6444","6445","6540","6541","6542","6543",
        "6544","6545","6641","6642","6643","6644","6645",
        "6741","6742","6840","6841","6842"]

def check_user(error=False):

    if USER == "利用者ID":
        if error:
            print("")
            print("  ====> ファイルAMD_Tools.pyの49-50行目に利用者IDとパスワードを指定してください。 <====")
            print("")
            exit(1)
        else:
            print("")
            print("  ====> WARNING: userID and password are not set <====")
            print("  ====> not all data sources can be used         <====")
            print("")

check_user()

def urljoin(xs):
    if len(xs) <= 1:
        return "".join(xs)
    if xs[0].startswith("http"):
        return xs[0].rstrip("/") + "/" + "/".join(xs[1:-1]) + "/" + xs[-1].lstrip("/")
    else:
        return join(*xs)

def ir(x): return int(round(x))

def nan2mv(a,val):
    a[a!=a] = val

def mv2nan(a,val):
    a[a==val] = np.nan

def ma2nan(a):
    if a.mask is not False:
        a[a.mask]=np.nan
    return np.array(a.data)

def lalo2mesh(lat,lon):
    lat = lat * 1.5
    lon = lon - 100
    lat1 = int(floor(lat))
    lat = 8*(lat - lat1)
    lon1 = int(floor(lon))
    lon = 8*(lon - lon1)
    lat2 = int(floor(lat))
    lat = 10*(lat - lat2)
    lon2 = int(floor(lon))
    lon = 10*(lon - lon2)
    lat3 = int(floor(lat))
    lon3 = int(floor(lon))
    return "".join([str(x) for x in [lat1,lon1,lat2,lon2,lat3,lon3]])

def mesh2lalo(code):
    assert len(code) == 8
    lat = int(code[:2])/1.5 + int(code[4])/12.0 + int(code[6])/120.0
    lon = int(code[2:4]) + 100 + int(code[5])/8.0 + int(code[7])/80.0
    return lat+1/240.0,lon+1/160.0


def timedom(tup):
    t1 = dt.strptime(tup[0], '%Y-%m-%d')
    t2 = dt.strptime(tup[1], '%Y-%m-%d')
    noda = (t2 - t1).days
    tr = [t1+td(days=oo) for oo in range(noda+1)]
    return np.array(tr)


def lalodom(tup):
    assert tup[0] < tup[1] and tup[2] < tup[3]
    div = 120.0
    nodi = floor(tup[1]*div) - floor(tup[0]*div)
    deg0 = floor(tup[0]*div)/div + 0.5/div
    lat = [deg0+oo/div for oo in range(nodi+1)]
    div = 80.0
    nodi = floor(tup[3]*div) - floor(tup[2]*div)
    deg0 = floor(tup[2]*div)/div + 0.5/div
    lon = [deg0+oo/div for oo in range(nodi+1)]
    return np.array(lat), np.array(lon)


def getFileContent(path):
    read = None
    for e in ['utf-8_sig','utf-8','cp932','euc-jp']:
        try:
            with open(path,'r',encoding=e) as f:
                read = f.read()
                break
        except UnicodeDecodeError:
            pass
    if read is None:
        print("CSVファイルを読み込めませんでした。ファイルのエンコーディングを確認してください。",path)
        exit(1)
    else:
        return read


def StartUnlink(path):
    t = Thread(target=UnlinkTryLoop, args=(path,))
    t.start()


def UnlinkTryLoop(path):
    while True:
        try:
            unlink(path)
            break
        except:
            continue


# これらは近く廃止＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
def timrange(tim1, tim2):
    t1 = dt.strptime(tim1, '%Y-%m-%d')
    t2 = dt.strptime(tim2, '%Y-%m-%d')
    noda = (t2 - t1).days
    tr = [t1+td(days=oo) for oo in range(noda+1)]
    return np.array(tr)

def latrange(deg1, deg2):
    assert deg1 < deg2
    div = 120.0
    nodi = floor(deg2*div) - floor(deg1*div)
    deg0 = floor(deg1*div)/div + 0.5/div
    deg = [deg0+oo/div for oo in range(nodi+1)]
    return np.array(deg)

def lonrange(deg1, deg2):
    assert deg1 < deg2
    div = 80.0
    nodi = floor(deg2*div) - floor(deg1*div)
    deg0 = floor(deg1*div)/div + 0.5/div
    deg = [deg0+oo/div for oo in range(nodi+1)]
    return np.array(deg)
# これらは近く廃止＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝


def tll_extract(dh,tmd,lld,element):
    time  = dh.variables['time']
    times = num2date(time[:], units=time.units)
    tr = tmd.restrict(times)
    #tim = times[tr]
    tim = np.array([dt(x.year,x.month,x.day) for x in times[tr]])
    latitude = dh.variables['lat']
    yr = lld.latrestrict(latitude[:])
    lat = latitude[yr]
    longitude = dh.variables['lon']
    xr = lld.lonrestrict(longitude[:])
    lon = longitude[xr]
    ncMet =  dh.variables[element]
    name = ncMet.long_name
    unit = ncMet.units
    dims = ncMet.dimensions
    fill = ncMet._FillValue
    if dims == ('time', 'lat', 'lon'):
        try:
            Me = np.array(ncMet[tr, yr, xr])
        except IndexError:
            Me = np.array([])
    else:
        vals = ncMet[:]
        tidx = dims.index("time")
        if tidx != 0:
            vals = np.swapaxes(vals,0,tidx)
            if tidx == 1:
                dims = (dims[1],dims[0],dims[2])
            else:
                dims = (dims[2],dims[1],dims[0])
        tidx = dims.index("lat")
        if tidx != 1:
            vals = np.swapaxes(vals,1,2)
        Me = vals[tr,:,:][:,yr,:][:,:,xr]
    Met = np.where(Me == fill, np.nan, Me)
    dh.close()
    #print("LAT0",lat)
    #print("MET0",Met[:3,:3,:3])
    if len(lat) and lat[0] > lat[-1]:
        lat = lat[::-1]
        Met = Met[:,::-1,:]
    #print("LAT1",lat)
    #print("MET1",Met[:3,-3:,:3])
    return tim,lat,lon,Met,name,unit

def xtll_extract(dh,tmd,lld,element):
    time  = dh.variables['time']
    #times = np.array([dt.fromtimestamp(x.item()/1000000000) for x in time.values])
    times = np.array([dt.utcfromtimestamp((x - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')) for x in time.values])
    tr = tmd.restrict(times)
    tim = times[tr]
    latitude = dh.variables['lat']
    yr = lld.latrestrict(latitude[:])
    lat = latitude[yr]
    longitude = dh.variables['lon']
    xxr = lld.lonrestrict(longitude[:])
    lon = longitude[xxr]
    ncMet = dh.variables[element]
    name = ncMet.attrs["long_name"]
    unit = ncMet.attrs["units"]
    try:
        Met = np.array(ncMet[tr, yr, xxr])
    except IndexError:
        Met = np.array([])
    dh.close()
    if len(lat) and lat[0] > lat[-1]:
        lat = lat[::-1]
        Met = Met[:,::-1,:]
    return tim,ma.masked_array(lat.values),ma.masked_array(lon.values),Met,name,unit

class Area:
    latIdx = {1:799, 2:879, 3:799, 4:479, 5:799, 6:799} # number of cell for each Areas
    lonIdx = {1:559, 2:479, 3:559, 4:639, 5:399, 6:799} # number of cell for each Areas
    def __init__(self,name,num,s,n,w,e):
        self.name = name
        self.num = num
        self.s,self.n,self.w,self.e = s,n,w,e
    def __str__(self):
        return self.name
    def __contains__(self, ll):
        return self.s <= ll.latmin and self.n >= ll.latmax and self.w <= ll.lonmin and self.e >= ll.lonmax
    def get_idx(self, ll):
        if self.s <= ll.latmin and ll.latmin <= self.n:
            sIdx = int((ll.latmin - self.s)/(2/3/8/10))
        else:
            sIdx = 0
        if self.s <= ll.latmax and ll.latmax <= self.n:
            nIdx = int((ll.latmax - self.s)/(2/3/8/10))
        else:
            nIdx = self.latIdx[self.num]
        if self.w <= ll.lonmin and ll.lonmin <= self.e:
            wIdx = int((ll.lonmin - self.w)/(1/8/10))
        else:
            wIdx = 0
        if self.w <= ll.lonmax and ll.lonmax <= self.e:
            eIdx = int((ll.lonmax - self.w)/(1/8/10))
        else:
            eIdx = self.lonIdx[self.num]
        return f'[{sIdx}:1:{nIdx}][{wIdx}:1:{eIdx}]'

AREAS = {
    "北海道" : Area("北海道",1, 118/3, 46.0, 139.0, 146.0),
    "東北" : Area("東北",2, 208/6, 42.0, 137.0, 143.0),
    "関東北陸" : Area("関東北陸",3, 32.0, 232/6, 135.0, 142.0),
    "西日本" : Area("西日本",4, 196/6, 220/6, 130.0, 138.0),
    "九州" : Area("九州",5, 172/6, 106/3, 128.0, 133.0),
    "西南諸島" : Area("西南諸島",6, 24.0, 88/3, 122.0, 132.0)
}

class LatLonDomain:
    def __init__(self,latmin,latmax,lonmin,lonmax,area=None):
        """2d-region: latmin,latmax,lonmin,lonmax"""
        self.latmin = latmin
        self.latmax = latmax
        self.lonmin = lonmin
        self.lonmax = lonmax
        self.area = area
        self.check()

    def __str__(self):
        return str((self.latmin,self.latmax,self.lonmin,self.lonmax))

    def check(self):
        if self.latmin > self.latmax:
            raise ValueError("South:" + str(self.latmin) + " North:" + str(self.latmax))
        if self.lonmin > self.lonmax:
            raise ValueError("West:" + str(self.lonmin) + " East:" + str(self.lonmax))
        if self.area is not None:
            if self.area not in ['Area1','Area2','Area3','Area4','Area5','Area6']:
                raise ValueError(f'Unknown area name: {self.area}')
            self.area = [x for x in AREAS.values() if x.num==int(self.area[-1])][0]

    def get_area(self,areas=None):
        if areas is None:
            areas = AREAS
        matches = [a for a in areas.values() if self in a]
        if not matches:
            raise ValueError("No area containing " +str(self) + " found.")
        self.area = matches[0]
        return "Area"+str(self.area.num)

    def latrestrict(self,a):
        if self.latmin != self.latmax:
            b = (a >= self.latmin) & (a <= self.latmax)
        else:
            c = np.abs(a-self.latmin)
            v = np.min(c)
            b = (c == v)
            for i in range(len(b)):
                if b[i]:
                    if i < len(b)-1:
                        b[i+1] = False
                    break
        return b

    def lonrestrict(self,a):
        if self.lonmin != self.lonmax:
            b = (a >= self.lonmin) & (a <= self.lonmax)
        else:
            c = np.abs(a-self.lonmin)
            v = np.min(c)
            b = (c == v)
            for i in range(len(b)):
                if b[i]:
                    if i < len(b)-1:
                        b[i+1] = False
                    break
        return b

    def geogrid(self):
        return ",".join([str(x) for x in [self.latmax,self.lonmin,self.latmin,self.lonmax]])

    def getIdx(self):
        if not self.area:
            self.get_area()
        return self.area.get_idx(self)

    def codes(self):
        lats = [str(x) for x in range(floor(self.latmin*3/2),ceil(self.latmax*3/2))]
        lons = [str(x) for x in range(floor(self.lonmin - 100),ceil(self.lonmax - 100))]
        if len(lats) == 0:
            lats = [str(floor(self.latmin*3/2))]
        if len(lons) == 0:
            lons = [str(floor(self.lonmin-100))]
        return [f"{lat}{lon}" for lat in lats for lon in lons if lat+lon in MESHLIST]

    def getCodeWithIdx(self):
        ret = []
        codes = self.codes()
        for c in codes:
            ret.append((c, self.getCodeIdx(c)))
        return ret

    def getCodeIdx(self, code):
        f,b = int(code[0:2]), int(code[2:])
        s = f/1.5
        n = s + 40/60
        w = 100 + b
        e = w + 1
        sIdx = self.calcCodeIdx(s, n, self.latmin, True, False)
        nIdx = self.calcCodeIdx(s, n, self.latmax, True, True)
        wIdx = self.calcCodeIdx(w, e, self.lonmin, False, False)
        eIdx = self.calcCodeIdx(w, e, self.lonmax, False, True)
        return f'[{sIdx}:1:{nIdx}][{wIdx}:1:{eIdx}]'

    def calcCodeIdx(self, vmin, vmax, v, lat=True, end=True):
        if lat:
            div = 40/60/80
        else:
            div = 1/80
        if vmin <= v <= vmax:
            if end:
                idx = ceil((v - vmin)/div)
            else:
                idx = floor((v - vmin)/div)
            if idx < 0: return 0
            if idx > 79: return 79
            return idx
        else:
            if end:
                return 79
            else:
                return 0

class TimeDomain:
    def __init__(self,t0,t1):
        """time range, t0,t1 dates in yyyy-mm-dd format"""
        if "-" in t0:
            self.beg = dt.strptime(t0,'%Y-%m-%d')
        elif "." in t0:
            self.beg = dt.strptime(t0,'%Y.%m.%d')
        elif "/" in t0:
            self.beg = dt.strptime(t0,'%Y/%m/%d')
        elif " " in t0:
            self.beg = dt.strptime(t0,'%Y %m %d')
        if "-" in t1:
            self.end = dt.strptime(t1,'%Y-%m-%d')
        elif "." in t1:
            self.end = dt.strptime(t1,'%Y.%m.%d')
        elif "/" in t1:
            self.end = dt.strptime(t1,'%Y/%m/%d')
        elif " " in t1:
            self.end = dt.strptime(t1,'%Y %m %d')
    def years(self):
        return self.end.year - self.beg.year + 1
    def yrange(self):
        return range(self.beg.year, self.end.year + 1)
    def restrict(self,a):
        b = (a >= self.beg) & (a <= self.end)
        return b

    def geogrid(self):
        a = (self.beg - TIMEZERO).days - 1
        b = (self.end - TIMEZERO).days + 1
        return '"' + str(a) + '&lt;time","time&lt;'+ str(b) + '"'

    def isleap(self, year):
        return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

    def getIdx(self):
        ret = []
        for y in range(self.beg.year, self.end.year+1):
            if y == self.beg.year:
                sidx = (self.beg - dt(y,1,1)).days
            else:
                sidx = 0
            if y == self.end.year:
                eidx = (self.end - dt(y,1,1)).days
            else:
                eidx = 365 if self.isleap(y) else 364
            ret.append((y,f'[{sidx}:1:{eidx}]'))
        return ret

class TimeDomainHourly(TimeDomain):
    def __init__(self,t0,t1):
        """time range for hourly data, t0,t1 datetime in %Y-%m-%dT%H:%M:%S format"""
        self.beg = self.parseDT(t0)
        if t0 == t1:
            self.end = self.parseDT(t1+"T24", True)
        else:
            self.end = self.parseDT(t1, True)

    def parseDT(self, t, isEnd=False):
        try:
            if "T" in t:
                x = t.split('T')
                d = self.parseD(x[0])
                return self.parseT(d, x[1], isEnd)
            else:
                d = self.parseD(t)
                return self.parseT(d, None, isEnd)
        except (NameError,ValueError):
            raise("Cannot parse TimeDomain.")

    def parseD(self, d):
        if "-" in d:
            dd = dt.strptime(d,'%Y-%m-%d')
        elif "." in d:
            dd = dt.strptime(d,'%Y.%m.%d')
        elif "/" in d:
            dd = dt.strptime(d,'%Y/%m/%d')
        elif " " in d:
            dd = dt.strptime(d,'%Y %m %d')
        return dd

    def parseT(self, d, t, isEnd=False):
        delta = td()
        if t is None:
            if isEnd:
                hour = 23
                delta += td(hours=1)
            else:
                hour = 1
        elif ":" in t:
            tt = t.split(":")
            if int(tt[1]) > 30:
                delta += td(hours=1)
            hour = int(tt[0])
        else:
            hour = int(t)
        if hour == 24:
            hour = 23
            delta += td(hours=1)
        return dt(d.year, d.month, d.day, hour) + delta

    def geogrid(self):
        a = (self.beg - TIMEZERO).days * 24 + self.beg.hour - 1
        b = (self.end - TIMEZERO).days * 24 + self.end.hour + 1
        return '"' + str(a) + '&lt;time","time&lt;'+ str(b) + '"'

    def getIdx(self):
        ret = []
        for y in range(self.beg.year, self.end.year+1):
            if y == self.beg.year:
                sidx = int((self.beg - dt(y,1,1)) / td(hours=1)) - 1
            else:
                sidx = 0
            if y == self.end.year:
                eidx = int((self.end - dt(y,1,1)) / td(hours=1))
            else:
                eidx = 8783 if self.isleap(y) else 8759
            ret.append((y,f'[{sidx}:1:{eidx}]'))
        return ret

def url2dh(url):
    url = url.replace("\\","/")
    a,b = url.split("?")
    no_pw = not url.startswith("https://amd.rd.naro.go.jp")

    if no_pw:
        try: #try local storage
            dh = Dataset(a)
            return dh,None
        except:
            pass

        try: #try local storage
            dh = Dataset(a[:-3])
            return dh,None
        except:
            pass

        try: #try old passwordless site
            dh = Dataset(a[:-3] + "?" + b)
            return dh,None
        except:
            pass
    else:
        check_user(error=True)
        for pw in PASSWORDS:
            #print("accessing URL",url)
            q = urllib.request.Request(url)
            q.add_header('User-Agent','curl/7.50.1')
            q.add_header('Accept','*/*')
            mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()  #
            #print("USER",USER,"PW",pw)
            mgr.add_password(None, url, USER, pw)

            auth_handler = urllib.request.HTTPBasicAuthHandler(mgr)
            if PROXY_IP:
                print(f"using proxy IP: {PROXY_IP} port:{PROXY_PORT}")
                proxy = urllib.request.ProxyHandler({'https':PROXY_IP+':'+PROXY_PORT})
                opener = urllib.request.build_opener(proxy,auth_handler)
            else:
                opener = urllib.request.build_opener(auth_handler)
            urllib.request.install_opener(opener)
            try:
                response = urllib.request.urlopen(q)
            except urllib.error.HTTPError:
                #print("WRONG PW",pw)
                continue
                #raise
            data = response.read()
            d = tempfile.gettempdir()
            p = "amd_" + str(randint(100000000,999999999))
            dfile = join(d,p)
            f = open(dfile,"wb")
            f.write(data)
            f.close()
            dh = Dataset(dfile)
            return dh,dfile

    raise ValueError("Network error found. Please check id and password or try again after a while. ",a)


def xurl2dh(url):
    url = url.replace("\\","/")
    a,b = url.split("?")
    no_pw = not url.startswith("https://amd.rd.naro.go.jp")

    if no_pw:
        try: #try local storage
            dh = xr.open_dataset(a, decode_timedelta=False) # masks fillvalues (see mask_and_scale)
            return dh,None
        except:
            pass

        try: #try local storage
            #print("trying",a[:-3])
            dh = xr.open_dataset(a[:-3], decode_timedelta=False)
            return dh,None
        except:
            pass

        try: #try old passwordless site
            dh = xr.open_dataset(a[:-3] + "?" + b, decode_timedelta=False)
            return dh,None
        except:
            pass
    else:
        check_user(error=True)
        for pw in PASSWORDS:
            #print("accessing URL",url)
            q = urllib.request.Request(url)
            q.add_header('User-Agent','curl/7.50.1')
            q.add_header('Accept','*/*')
            mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()  #
            #print("USER",USER,"PW",pw)
            mgr.add_password(None, url, USER, pw)

            auth_handler = urllib.request.HTTPBasicAuthHandler(mgr)
            if PROXY_IP:
                print(f"using proxy IP: {PROXY_IP} port:{PROXY_PORT}")
                proxy = urllib.request.ProxyHandler({'https':PROXY_IP+':'+PROXY_PORT})
                opener = urllib.request.build_opener(proxy,auth_handler)
            else:
                opener = urllib.request.build_opener(auth_handler)
            urllib.request.install_opener(opener)
            try:
                response = urllib.request.urlopen(q)
            except urllib.error.HTTPError:
                #print("WRONG PW",pw)
                continue
                #raise
            data = response.read()
            d = tempfile.gettempdir()
            p = "amd_" + str(randint(100000000,999999999))
            dfile = join(d,p)
            f = open(dfile,"wb")
            f.write(data)
            f.close()
            #dh = xr.open_dataset(dfile)
            with xr.open_dataset(dfile, decode_timedelta=False) as dh:
                pass
            return dh,dfile
    raise ValueError("Network error found. Please check id and password or try again after a while. ",a)



def GetMetData_Area(element, timedomain, lalodomain, area=None,
               cli=False, namuni=False, url='https://amd.rd.naro.go.jp/opendap/AMD/'):
    """
概要：
    メッシュ農業気象データを、気象データをデータ配信サーバーまたはローカルファイルから取得する関数。
書式：
　GetMetData(element, timedomain, lalodomain, area=None, cli=False, namuni=False, url='https://amd.rd.naro.go.jp/opendap/AMD')
引数(必須)：
    element：気象要素記号で、'TMP_mea'などの文字列で与える
    timedomain：取得するデータの時間範囲で、['2008-05-05', '2008-05-05']
                のような文字列の2要素リストで与える。特定の日のデータを
                取得するときは、二カ所に同じ日付を与える。
    lalodomain：取得するデータの緯度と経度の範囲で、
                [36.0, 40.0, 130.0, 135.0] のように緯度,緯度,経度,経度の順で指定する。
                特定地点のデータを取得するときは、緯度と経度にそれぞれ同じ値を与える。
引数(必要に応じ指定)：
    cli:True => 平年値が返される。
        False => 観測値が返される。
    namuni:True => 気象要素の正式名称と単位を取り出す。戻り値の数は2つ増えて6つになる。
        False => 気象要素の正式名称を取り出さない。戻り値の数は4つ(気象値、時刻、緯度、経度)。
    area:データを読み出すエリア(Area1～Area6)を指定する。省略した場合は自動的に選ばれる。
    url:データファイルの場所を指定する。省略した場合はデータ配信サーバーに読みに行く。
        ローカルにあるファイルを指定するときは、AreaN(N=1～6)の直上(通常は"・・・/AMD")を指定する。
戻り値：
    第1戻り値：指定した気象要素の三次元データ。[時刻、緯度、経度]の次元を持つ。
    第2戻り値：切り出した気象データの時刻の並び。Pythonの時刻オブジェクトの
                一次元配列である。
    第3戻り値：切り出した気象データの緯度の並び。実数の一次元配列である。
    第4戻り値：切り出した気象データの経度の並び。実数の一次元配列である。
　 第5戻り値(namuni=Trueのときのみ)：気象データの正式名称。文字列である。
    第6戻り値(namuni=Trueのときのみ)：気象データの単位。文字列である。
使用例：北緯35度、東経135度の地点の2008年1月1日～2012年12月31日の日最高気温を取得する場合。
    import AMD_Tools3 as AMD
    timedomain = ['2008-01-01', '2012-12-31']
    lalodomain = [35.0,  35.0, 135.0, 135.0]
    Tm, tim, lat, lon = AMD.GetMetData('TMP_max', timedomain, lalodomain)
    """

    lld = LatLonDomain(*lalodomain)
    if area is None:
        area = lld.get_area()
    td = TimeDomain(*timedomain)
    filename = 'AMD_' + area +'_' + ('Cli_' if cli else '') + element + '.nc.nc'
    dhpaths = []

    years = td.getIdx()
    tim,lat,lon,Met,name,unit = None,None,None,None,None,None
    for year, tidx in years:
        opendap_source = urljoin([url,area,str(year),filename]) + '?' + element + tidx + lld.getIdx()
        dh,dhpath = url2dh(opendap_source)
        if dhpath is not None:
            dhpaths.append(dhpath)
        ti,la,lo,nm,name,unit = tll_extract(dh,td,lld,element)
        if tim is None:
                tim = ti
                lat = la
                lon = lo
                Met = nm
        else:
            tim = np.concatenate((tim,ti))
            Met = np.concatenate((Met,nm))

    print(('Cli_' if cli else '')+element, Met.shape)
    for dhpath in dhpaths:
        StartUnlink(dhpath)

    if namuni:
        return Met, tim, lat, lon, name, unit
    else:
        return Met, tim, lat, lon

def GetSceData_Area(element, timedomain, lalodomain, model, scenam, area=None,
               cli=False, namuni=False, url='https://amd.rd.naro.go.jp/opendap/AMS'):
    """
概要：
    気候予測シナリオデータを、気象データをデータ配信サーバーまたはローカルファイルから取得する関数。
書式：
　GetSceData(element, timedomain, lalodomain, model, scenam, area=None, cli=False, namuni=False, url='https://amd.rd.naro.go.jp/opendap/AMS')
引数(必須)：
    element：気象要素記号で、'TMP_mea'などの文字列で与える
    timedomain：取得するデータの時間範囲で、['2008-05-05', '2008-05-05']
                のような文字列の2要素リストで与える。特定の日のデータを
                取得するときは、二カ所に同じ日付を与える。
    lalodomain：取得するデータの緯度と経度の範囲で、
                [36.0, 40.0, 130.0, 135.0] のように緯度,緯度,経度,経度の順で指定する。
                特定地点のデータを取得するときは、緯度と経度にそれぞれ同じ値を与える。
    model：気候モデルの記号で、'MIROC5'などの文字列で与える
    scenam：排出シナリオ等の記号で、'RCP8.5'などの文字列で与える
引数(必要に応じ指定)：
    cli:True => 平年値が返される。
        False => 観測値が返される。
    namuni:True => 気象要素の正式名称と単位を取り出す。戻り値の数は2つ増えて6つになる。
        False => 気象要素の正式名称を取り出さない。戻り値の数は4つ(気象値、時刻、緯度、経度)。
    area:データを読み出すエリア(Area1〜Area6)を指定する。省略した場合は自動的に選ばれる。
    url:データファイルの場所を指定する。省略した場合はデータ配信サーバーに読みに行く。
        ローカルにあるファイルを指定するときは、AreaN(N=1〜6)の直上(通常は"・・・/AMS")を指定する。
戻り値：
    第1戻り値：指定した気象要素の三次元データ。[時刻、緯度、経度]の次元を持つ。
    第2戻り値：切り出した気象データの時刻の並び。Pythonの時刻オブジェクトの
                一次元配列である。
    第3戻り値：切り出した気象データの緯度の並び。実数の一次元配列である。
    第4戻り値：切り出した気象データの経度の並び。実数の一次元配列である。
　 第5戻り値(namuni=Trueのときのみ)：気象データの正式名称。文字列である。
    第6戻り値(namuni=Trueのときのみ)：気象データの単位。文字列である。
使用例：MIROC5モデルで予測したRCP8.5シナリオにおける、北緯35度、東経135度の地点の
　　　2020年〜2030年の日最高気温を取得する場合。
    import AMD_Tools3 as AMD
    model  = 'MROC5'
    scenario = 'RCP8.5'
    timedomain = ['2020-01-01', '2030-12-31']
    lalodomain = [35.0,  35.0, 135.0, 135.0]

    Tm, tim, lat, lon = AMD.GetSceData('TMP_max', timedomain, lalodomain,model, scenario)
    """
    lld = LatLonDomain(*lalodomain)
    if area is None:
        area = lld.get_area()
    td = TimeDomain(*timedomain)
    filename = 'AMS_' + area +'_' + ('Cli_' if cli else '') + element + '.nc.nc'
    dhpaths = []
    years = td.getIdx()
    tim,lat,lon,Met,name,unit = None,None,None,None,None,None
    for year, tidx in years:
        opendap_source = urljoin([url,model,scenam,area,str(year),filename]) + '?' + element + tidx + lld.getIdx()
        dh,dhpath = url2dh(opendap_source)
        if dhpath is not None:
            dhpaths.append(dhpath)
        ti,la,lo,nm,name,unit = tll_extract(dh,td,lld,element)
        if tim is None:
            tim = ti
            lat = la
            lon = lo
            Met = nm
        else:
            tim = np.concatenate((tim,ti))
            Met = np.concatenate((Met,nm))
    print(element, Met.shape)
    for dhpath in dhpaths:
        StartUnlink(dhpath)
    if namuni:
        return Met, tim, lat, lon, name, unit
    else:
        return Met, tim, lat, lon

def GetGeoData_Area(element, lalodomain, area=None,
               namuni=False, url='https://amd.rd.naro.go.jp/opendap/AMD/'):
    """
概要：
    土地利用区分等の地理情報をデータ配信サーバーまたはローカルファイルから取得する関数。
書式：
    GetGeoData(element, lalodomain, area=None, namuni=False, url='https://amd.rd.naro.go.jp/opendap/AMD/')
引数(必須)：
    element：地理情報記号で、'altitude'などの文字列で与える
    lalodomain：取得するデータの緯度と経度の範囲で、
        [36.0, 40.0, 130.0, 135.0] のように緯度,緯度,経度,経度の順で指定する。
        特定地点のデータを取得するときは、緯度と経度にそれぞれ同じ値を与える。
引数(必要に応じ指定)：
    namuni:True => 気象要素の正式名称と単位を取り出す。戻り値の数は2つ増えて6つになる。
        False => 気象要素の正式名称を取り出さない。戻り値の数は4つ(気象値、時刻、緯度、経度)。
    area:データを読み出すエリア(Area1～Area6)を指定する。省略した場合は自動的に選ばれる。
    url:データファイルの場所を指定する。省略した場合はデータ配信サーバーに読みに行く。
        ローカルにあるファイルを指定するときは、AreaN(N=1～6)の直上(通常は"・・・/AMD")を指定する。
戻り値：
    第1戻り値：指定した地理情報値の二次元データ。[緯度、経度]の次元を持つ。
    第2戻り値：切り出した地理情報値の緯度の並び。実数の一次元配列である。
    第3戻り値：切り出した地理情報値の経度の並び。実数の一次元配列である。
    第4戻り値(namuni=Trueのときのみ)：地理情報の正式名称。文字列である。
    第5戻り値(namuni=Trueのときのみ)：地理情報の単位。文字列である。

使用例：
    北緯35～36、東経135～136度の範囲にある各メッシュの水田面積比率の分布を取得する場合。
    import AMD_Tools3 as AMD
    lalodomain = [35.0, 36.0, 135.0, 136.0]
    Ppad, lat, lon = AMD.GetGeoData('landluse_H210100', lalodomain)
    """
    lld = LatLonDomain(*lalodomain,area)
    if area is None:
        area = lld.get_area()
    filename = 'AMD_' + area +'_Geo_' + element + '.nc.nc'
    opendap_source = urljoin([url,area,'GeoData',filename]) + '?' + element + lld.getIdx()
    dh,dhpath = url2dh(opendap_source)
    latitude = dh.variables['lat'][:]
    yr = lld.latrestrict(latitude)
    lat = latitude[yr]
    longitude = dh.variables['lon'][:]
    xr = lld.lonrestrict(longitude)
    lon = longitude[xr]
    ncGeo =  dh.variables[element]
    name = ncGeo.long_name
    unit = ncGeo.units
    Ge  = np.array(ncGeo[yr, xr])
    Geo = np.where(Ge==ncGeo._FillValue,np.nan,Ge)
    dh.close()
    print(element, Geo.shape)
    if dhpath is not None:
        StartUnlink(dhpath)

    if len(lat) and lat[0] > lat[-1]:
        lat = lat[::-1]
        Geo = Geo[::-1,:]

    if namuni:
        return Geo, lat, lon, name, unit
    else:
        return Geo, lat, lon

def GetCSV_Table(filename,delimiter=","):
    """
概要：
    CSVファイルのように、カンマ等で区切られた表を内容とするテキストファイルを読み込み、
    それを文字列のリストとして返す関数。数表先頭の１行はヘッダーとみなす。さらに、この
　行のフィールド数をこの数表のフィールド数とみなし、これと異なるフィールド数のレコード
　は戻さない。
書式：
    GetCSV_Table(filename,delimiter=",")
引数(必須)：
    filename：読み込むべきCSVファイルの名前。
引数(必要に応じ指定)：
    delimiter：区切り文字を引用符で括って指定する。”\t”としてすればタブ区切りの表を読むことができる。
            指定を省略した場合はカンマが区切り文字に設定され、CSVファイルを読むことができる。
戻り値：
    header：見出し行の内容(文字列のリスト)
    body：フィールドを要素に持つリストのリスト
備考：
    例として、表の第3フィールドの見出しはheader[2]に格納される。また、5番目のデータの
    第3フィールドの内容はbody[4][2]に文字列で格納される。

    """
    #fh = open(filename,'r',encoding="utf-8_sig")
    fh = getFileContent(filename)
    lines = fh.replace('"','').split('\n')  #二重引用符を除去して改行で区切りる。
    nore = len(lines)
    header = lines[0].split(delimiter)  #カンマで区切る。
    nofi = len(header)
    body = []
    for i in range(nore):
        record = lines[i].split(delimiter)
        if len(record) == nofi:
            body.append(record)
    nore = len(body)
    return header, body[1:nore]



def GetCSV_Map(filename, skiprow=0,upsidedown=True):
    """
概要：
    CSV形式のテキストファイルを読み込み、それを浮動小数のNumpy Array Objectとして返す関数。
    配列の列数は、取り込み範囲の先頭行で判別する。行数は、EOFまでの行数から判別する。文字列
    「nan」は、numpy.nanとして理解される。これ以外の文字列が検出されると警告文を表示し、
    nanを返す。
書式：
    GetCSV_Map(filename, skiprow=0,upsidedown=True)
引数(必須)：
    filename：読み込むべきCSVファイルの名前。
引数(必要に応じ指定)：
    skiprow：余白や見出しなどに使用されていて読み込み対象としない行の数。
            指定を省略した場合は0に設定される。
    upsidedown：余白や見出しなどに使用されていて読み込み対象としない行の数。
            指定を省略した場合は0に設定される。
戻り値：
    numpyの浮動小数点配列。数値以外には無効値(nan)が与えられる。
    """
    llen = {}
    records = []
    f = getFileContent(filename)
    lines = f.strip().split("\n")
    for i,x in enumerate(lines):
        if "," in x:
            ws = [y.strip() for y in x.split(",")]
        else:
            ws = x.split()
        line = []
        for y in ws:
            try:
                y = float(y)
            except:
                y = np.nan
            line.append(y)
        records.append(line)
        if i >= skiprow:
            n = len(line)
            llen[n] = llen.get(n,0)+1
            if np.nan in line:
                print("Warning: possibly illegal elements")
                print("Row " + str(i+1) + ":",x)
    if len(llen.keys()) >= 2: # check if all lines have the same number of items
        print("Error: variable number of items per line",llen)
    Var = np.array(records[skiprow:], dtype=np.float32)
    if upsidedown :
        Var = Var[::-1,:]
    print(filename+": ", Var.shape)
    return Var



def PutCSV_TS(Var, tim, header=None, filename='result.csv'):
    """
概要：
    時系列のデータをCSV形式のファイルで出力する関数
書式：
    PutCSV_TS(Var, tim, header=None, filename='result.csv')
引数(必須)：
     Var：時系列の1次元配列データ。ただし、「Ver=np.array([V1,V2,..,Vn])」とすれ
        ば、n個の1次元配列V1,V2,..,Vnを一度に出力することができる。
     tim：時刻の見出しとして使用される一次元配列。第1列に行方向に出力される。
引数(必要に応じ指定)：
     header：引数に「header='moji,retsu'」として文字列を与えると、出力ファイルの
         第1行目にこの文字列が出力される。
     filename：引数に「filename='fairun_no_namae.csv'」として文字列を与えると、
         ファイルがこの名前で出力される。これを指定しない場合は、デフォルトの
         ファイル名「result.csv」が用いられる。
戻り値：なし
    """
    dim = len(Var.shape)
    assert dim in (1,2)

    with open(filename, 'wt') as f:
        if header != None:
            f.write(header + '\n')
        if dim == 1:
            for i in range(0, Var.shape[0]):
                f.write(str(tim[i]) + "," + str(Var[i]) + "\n")
        else:
            for i in range(0, Var.shape[1]):
                f.write(str(tim[i]) + "," + ",".join([str(y) for y in Var[:,i]]) + "\n")




def PutCSV_Map(Var, lat, lon, filename='result.csv'):
    """
2次元の浮動小数点配列をCSV形式のファイルで出力する関数。
書式：
    PutCSV_Map(Var, lat, lon, filename='result.csv')
引数(必須)：
     Var：出力すべき2次元の配列
     lat：緯度の見出しとして使用される一次元配列。
     lon：経度の見出しとして使用される一次元配列。
引数(必要に応じ指定)：
     filename：引数に「filename='fairun_no_namae.csv'」として文字列を与えると、
         ファイルがこの名前で出力される。これを指定しない場合は、デフォルトの
         ファイル名「result.csv」が用いられる。
 戻り値：なし
    """
    with open(filename, 'wt') as f:
        exVar = np.zeros((Var.shape[0]+1, Var.shape[1]+1))
        exVar[0,1:] = lon
        exVar[1:,0] = lat[::-1]
        exVar[1:,1:] = Var[::-1,:]
        for i in range(0, exVar.shape[0]):
            f.write(",".join([str(x) for x in exVar[i,:]]) + "\n")


def PutCSV_MT(Var, lat, lon, addlalo=False, header=None, filename='result.csv', removenan=True, delimiter=','):
    """
概要：
    3次元の配列を、基準3次メッシュコードをキーとするテーブルをCSVファイルで出力する関数。
    メッシュ農業気象データのメッシュは、基準国土3次メッシュと一致しているので、あるメッシュに
    おける値をメッシュコードを行見出しとするテーブルにすることができる。３次メッシュポリゴン
    データ(三次メッシュコードを属性に持つ)を持つGISを用意し、このテーブルをインポートして三次
    メッシュコードをキーにして連結すると、テーブルの値をGISで表示することができる。
書式：
    PutCSV_MT(Var, lat, lon, addlalo=False, header=None, filename='result.csv',
        removenan=True, delimiter=',')
引数(必須)：
    Var:内容を書き出す配列変数。第0次元の内容を添え字の順に記号で区切って出力する。
    lat:配列Varの各行が位置する緯度値が格納されている配列。Varの第1次元の要素数と一致していなくてはならない。
    lon:配列Varの各列が位置する経度値が格納されている配列。Varの第2次元の要素数と一致していなくてはならない。
引数(必要に応じ指定)：
    addlalo:これをTrueにすると、3次メッシュ中心点の緯度と経度が第2フィールドと第3フィールド追加挿入される。デフォルトはFalseであり挿入されない。
    header:一行目に見出しやタイトルなど何か書き出すときはここに「header='文字列'」として指定する。
    filename：出力されるファイルの名前。デフォルト値は'result.csv'。
    removenan:無効値だけのレコードを削除するかを指定するキーワード。
        True=>無効値だけのレコードを削除する。水域を含む領域を出力するときに削除すると無駄なレコードが出ない。
        False=>Varに含まれるメッシュコードをすべて出力する。
    delimiter:フィールドの区切り文字。デフォルト値は','。すなわち、CSVファイルとなる。
 戻り値：なし。
    """
    if len(Var.shape) == 2:     #2次元配列の場合は3次元配列にする。
#        Var = np.ma.array(Var, ndmin=3)
        Var = np.array(Var, ndmin=3)
    Var = np.where(Var == 9.96921E+36, np.nan, Var)  #NCLにおけるmissing　value
    noti = Var.shape[0]
    nola = Var.shape[1]
    nolo = Var.shape[2]
    #配列要素数のチェック。
    if nola != len(lat) or nolo != len(lon):
        print('エラー：緯度/経度の情報が整合していないのでメッシュコードを生成できません。')
    fh = open(filename, 'wt')
    #ヘッダが指定されていたらそれを書き出す。
    if header != None:
        fh.write( header+'\n' )
    for y in range(nola):
        for x in range(nolo):
            if any([not np.isnan(v) for v in Var[:,y,x]]) or not removenan:
                line = [lalo2mesh(lat[y],lon[x])]
                if addlalo == True:
                    line += [str(lat[y]),str(lon[x])]
                line += [str(Var[t,y,x]) for t in range(noti)]
                fh.write(delimiter.join(line) + '\n')
    fh.close()




def PutNC_Map( Var, lat, lon, description='Variable', symbol='Var', unit='--', fill=9.96921e+36, filename='result.nc'):
    """
概要：
    2次元(空間分布)の気象変量をnetCDF形式のファイルで出力する関数。
書式：
     PutNC_Map( Var, lat, lon, description='Variable', symbol='Var', unit='--', fill=9.96921e+36, filename='result.nc')
引数(必須)：
    Var:内容を書き出す配列変数。
    lat:配列Varの各行が位置する緯度値が格納されている配列。Varの第1次元の要素数と一致していなくてはならない。
    lon:配列Varの各列が位置する経度値が格納されている配列。Varの第2次元の要素数と一致していなくてはならない。
引数(必要に応じ指定)：
    description:データの正式名称などデータを説明する文。
    symbol:データに対して用いる記号
    unit:データの数値が従う単位
    fill：無効値として取り扱う数値。デフォルト値は、9.96921e+36
    filename：出力されるファイルの名前。デフォルト値は'result.nc'。
 戻り値：なし。
    """
    #----------------------------------------------------------------------------------
    fh = Dataset(filename,"w")
    fh.title = 'Map data crated by PutNC_Map function.'
    fh.source = "Agro-Meteorological Grid Square Data System NIAES, NARO"
    fh.Conventions = "None"
    fh.creation_date = "Created " + dt.now().strftime("%Y/%m/%d %H:%M:%S JST")
    fh.createDimension("lat",  len(lat))
    fh.createDimension("lon",  len(lon))
    nc_data = fh.createVariable(symbol, "f4", ("lat", "lon",), fill_value=fill)
    nc_lat  = fh.createVariable("lat" , "f4", ("lat",))
    nc_lon  = fh.createVariable("lon" , "f4", ("lon",))
    nc_lat.long_name = "latitude"
    nc_lat.units = "degrees_north"
    nc_lon.long_name = "longitude"
    nc_lon.units = "degrees_east"
    nc_data.long_name =  description + '(' + symbol + ')'
    nc_data.units = unit
    nc_data[:,:] = Var
    nc_lat[:]  = lat
    nc_lon[:]  = lon
    fh.close()


def PutNC_3D( Var, tim, lat, lon, description='None', symbol='Var', unit='--', fill=9.96921e+36, filename='result.nc' ):
    """
概要：
    3次元(空間分布×時間変化)の気象変量をnetCDF形式のファイルで出力する関数。
書式：
     PutNC_3D( Var, tim, lat, lon, description='None', symbol='Var', unit='--', fill=9.96921e+36, filename='result.nc' )
引数(必須)：
    Var:内容を書き出す3次元配列変数。第0次元は日付、第1次元は緯度、第2次元は経度とする。
    tim：配列Varがカバーする日付の配列。日付はPythonの日付オブジェクトで与える。Varの第0次元の要素数と一致していなくてはならない。
    lat:配列Varの各行が位置する緯度値が格納されている配列。Varの第1次元の要素数と一致していなくてはならない。
    lon:配列Varの各列が位置する経度値が格納されている配列。Varの第2次元の要素数と一致していなくてはならない。
引数(必要に応じ指定)：
    description:データの正式名称などデータを説明する文。
    symbol:データに対して用いる記号
    unit:データの数値が従う単位
    fill：無効値として取り扱う数値。デフォルト値は、9.96921e+36
    filename：出力されるファイルの名前。デフォルト値は'result.nc'。
 戻り値：なし。
    """
    fh = Dataset(join(filename),"w")
    fh.title = '3-D data crated by PutNC_3D function.'
    fh.source = "Agro-Meteorological Grid Square Data System NIAES, NARO"
    fh.Conventions = "None"
    fh.creation_date = "Created " + dt.now().strftime("%Y/%m/%d %H:%M:%S JST")
    fh.createDimension("time", len(tim))
    fh.createDimension("lat",  len(lat))
    fh.createDimension("lon",  len(lon))
    nc_Var = fh.createVariable(symbol, "f4", ("time", "lat", "lon",), fill_value=fill)
    nc_time = fh.createVariable("time", "f8", ("time",), fill_value=9.969209968386869E36)
    nc_lat  = fh.createVariable("lat" , "f4", ("lat",))
    nc_lon  = fh.createVariable("lon" , "f4", ("lon",))
    nc_time.calendar   = "standard"
    nc_time.units      = "days since 1900-1-1 00:00:0.0"
    nc_lat.long_name = "latitude"
    nc_lat.units = "degrees_north"
    nc_lon.long_name = "longitude"
    nc_lon.units = "degrees_east"
    nc_Var.long_name = description + '(' + symbol + ')'
    nc_Var.units = unit
    nc_Var[:,:,:] = Var
    nc_time[:] = date2num(tim[:], units=nc_time.units )
    nc_lat[:]  = lat
    nc_lon[:]  = lon
    fh.close()


#以下は make_kmのための関数
def fig_ax(lon0, lat0, lon1, lat1, pixels=1024, asp=None):
    "matplotlib `fig` and `ax` handles"
    if asp:
        aspect = asp
    else:
        aspect = np.cos(np.mean([lat0, lat1]) * np.pi/180.0)
    xsize = np.ptp([lon1, lon0]) * aspect
    ysize = np.ptp([lat1, lat0])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        plt.ioff()
    fig = plt.figure(figsize=figsize,frameon=False,dpi=pixels//10)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(lon0, lon1)
    ax.set_ylim(lat0, lat1)
    return fig, ax

def make_html(lon0, lat0, lon1, lat1, figs, colorbar, htmlfile, name):
    txt = """
<!DOCTYPE html>
<html style="height: 100%; width: 100%;">
<head>
<meta charset="UTF-8">
<title>"""+name+"""</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.2.0/dist/leaflet.css" />
<script src="https://unpkg.com/leaflet@1.2.0/dist/leaflet.js"></script>
<style type="text/css">
<!--
input[type=range]::-ms-tooltip {
    display: none;
}
-->
</style>
</head>
<body style="padding: 0; margin: 0; height: 100%; width: 100%;">
<div id="map" style="height: 100%; width: 100%;">
<script>
var lat0 = """+str(lat0)+""";
var lat1 = """+str(lat1)+""";
var lon0 = """+str(lon0)+""";
var lon1 = """+str(lon1)+""";
var maxrange = Math.max([lon1-lon0, lat1-lat0]);
var bmap_opacity = 1;
var omap_opacity = 0.6;
var maxrange = Math.max(lon1-lon0, lat1-lat0);
var z = 5;
for(var i=1; i<5; i++) if(maxrange<5*i) z+=1;
var mapname = \""""+name+"""\";
var overlay = \""""+figs[0]+"""\";
var legend  = \""""+colorbar+"""\";
var pre_map = "国土地理院 ";
var bslider = "<br><input id='bslider' type='range' min='0' max='1' step='0.01' value='1' style='width:100%;'></input>";
var oslider = "<br><input id='oslider' type='range' min='0' max='1' step='0.01' value='0.6' style='width:100%;'></input>";
var gsi_attr = "<a href='http://www.gsi.go.jp/kikakuchousei/kikakuchousei40182.html' target='_blank'>GSI</a>";
var basemap = {"標準地図":L.tileLayer("http://cyberjapandata.gsi.go.jp/xyz/std/{z}/{x}/{y}.png",{attribution:gsi_attr, opacity:bmap_opacity}),
               "淡色地図":L.tileLayer("http://cyberjapandata.gsi.go.jp/xyz/pale/{z}/{x}/{y}.png",{attribution:gsi_attr, opacity:bmap_opacity}),
               "白地図":L.tileLayer("http://cyberjapandata.gsi.go.jp/xyz/blank/{z}/{x}/{y}.png",{attribution:gsi_attr, opacity:bmap_opacity}),
               "写真":L.tileLayer("http://cyberjapandata.gsi.go.jp/xyz/seamlessphoto/{z}/{x}/{y}.jpg",{attribution:gsi_attr, opacity:bmap_opacity}),
               "標高段彩":L.tileLayer("http://cyberjapandata.gsi.go.jp/xyz/relief/{z}/{x}/{y}.png",{attribution:gsi_attr, opacity:bmap_opacity}),
               "英語版":L.tileLayer("http://cyberjapandata.gsi.go.jp/xyz/english/{z}/{x}/{y}.png",{attribution:gsi_attr, opacity:bmap_opacity})};
var map = L.map("map", {center: [(lat0+lat1)/2, (lon0+lon1)/2], zoom: z, maxZoom: 18, layers: basemap["白地図"]});
var temp = L.imageOverlay(overlay, [[lat0, lon0], [lat1, lon1]], {opacity:omap_opacity});
temp.addTo(map);
L.Control.Watermark = L.Control.extend({
    onAdd: function(map) {
        var img = L.DomUtil.create("img");
        img.src = legend;
        img.style.width = "100px";
        return img;
    }
});
L.control.watermark = function(opts) {return new L.Control.Watermark(opts);}
L.control.watermark({ position: "bottomleft" }).addTo(map);
var baseLayers = {};
for (var i=0; i<Object.keys(basemap).length; i++) {
    if (i+1 == Object.keys(basemap).length) baseLayers[pre_map+Object.keys(basemap)[i]] = basemap[Object.keys(basemap)[i]];
    else baseLayers[pre_map+Object.keys(basemap)[i]] = basemap[Object.keys(basemap)[i]];
}
var overlays = {};
overlays[mapname] = temp;
L.control.layers(baseLayers, overlays).addTo(map);

var blabel_element  = document.getElementsByClassName("leaflet-control-layers-base")[0];
var blabel = document.createElement("label");
blabel.innerHTML = bslider;
blabel_element.appendChild(blabel);

var olabel_element  = document.getElementsByClassName("leaflet-control-layers-overlays")[0];
var olabel = document.createElement("label");
olabel.innerHTML = oslider;
olabel_element.appendChild(olabel);

document.getElementById('oslider').addEventListener('input',  function() {temp.setOpacity(+document.getElementById('oslider').value)});
document.getElementById('oslider').addEventListener('change',  function() {temp.setOpacity(+document.getElementById('oslider').value)});
document.getElementById('bslider').addEventListener('input',  function() {
    for (var i=0; i<Object.keys(basemap).length; i++) {
         basemap[Object.keys(basemap)[i]].setOpacity(+document.getElementById('bslider').value);
    }
});
document.getElementById('bslider').addEventListener('change',  function() {
    for (var i=0; i<Object.keys(basemap).length; i++) {
         basemap[Object.keys(basemap)[i]].setOpacity(+document.getElementById('bslider').value);
    }
});
</script>
</body>
</html>
"""
    with codecs.open(htmlfile, "w", "utf8") as f:
        f.write(txt)

def map_figs(data, lat, lon, label, cmapstr, minmax, filename, overlay, legend, asp=None):
    lat, lon = np.meshgrid(lat, lon)
    pixels = 1024 * 10
    # data = data.transpose() # pcolormesh needs nans converted to mask
    if data.dtype == np.dtype('<M8[D]'):
        data = ma.array(data,mask=data==np.datetime64("1-01-01", "D")).transpose() # pcolormesh needs nans converted to mask
        fig, ax = fig_ax(lon0=lon.min(),lat0=lat.min(),lon1=lon.max(),lat1=lat.max(),pixels=pixels,asp=asp)
        if label is None:
            label = filename
        if cmapstr is None:
            # cmap = colorbrewer.get_map('RdYlGn', 'diverging', 11, reverse=True).mpl_colormap
            #cmap = cm.RdYlGn_r
            cmap = copy.copy(cm.get_cmap("RdYlGn_r"))
        else:
            #cmap = eval("cm."+cmapstr)
            cmap = copy.copy(cm.get_cmap(cmapstr))
        if minmax is None:
            sclint = 1  #何日ごとに色分けするか
            sclmin = data.min()     #何月何日から色を付けるか
            sclmax = data.max()     #何月何日まで色を付けるか
            levels = np.arange(sclmin, sclmax+np.timedelta64(sclint,'D')+sclint, sclint)
        else:
            sclint = 1  #何日ごとに色分けするか
            sclmin = minmax[0]     #何月何日から色を付けるか
            sclmax = minmax[1]     #何月何日まで色を付けるか
            levels = np.arange(sclmin, sclmax+np.timedelta64(sclint,'D')+sclint, sclint)
        cs = ax.contourf(lon, lat, data, levels, cmap=cmap)
    else:
        data = ma.array(data,mask=np.isnan(data)).transpose() # pcolormesh needs nans converted to mask
        fig, ax = fig_ax(lon0=lon.min(),lat0=lat.min(),lon1=lon.max(),lat1=lat.max(),pixels=pixels,asp=asp)
        if label is None:
            label = filename
        if cmapstr is None:
            # cmap = colorbrewer.get_map('RdYlGn', 'diverging', 11, reverse=True).mpl_colormap
            #cmap = cm.RdYlGn_r
            cmap = copy.copy(cm.get_cmap("RdYlGn_r"))
        else:
            #cmap = eval("cm."+cmapstr)
            cmap = copy.copy(cm.get_cmap(cmapstr))
        if minmax is None:
            cs = ax.pcolormesh(lon, lat, data, cmap=cmap, shading="auto")
        else:
            cs = ax.pcolormesh(lon, lat, data, vmin=minmax[0], vmax=minmax[1], cmap=cmap, shading="auto")
    ax.set_axis_off()
    fig.savefig(overlay, transparent=True, format='png', dpi=100)
    plt.close()

    fig = plt.figure(figsize=(1.3, 4.0), facecolor=None, frameon=False)
    ax = fig.add_axes([0.75, 0.05, 0.2, 0.9])
    if data.dtype == np.dtype('<M8[D]'):
        cmap.set_over('w', 1.0)   #上限を超えたときは白色
        cmap.set_under('k', 1.0)  #下限を超えたときは黒色
        cmap.set_bad('w', 1.0)
        cb = fig.colorbar(cs, cax=ax, format=DateFormatter('%b %d'))
    else:
        cb = fig.colorbar(cs, cax=ax)
    cb.ax.yaxis.set_label_position('left')
    cb.ax.yaxis.set_ticks_position('left')
    cb.ax.tick_params(labelsize=12)
    cb.set_label(label, rotation=90, color='k', labelpad=5, fontsize=12)
    fig.savefig(legend, transparent=False, bbox_inches='tight', pad_inches=0, format='png', dpi=100)
    plt.close()

def PutGSI_Map(data, lat, lon, label=None, cmapstr=None, minmax=None, filename="result", outdir="."):
    """
概要：
    2次元(空間分布)の配列を地理院地図オーバーレイ用のHTMLファイルで出力する関数。
書式：
    PutGSI_Map(data, lat, lon, label=None, cmapstr=None, minmax=None, filename="result", outdir="."):
引数(必須)：
    data：表示させるデータ（2次元numpyアレイ）
    lat：緯度（1次元numpyアレイ）
    lon：経度（1次元numpyアレイ）
引数(必要に応じ指定)：
    label：凡例のタイトルの文字列
    cmapstr：カラーマップを指定（詳細は後述）
    minmax：カラースケールを指定（[min,max]のリスト型）
    filename：出力ファイル名
    outdir:出力フォルダ名
カラーマップについて：
    カラーマップには名称があるのでこれを文字列で("で囲んで)指定する。
    例)　レインボーカラー:rainbow、黄色-オレンジ-赤の順で変化:YlOrRdなど
    色の順序をを反転させたい場合は、rainbow_rのよう名称の後ろに"_r"を付加する。
    詳細は下記URLを参照。
        http://matplotlib.org/examples/color/colormaps_reference.html
注意：
    この関数は、numpy.datetime64[D]型(日時のを格納)の配列も表示することができるので、日付の
    分布図を描画できます。但し、等値色の塗り方はメッシュ単位ではなく、メッシュ中心点の値を等高線
    で結ぶような描画になります。
　　　また、オプション引数minmaxで配色の下限と上限を指定する場合は、下記のようにして日付値を指定
    してください(datetimeオブジェクトで指定すると誤作動します)。
    minmax = [np.datetime64('2018-08-25','D'),np.datetime64('2018-09-05','D')]

使用例：
北緯36～38.5度、東経137.5～141.5度の範囲における2016年1月1日の平均気温分布図のHTMLファイルを作成する。
import AMD_Tools3 as AMD
element = 'TMP_mea'
timedomain = [ "2016-01-01", "2016-01-01" ]
lalodomain = [ 36.0, 38.5, 137.5, 141.5]
Msh,tim,lat,lon,nam,uni = AMD.GetMetData(element, timedomain, lalodomain,namuni=True)
dat = Msh[0,:,:]
AMD.PutGSI_Map(dat,lat,lon,label=nam+" ["+uni+"]", cmapstr="rainbow",minmax=None, filename=element)
    """
    if not exists(outdir):
        print("DirectoryFolder",outdir,"does not exists")
        return
    if not isdir(outdir):
        print("Path",outdir,"is not a directory")
        return
    from pyproj import Transformer
    # P3857 = Proj(init="epsg:3857")
    # P4326 = Proj(init="epsg:4326")
    transformer = Transformer.from_crs(4326, 3857, always_xy=True)
    pts = [(x,y) for x,y in zip(lon,lat)]
    if len(lon)>len(lat):
        pts.extend([(x,lat.mean()) for x in lon[len(lat):]])
    elif len(lat)>len(lon):
        pts.extend([(lon.mean(),y) for y in lat[len(lon):]])
    trans_pts = [pt for pt in transformer.itransform(pts)]
    lon2 = np.array([pt[0] for pt in trans_pts][:len(lon)])
    lat2 = np.array([pt[1] for pt in trans_pts][:len(lat)])
    # lon2 = np.array([transform(P4326, P3857, x, lat.mean())[0] for x in lon])
    # lat2 = np.array([transform(P4326, P3857, lon.mean(), y)[1] for y in lat])
    hdlat = 0.5 if len(lat) < 2 else (lat[1] - lat[0]) * 0.5
    hdlon = 0.5 if len(lon) < 2 else (lon[1] - lon[0]) * 0.5
    overlay = join(outdir,filename+"_o.png")
    legend  = join(outdir,filename+"_l.png")
    html    = join(outdir,filename+".html")
    map_figs(data, lat2, lon2, label, cmapstr, minmax, filename, overlay, legend, asp=1)
    make_html(lon0=float(lon.min()-hdlat), lat0=float(lat.min()-hdlon),lon1=float(lon.max()), lat1=float(lat.max()),
              figs=[basename(overlay)], colorbar=basename(legend), htmlfile=html, name=filename)



def GetGeoData_Tile(element, lalodomain, namuni=False,
               url='https://amd.rd.naro.go.jp/opendap/AMD/'):
    lld = LatLonDomain(*lalodomain)
    dh,dhpath = {},{}
    for code, cIdx in lld.getCodeWithIdx():
        filename = f'AMDy____p{code}g{element}.nc.nc'
        opendap_source = urljoin([url,'geodata',f'g{element}',filename]) + '?' + element + cIdx
        dh[code],dhpath[code] = xurl2dh(opendap_source)

    if len(dh.keys()) == 0:
        print("No data to retrieve. Please check lat-lon or time domain.")
        if namuni:
            return None, None, None, None, None, None
        else:
            return None, None, None, None
    m = xr.merge(dh.values())

    latitude = m.variables['lat'][:]
    yrest = lld.latrestrict(latitude)
    lat = latitude[yrest]
    longitude = m.variables['lon'][:]
    xrest = lld.lonrestrict(longitude)
    lon = longitude[xrest]
    ncGeo =  m.variables[element]
    name = ncGeo.attrs["long_name"]
    unit = ncGeo.attrs["units"]
    Geo = np.array(ncGeo[yrest, xrest])

    for ds in dh.values():
        ds.close()
    #print(element, Geo.shape)

    for path in dhpath.values():
        if path is not None:
            StartUnlink(path)

    if len(lat) and lat[0] > lat[-1]:
        lat = lat[::-1]
        Geo = Geo[::-1,:]

    if namuni:
        return Geo, lat, lon, name, unit
    else:
        return Geo, lat, lon



def GetMetData_Tile(element, timedomain, lalodomain,
                 cli=False, namuni=False, url='https://amd.rd.naro.go.jp/opendap/AMD/'):
    lld = LatLonDomain(*lalodomain)
    td = TimeDomain(*timedomain)
    dh,dhpath = {},{}
    ec = "c" if cli else "e"
    for year, tidx in td.getIdx():
        for code, cidx in lld.getCodeWithIdx():
            filename = f'AMDy{year}p{code}{ec}{element}.nc.nc'
            opendap_source = urljoin([url,f'{year}',f'{ec}{element}',filename]) + '?' + element + tidx + cidx
            #print(opendap_source)
            dh[code,year],dhpath[code,year] = xurl2dh(opendap_source)
    #print(dh.keys())
    if len(dh.keys()) == 0:
        print("No data to retrieve. Please check lat-lon or time domain.")
        if namuni:
            return None, None, None, None, None, None
        else:
            return None, None, None, None
    m = xr.merge(dh.values())
    tim,lat,lon,Met,name,unit = xtll_extract(m,td,lld,element)
    for path in dhpath.values():
        if path is not None:
            StartUnlink(path)

    if namuni:
        return Met, tim, lat, lon, name, unit
    else:
        return Met, tim, lat, lon

def GetSceData_Tile(element, timedomain, lalodomain, model, scenam,
               cli=False, namuni=False, url='https://amd.rd.naro.go.jp/opendap/AMS'):
    lld = LatLonDomain(*lalodomain)
    td = TimeDomain(*timedomain)
    dh,dhpath = {},{}
    for year, tidx in td.getIdx():
        for code, cidx in lld.getCodeWithIdx():
            filename = f'AMSy{year}p{code}e{element}.nc.nc'
            opendap_source = urljoin([url,model,scenam,f'{year}',f'e{element}',filename]) + '?' + element + tidx + cidx
            #print(opendap_source)
            dh[code,year],dhpath[code,year] = xurl2dh(opendap_source)
    #print(dh.keys())
    if len(dh.keys()) == 0:
        print("No data to retrieve. Please check lat-lon or time domain.")
        if namuni:
            return None, None, None, None, None, None
        else:
            return None, None, None, None
    m = xr.merge(dh.values())
    tim,lat,lon,Met,name,unit = xtll_extract(m,td,lld,element)
    for path in dhpath.values():
        if path is not None:
            StartUnlink(path)

    if namuni:
        return Met, tim, lat, lon, name, unit
    else:
        return Met, tim, lat, lon


def GetMetData(*args, **kargs):
    if 'isTile' in kargs:
        if kargs['isTile']:
            kargs.pop('isTile')
            return GetMetData_Tile(*args, **kargs)
        else:
            kargs.pop('isTile')
    return GetMetData_Area(*args, **kargs)


def GetGeoData(*args, **kargs):
    if 'isTile' in kargs:
        if kargs['isTile']:
            kargs.pop('isTile')
            return GetGeoData_Tile(*args, **kargs)
        else:
            kargs.pop('isTile')
    return GetGeoData_Area(*args, **kargs)


def GetSceData(*args, **kargs):
    if 'isTile' in kargs:
        if kargs['isTile']:
            kargs.pop('isTile')
            return GetSceData_Tile(*args, **kargs)
        else:
            kargs.pop('isTile')
    return GetSceData_Area(*args, **kargs)


def GetMetDataHourly(element, timedomain, lalodomain,
                 cli=False, namuni=False, url='https://amd.rd.naro.go.jp/opendap/AMD_Hourly/'):
    lld = LatLonDomain(*lalodomain)
    td = TimeDomainHourly(*timedomain)
    dh,dhpath = {},{}
    lgg = lld.geogrid()
    tgg = td.geogrid()
    ec = "c" if cli else "e"

    for year, tidx in td.getIdx():
        for code, cidx in lld.getCodeWithIdx():
        #print("CY",code,year)
            filename = f'AMDy{year}p{code}{ec}_h_{element}.nc.nc'
            opendap_source = urljoin([url,f'{year}',f'{ec}{element}',filename]) + '?' + element + tidx + cidx
            #print(opendap_source)
            dh[code,year],dhpath[code,year] = xurl2dh(opendap_source)
    #print(dh.keys())
    if len(dh.keys()) == 0:
        print("No data to retrieve. Please check lat-lon or time domain.")
        if namuni:
            return None, None, None, None, None, None
        else:
            return None, None, None, None
    m = xr.merge(dh.values())
    tim,lat,lon,Met,name,unit = xtll_extract(m,td,lld,element)
    for path in dhpath.values():
        if path is not None:
            StartUnlink(path)

    if namuni:
        return np.round(Met,decimals=1), tim, lat, lon, name, unit
    else:
        return np.round(Met,decimals=1), tim, lat, lon


if __name__ == "__main__":
    print("AMD_Tools3.pyは、実行したいプログラムのあるフォルダに配置してください。")
    print("これ自身を実行しても何も起こりません。配置だけしてください。")
