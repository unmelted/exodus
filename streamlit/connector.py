import os
import pathlib
import genesis as gn

class BaseData(object) : 
    instance = None
    cur_path = os.getcwd() +'/'
    os.chdir('../')
    prj_path = os.getcwd() +'/'
    img_path = prj_path + 'image/new/'
    libname = prj_path + 'libexodus.dylib'
    print(cur_path, prj_path, img_path)

    @staticmethod
    def getInstance():
        if BaseData.instance == None:
            BaseData.instance = BaseData()
        return BaseData.instance 

    def getImageList(self):
        newlist = []
        imglist = os.listdir(self.img_path)
        imglist.sort()
        for i in imglist :
            if i == '.DS_Store':
                continue
            newlist.append(self.img_path + i)

        return newlist

    def getGroundType(self):
        gr_type = { "BaseballHome" : 1,
                    "BaseballGround" : 2,
                    "BasketballHalf" : 3,
                    "BasketballGround" : 4,
                    "Boxing" : 5,
                    "IceLinkHalf" : 6,
                    "IceLink" : 7,
                    "SoccerHalf" : 8,
                    "Soccer" : 9,
                    "Taekwondo" : 10,
                    "TennisHalf" : 11 ,
                    "Tennis" : 12,
                    "Ufc" : 13,
                    "VolleyballHalf" : 14,
                    "VolleyballGround" : 15,
                    "Football" : 16 }
        return gr_type


   
class Handler(object):
    instance = None
    data = None
    instance = None
    ground = None
    imageset = None
    region = []
    dim = None
    bd = BaseData.getInstance()


    @staticmethod
    def getInstance():
        if Handler.instance == None :
            Handler.instance = Handler()
        return Handler.instance

    def __init__(self) :
        gn.Calibrator.getInstance().setLib(self.bd.libname)

    def ExecuteExtract(self) :
        temp = []
        temp.append(self.dim)
        print("Execute Extract", len(self.region))
        for i in self.region :
            temp.append(i)
        gn.Calibrator.getInstance().extract(self.dim, temp, self.bd.img_path)

    def setGround(self, ground):
        self.ground = ground
        print("setGround is called ", self.ground)

    def setRegion(self, rg):
        for i in rg:
            self.region.append(int(i))

        self.dim = int(len(self.region) / 3)
        print(self.dim, self.region)

    def setImgData(self, imageset):
        self.imageset = imageset




