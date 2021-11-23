import os
import pathlib
import exodus as gn

class BaseData(object) : 
    instance = None
    cur_path = os.getcwd() +'/'
    os.chdir('../')
    prj_path = os.getcwd() +'/'
    img_path = prj_path + 'image/new/'
    libname = prj_path + 'libexodus.dylib'
    gr_img = prj_path + "football.png"

    print(libname)

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
    gr_line = []
    img_line = []    
    dim = None
    bd = BaseData.getInstance()


    @staticmethod
    def getInstance():
        if Handler.instance == None :
            Handler.instance = Handler()
        return Handler.instance

    def ExecuteExtract(self) :
        gn.Calibrator.getInstance().setLib(self.bd.libname)        
        temp = []
        temp.append(self.dim)
        print("Execute Extract", self.dim)
        for i in self.gr_line :
            temp.append(int(i[0]))
            temp.append(int(i[1]))
            temp.append(int(i[2]))           
            temp.append(int(i[3]))
        for i in self.img_line :
            temp.append(int(i[0]))
            temp.append(int(i[1]))
            temp.append(int(i[2]))           
            temp.append(int(i[3]))

        print(temp)
        gn.Calibrator.getInstance().extract(self.dim, temp)

    def setGround(self, ground):
        self.ground = ground
        print("setGround is called ", self.ground)

    def setRegion(self, gr_line, img_line):
        print("Set Region is called ")
        self.gr_line.clear()
        self.img_line.clear()

        for i in gr_line:
            self.gr_line.append(i)
        for j in img_line:
            self.img_line.append(j)

        print(self.gr_line)
        print(self.img_line)

        self.dim = len(self.gr_line)
        print(self.dim)

    def setImgData(self, imageset):
        self.imageset = imageset




