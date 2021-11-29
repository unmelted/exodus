
'''
*****************************************************************************
*                                                                            *
*                            connector         								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : connector.py
    Author(S)       : Me Eunkyung
    Created         : 23 Nov 2021

    Description     : connector.py
    Notes           : data handler fomr ui.
'''

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
        temp1 = []
        temp1.append(self.dim)
        temp2 = []
        temp2.append(self.dim)

        for i in self.gr_line :
            temp1.append(int(i[0]))
            temp1.append(int(i[1]))
            temp1.append(int(i[2]))           
            temp1.append(int(i[3]))
        for i in self.img_line :
            temp2.append(int(i[0]))
            temp2.append(int(i[1]))
            temp2.append(int(i[2]))           
            temp2.append(int(i[3]))

        print(temp1)
        print(temp2)
        gn.Calibrator.getInstance().extract(self.dim, temp1, temp2)


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

    def setImgData(self, imageset):
        self.imageset = imageset




