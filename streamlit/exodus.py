
'''
*****************************************************************************
*                                                                            *
*                            Exodus          								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : exodus.py
    Author(S)       : Me Eunkyung
    Created         : 23 Nov 2021

    Description     : exodus.py
    Notes           : Python - C  connector in python side.
'''

import os
import cv2
import ctypes as c
import pathlib
import numpy as np


""" imglist = os.listdir(datapath)
print(imglist)

for i in imglist :
    if i == ".DS_Store" :
        continue
        
    filename = datapath + i
    imgsrc = cv2.imread(filename)
    testarray1 = np.frombuffer(imgsrc, np.uint8)
    testarray2 = np.reshape(testarray1, (1024, 688, 3))
    framearray = testarray2.tobytes()
    print("CAll clib")
    clib.Feature( framearray)
 """

class Calibrator(object) :
    instance = None
    clib = None

    @staticmethod
    def getInstance():
        if Calibrator.instance == None:
            Calibrator.instance = Calibrator()
        return Calibrator.instance

    def setLib(self, libname):
        self.clib = c.CDLL(libname)
        print(self.clib)


    def extract(self, dim, gr_line, im_line):
        print("exodus.py extract is called")
        d = (c.c_int * ( (dim* 4) + 1))(*gr_line)
        e = (c.c_int * ( (dim* 4) + 1))(*im_line)
        #ch = c.create_string_buffer(img_path)
        print("Calibrator is called ")        
        self.clib.Extract(c.byref(d), c.byref(e))

