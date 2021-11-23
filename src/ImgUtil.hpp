
/*****************************************************************************
*                                                                            *
*                            ImgUtil          								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : ImgUtilr.hpp
    Author(S)       : Me Eunkyung
    Created         : 23 Oct 2021

    Description     : ImgUtil.Hpp
    Notes           : Utility related image.
*/
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>
#include "DefData.hpp"
#include "MtrxUtil.hpp"

using namespace std;
using namespace cv;

class ImgUtil {

    public :
        MtrxUtil mtrx;        
        string fname;

        ImgUtil();

        void SaveImage(SCENE* sc, int type = 0, SCENE* sc2 = 0, PARAM* p = 0, int opt = -1);
        vector<Mat>LoadImages(const string& path, vector<string>* dsc_id);     
        void SimpleColorBalance(Mat& in, Mat& out, float percent);
        int AdjustImage(SCENE* sc, ADJST adj);
        void ColorCorrection(Mat& ref, Mat& src, Mat& out);
        void SetLogFilename(std::string fn) { fname = fn; };

};