  
/*****************************************************************************
*                                                                            *
*                            MtrxUtil         								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : MtrxUtil.hpp
    Author(S)       : Me Eunkyung
    Created         : 04 Oct 2021

    Description     : MtrxUtil.hpp
    Notes           : Matrix calculate uility
*/
#pragma once
#include <cmath>
#include "DefData.hpp"


using namespace std;
using namespace cv;


class MtrxUtil {

public:

    Point2f GetRotatePoint(Point2f ptCenter, Point2f ptRot, double dbAngle);
    FPt GetRotatePoint(FPt center, FPt target, double angle);
    Mat GetRotationMatrix(float rad, float cx, float cy);
    Mat GetScaleMatrix(float scalex, float scaley, float cx, float cy);
    Mat GetScaleMatrix(float scalex, float scaley);
    Mat GetTranslationMatrix(float tx, float ty);
    Mat GetRotationMatrix(float rad);

    int TransformPtsbyHomography(Pt* in, Mat& homography, int cnt);
    Pt TransformPtbyHomography(Pt* in, Mat& homography);
    FPt TransformPtbyHomography(FPt in, Mat& homography);
    Point2f TransformPtbyHomography(Point2f* in, Mat& homography);

    FPt TransformPtbyAffine(FPt in, Mat& aff);
    int Hamming(uchar* arr1, uchar* arr2, int size);
    int SignValue(double x);
};