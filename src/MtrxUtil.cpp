  
/*****************************************************************************
*                                                                            *
*                            MtrxUtil         								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : MtrxUtil.cpp
    Author(S)       : Me Eunkyung
    Created         : 04 Oct 2021

    Description     : MtrxUtil.Cpp
    Notes           : Matrix calculate uility
*/
#include "MtrxUtil.hpp"

#if defined _WIN_ || _WINDOWS
#include <nmmintrin.h>
#endif

Point2f MtrxUtil::GetRotatePoint(Point2f ptCenter, Point2f ptRot, double dbAngle)
{
    Point2f ptResult =Point2f();

    // 회전 중심좌표와의 상대좌표
    ptRot.x = ptRot.x - ptCenter.x;
    ptRot.y = -(ptRot.y - ptCenter.y);

    double cosx = cos(dbAngle);
    double sinx = sin(dbAngle);

    ptResult.x = (float)(ptRot.x * cosx - ptRot.y * sinx);
    ptResult.y = (float)(ptRot.x * sinx + ptRot.y * cosx);

    ptResult.x = ptResult.x + ptCenter.x;
    ptResult.y = -(ptResult.y - ptCenter.y);

    return ptResult;

}

FPt MtrxUtil::GetRotatePoint(FPt center, FPt target, double angle) {

    FPt rpt;
    float dx = target.x - center.x;
    float dy = -(target.y - center.y);

    double cosx = cos(angle);
    double sinx = sin(angle);
    rpt.x = (float)(dx * cosx - dy * sinx);
    rpt.y = (float)(dx * sinx - dy * cosx);

    rpt.x = rpt.x + center.x;
    rpt.y = -(rpt.y - center.y);

    return rpt;
}


Mat MtrxUtil::GetRotationMatrix(float rad, float cx, float cy) {

    Mat m = Mat::eye(3, 3, CV_32FC1);
    if( rad == 0)
        return m;
    
    Mat mtrana, mtranb, mrot;
    mtrana = GetTranslationMatrix(-cx, -cy);
    mrot = GetRotationMatrix(rad);
    mtranb = GetTranslationMatrix(cx, cy);

    m = mtranb * mrot * mtrana;
    return m;

}

Mat MtrxUtil::GetScaleMatrix(float scalex, float scaley, float cx, float cy) {

    Mat m = Mat::eye(3, 3, CV_32FC1);

    if( scalex == 0 && scaley == 0)
        return m;

    Mat mtrana, mtranb, mscale;

    mtrana = GetTranslationMatrix(-cx, -cy);
    mscale = GetScaleMatrix(scalex, scaley);
    mtranb = GetTranslationMatrix(cx, cy);

    m = mtranb * mscale * mtrana;
    return m;
}

Mat MtrxUtil::GetScaleMatrix(float scalex, float scaley) {

    Mat m = Mat::eye(3, 3, CV_32FC1);
    if (scalex == 0 && scaley == 0)
        return m;
    
    m.at<float>(0, 0) = scalex;
    m.at<float>(1, 1) = scaley;
    return m;

}

Mat MtrxUtil::GetTranslationMatrix(float tx, float ty) {

    Mat m = Mat::eye(3, 3, CV_32FC1);
    if (tx == 0 && ty == 0)
        return m;
    
    m.at<float>(0, 2) = tx;
    m.at<float>(1, 2) = ty;
    return m;
}

Mat MtrxUtil::GetRotationMatrix(float rad) {

    Mat m = Mat::eye(3, 3, CV_32FC1);
    if (rad == 0)
        return m;
    
    m.at<float>(0, 0) = (float)cos(rad);
    m.at<float>(0, 1) = (float)-sin(rad);
    m.at<float>(1, 0) = (float)sin(rad);
    m.at<float>(1, 1) = (float)cos(rad);

    return m;

}

FPt MtrxUtil::TransformPtbyHomography(FPt in, Mat& homography) {

/*     Logger("TransformPt FPT start ");
    Logger("... %f %f ", in.x, in.y);
 */
/*     for (int i = 0; i < homography.rows; i++)
        for (int j = 0; j < homography.cols; j++)
            Logger("Check homo in transform [%d][%d] %f ", i, j, homography.at<double>(i, j));

 */    
    Mat mcenter(3, 1, CV_64F);
    mcenter.at<double>(0) = (double)in.x;
    mcenter.at<double>(1) = (double)in.y;    
    mcenter.at<double>(2) = 1;
    Mat mret = homography * mcenter;

    double newx = mret.at<double>(0) / mret.at<double>(2);
    double newy = mret.at<double>(1) / mret.at<double>(2);    

/*     double newx = mret.at<double>(0);
    double newy = mret.at<double>(1);
 */    
    FPt newpt = FPt((float)newx, (float)newy);

    return newpt;
}

Pt MtrxUtil::TransformPtbyHomography(Pt* in, Mat& homography) {

/*     Logger("TransformPt PT start ");
    Logger("in %d %d ", in->x, in->y);
    Logger("homography ");
    for (int i = 0; i < homography.rows; i++)
        for (int j = 0; j < homography.cols; j++)
            Logger("[%d][%d] %f ", i, j, homography.at<double>(i, j));
 */

    Mat mcenter(3, 1, CV_64F);
    mcenter.at<double>(0) = (double)in->x;
    mcenter.at<double>(1) = (double)in->y;    
    mcenter.at<double>(2) = 1;
    Mat mret = homography * mcenter;

    double newx = mret.at<double>(0) / mret.at<double>(2);
    double newy = mret.at<double>(1) / mret.at<double>(2);    
    Pt newpt = Pt((int)newx, (int)newy);
    return newpt;
}

Point2f MtrxUtil::TransformPtbyHomography(Point2f* in, Mat& homography) {

/*     Logger("TransformPt Point2f  start ");
    Logger("in %f %f ", in->x, in->y);
 */
    Mat mcenter(3, 1, CV_64F);
    mcenter.at<double>(0) = in->x;
    mcenter.at<double>(1) = in->y;    
    mcenter.at<double>(2) = 1;
    Mat mret = homography * mcenter;

    double newx = mret.at<double>(0) / mret.at<double>(2);
    double newy = mret.at<double>(1) / mret.at<double>(2);    
    Point2f newpt = Point2f((float)newx, (float)newy);
    Logger(" TransformPt Point2f end %f %f ", newx, newy);
    return newpt;
}

int MtrxUtil::TransformPtsbyHomography(Pt* in, Mat& homography, int cnt) {
    
    Logger(" Start Transform PTS ");
    for(int i = 0 ; i < cnt ; i ++) {
        Logger(" check .. in[i] %d %d ", i, in[i].x, in[i].y );
        Pt tp = TransformPtbyHomography(&in[i], homography);
        in[i].x = tp.x;
        in[i].y = tp.y;
    }

    return ERR_NONE;
}

FPt MtrxUtil::TransformPtbyAffine(FPt in, Mat& aff) {

    vector<Point2f>out;
    vector<Point2f> ps;    
    Point2f p = Point2f((float)in.x, (float)in.y);
    ps.push_back(p);

    transform(ps, out, aff);
/*     double newx = mret.at<double>(0) / mret.at<double>(2);
    double newy = mret.at<double>(1) / mret.at<double>(2);    
 */    
    FPt newpt;
    for(int i = 0 ; i < out.size(); i ++) {
        Logger(" affine pt [%d] %f %f", i, out[i].x, (float)out[i].y);
        newpt = FPt((float)out[i].x, (float)out[i].y);
    }
    return newpt;
}

int MtrxUtil::Hamming(uchar* arr1, uchar* arr2, int size) {
    int distance = 0;

    for(int i = 0; i < size ; i ++) {
/*         int result = arr1[i]^arr2[i];
        int count = 0;

        while(result > 0){
            count += result & 1;
            result >> =1;
        }
 */        

#if defined _MAC_
        distance += __builtin_popcount(arr1[i]^arr2[i]); //gcc
#endif
#if defined _WIN_ ||  _WINDOWS
        //Logger("windows popcnt ! ");
        distance += _mm_popcnt_u64(arr1[i]^arr2[i]); // for msvc <nmmintrin.h>

#endif
    }

    return distance;    
}

int MtrxUtil::SignValue(double x) {

  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return x;
}