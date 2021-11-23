/*****************************************************************************
*                                                                            *
*                            ESMonnector.cpp         								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : Connector
    Author(S)       : Me Eunkyung
    Created         : 03 Nov 2021

    Description     : Data Coverter from both direction : ESM <-> Genesis
    Notes           : 
*/


#include "ESMConnector.h"
#include "DefData.hpp"

ESMConnector::ESMConnector(){

#if defined _WIN_ || _WINDOWS
    string PrjDir = "";
    string LogDir = "";
    //int stat = mkdir(PrjDir.c_str());
    //System("mkdir C:\\NewDir");
#endif    
    dl = Dlog();
}

ESMConnector::~ESMConnector() {

}

int ESMConnector::CheckDataValidity(MppcPoints* mpp, std::string dsc_id, int _scale) {
    dl.SetLogFilename(dsc_id);
    dl.Logger("ESMConnector CheckDataValidity start. %s ", dsc_id.c_str());
    int result = -1;
    int _W = 3840;
    int _H = 2160;
    scale = _scale;
    if( mpp->Width == 1920) {
        _W = 1920; _H = 1080;
    }

    if (mpp->framenum < 0) {
        dl.Logger("ESMConnector CheckDataValidity framenum < 0 %s ", dsc_id.c_str());
        return result;
    }
    //if (mpp->Width != _W || mpp->Height != _H) {
   //     dl.Logger("ESMConnector CheckDataValidity Width %d  Height %d ?  %s ", mpp->Width, mpp->Height, dsc_id.c_str());
   //     return result;
    //}
    if (mpp->pts_3d.Point1.X/ scale < 0 || mpp->pts_3d.Point1.X / scale > _W
        || mpp->pts_3d.Point1.Y / scale < 0 || mpp->pts_3d.Point1.Y / scale > _H
        || mpp->pts_3d.Point2.X / scale < 0 || mpp->pts_3d.Point2.X / scale > _W
        || mpp->pts_3d.Point2.Y / scale < 0 || mpp->pts_3d.Point2.Y / scale > _H
        || mpp->pts_3d.Point3.X / scale < 0 || mpp->pts_3d.Point3.X / scale > _W
        || mpp->pts_3d.Point3.Y / scale < 0 || mpp->pts_3d.Point3.Y / scale > _H
        || mpp->pts_3d.Point4.X / scale < 0 || mpp->pts_3d.Point4.X / scale > _W
        || mpp->pts_3d.Point4.Y / scale < 0 || mpp->pts_3d.Point4.Y / scale > _H
        ) {
        dl.Logger("ESMConnector CheckDataValidity mpp->pts_3d.Point odd ?  %s ", dsc_id.c_str());
        dl.Logger(" %d %d %d %d %d %d %d %d ", mpp->pts_3d.Point1.X,
            mpp->pts_3d.Point1.Y,
            mpp->pts_3d.Point2.X,
            mpp->pts_3d.Point2.Y,
            mpp->pts_3d.Point3.X,
            mpp->pts_3d.Point3.Y,
            mpp->pts_3d.Point4.X,
            mpp->pts_3d.Point4.Y);
        return result;
    }
        
    return ERR_NONE;
}
#if defined GPU
int ESMConnector::Recalibration(cv::cuda::GpuMat& ref, cv::cuda::GpuMat& cur, MppcPoints* mpp, string dsc_id) 
#else
int ESMConnector::Recalibration(Mat& ref, Mat& cur, MppcPoints* mpp, string dsc_id) 
#endif
{
    dl.SetLogFilename(dsc_id);
    dl.Logger("ESMConnector Recalibration Start. ");    
    int result = 0;
    FPt in_pt[4]; 
    FPt out_pt[4];  
    ConvertFromESMtoGenesis(mpp, in_pt);
    try {
        result = genesis->ExecuteClient(ref, cur, in_pt, out_pt, dsc_id);
    } catch (std::exception e) {        
        //return EXECUTE_CLIENT_EXCEPTION;
        dl.Logger("ExecuteClient error. result %d ", result);
    }

    if (result > 0 ) {
        ConvertFromGenesistoESM(out_pt, mpp);
        dl.Logger("Retruend score %d ", result);
        mpp->recalibration_status = 1;

        mpp->recalibration_score = (double)result;
        dl.Logger("Retruend score %f ", mpp->recalibration_score);
    } else {
        dl.Logger("Recalibraiton Execute Error  %d ", result);
        mpp->recalibration_status = -1;  
        mpp->recalibration_score = result;
    }

    return result;;
}


int ESMConnector::Start(int width) {
    dl.Logger("ESMConnector Start() - Extractor construct call with Width %d ", width);        
    genesis = new Extractor(width, true);
    return ERR_NONE;    
}

int ESMConnector::Finish() {
    dl.Logger("ESMConnector Finish. Call Extractor Destroyer ");
    genesis->~Extractor();
    return ERR_NONE;
}

int ESMConnector::ConvertFromESMtoGenesis(MppcPoints* mpp, FPt* in_pt){ 

    in_pt[0].x = mpp->pts_3d.Point1.X / scale;
    in_pt[0].y = mpp->pts_3d.Point1.Y / scale;
    in_pt[1].x = mpp->pts_3d.Point2.X / scale;
    in_pt[1].y = mpp->pts_3d.Point2.Y / scale;
    in_pt[2].x = mpp->pts_3d.Point3.X / scale;
    in_pt[2].y = mpp->pts_3d.Point3.Y / scale;
    in_pt[3].x = mpp->pts_3d.Point4.X / scale;
    in_pt[3].y = mpp->pts_3d.Point4.Y / scale;

    for(int i = 0; i < 4 ; i ++) {
        dl.Logger("Input in_pt in Connector P[%d] %f %f %",
        i, in_pt[i]. x, in_pt[i].y);
    }
    return ERR_NONE;

}

int ESMConnector::ConvertFromGenesistoESM(FPt* out, MppcPoints* mpp) {

    mpp->pts_3d.Point1.X = out[0].x * scale;
    mpp->pts_3d.Point1.Y = out[0].y * scale;
    mpp->pts_3d.Point2.X = out[1].x * scale;
    mpp->pts_3d.Point2.Y = out[1].y * scale;
    mpp->pts_3d.Point3.X = out[2].x * scale;
    mpp->pts_3d.Point3.Y = out[2].y * scale;
    mpp->pts_3d.Point4.X = out[3].x * scale;
    mpp->pts_3d.Point4.Y = out[3].y * scale;

    for(int i = 0; i < 4 ; i ++) {
        dl.Logger("Retruend out_pt in Connector P[%d] %f %f %",
        i, out[i].x * scale, out[i].y * scale);
    }

    return ERR_NONE;

}