
/*****************************************************************************
*                                                                            *
*                            Covenant      							    	 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : Covenant.Hpp
    Author(S)       : Me Eunkyung
    Created         : 24 Nov 2021

    Description     : Covenant.Hpp
    Notes           : Auto calibration concept main procedure.
*/
#define _USE_MATH_DEFINES
#include <cmath>
#include "DefData.hpp"
#include "common/Pip.hpp"
#include "MtrxUtil.hpp"
#include "ExpUtil.hpp"
#include "ImgUtil.hpp"


using namespace std;
using namespace cv;

class Covenant {

public :
    Covenant(int width, int cnt , int* gline, int* iline);
    Covenant(int width = 3840, bool _use_gpu = false);    
    ~Covenant();
    int Execute();
    
    int Process(int step, int* buffers);

    PARAM* p;
    MtrxUtil mtrx;
    ExpUtil genutil;
    ImgUtil imgutil;
    Dlog dl;
    
    vector<string>dsc_id;
    vector<Mat>imgs;
    vector<SCENE>cal_group;


private :
    string imgset;
    TIMER* t;
    bool verify_mode = false;
    bool use_gpu = false;
    SCENE* ground;
    SCENE* base;

    SCENE* cur_train = 0;
    SCENE* cur_query = 0;

    void InitializeData(int width = 3840, int cnt = 0, int* gline = 0, int* iline = 0);
    
//    Mat ProcessImages(Mat& img);
    int ProcessImages(SCENE* sc);
    int ImageMasking(SCENE* sc);
    int GetFeature(SCENE* sc);

    vector<KeyPoint> KeypointMasking(vector<KeyPoint>* oip);

    void SetCurTrainScene(SCENE* sc) { cur_train = sc; };
    void SetCurQueryScene(SCENE* sc) { cur_query = sc; };
    int Match();    
    int MatchPlain();
    int MatchPyramid();    
    int MatchSplit(vector<Point2f> m_train, vector<Point2f>m_query);
    int MatchVerify();

    vector<DMatch> RefineMatch(vector<DMatch> good);
    vector<DMatch> RemoveOutlier(vector<DMatch> matches);

    int FindBaseCoordfromWd(int mode = 0);
    int FindHomographyP2P(); 

    int PostProcess();

    int DecomposeHomography();
    ADJST CalAdjustData();


    int WarpingStep1();

};
