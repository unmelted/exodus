
/*****************************************************************************
*                                                                            *
*                            Extractor      								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : Extractor.Hpp
    Author(S)       : Me Eunkyung
    Created         : 17 Sep 2021

    Description     : Extractor.Hpp
    Notes           : Feature Extractor From Image.
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

class Extractor {

public :
    Extractor(string& imgset, int cnt , int* roi);
    Extractor(int width = 3840, bool _use_gpu = false);    
    ~Extractor();
    int Execute();
    
//    int ExecuteSever(string ref_path, string cur_path, string ref_pts_path, string& out_pts_path);
#if defined GPU
    int ExecuteClient(cv::cuda::GpuMat ref_file, cv::cuda::GpuMat cur_file, FPt* in_pt, FPt* out_pt, string dsc_id);
#else
    int ExecuteClient(Mat ref_file, Mat cur_file, FPt* in_pt, FPt* out_pt, string dsc_id);
#endif
    PARAM* p;
    MtrxUtil mtrx;
    ExpUtil genutil;
    ImgUtil imgutil;
    Dlog dl;
    
    vector<string>dsc_id;
    vector<Mat>imgs;
    vector<SCENE>cal_group;


private :
    TIMER* t;
    bool verify_mode = false;
    bool use_gpu = false;
    SCENE* cur_train = 0;
    SCENE* cur_query = 0;

    int LoadConfig();
    int UpdateConfig();    
    void InitializeData(int width = 3840, int cnt = 0, int* roi = 0);
    
//    Mat ProcessImages(Mat& img);
    int ProcessImages(SCENE* sc);
    int ImageMasking(SCENE* sc);
    int GetFeature(SCENE* sc);
    int CreateFeature(SCENE* sc, bool train = false, bool query = false, int step = -1);    

    vector<KeyPoint> KeypointMasking(vector<KeyPoint>* oip);

    void SetCurTrainScene(SCENE* sc) { cur_train = sc; };
    void SetCurQueryScene(SCENE* sc) { cur_query = sc; };
    int Match();    
    int MatchPlain();
    int MatchPyramid();    
    int MatchSplit(vector<Point2f> m_train, vector<Point2f>m_query);
    int MatchVerify();
    float ncc(int max_index, Mat _h);

    vector<DMatch> RefineMatch(vector<DMatch> good);
    vector<DMatch> RemoveOutlier(vector<DMatch> matches);

    int FindBaseCoordfromWd(int mode = 0);
    int FindHomographyP2P(); 

    int PostProcess();

    void NormalizePoint(SCENE* sc, int maxrange);
    int DecomposeHomography();
    ADJST CalAdjustData();


    int WarpingStep1();

};
