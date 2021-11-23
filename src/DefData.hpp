  
/*****************************************************************************
*                                                                            *
*                            DefData         								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : DefData.hpp
    Author(S)       : Me Eunkyung
    Created         : 17 Sep 2021

    Description     : Data structure definition for Kelly Genesis
    Notes           : 
*/

#pragma once

#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <opencv2/calib3d/calib3d.hpp>

#include <string>
#include <vector>
#include <stdio.h>
#include <string.h>
#include "common/TimeUtil.hpp"
#include "common/os_api.hpp"
#if defined _WIN_ || _WINDOWS
#define _IMGDEBUG
#define _DEBUG
#define GPU

#include <opencv2/core/cuda.hpp>
#include <opencv2/cudawarping.hpp>
#include <opencv2/cudaimgproc.hpp>
#endif
using namespace cv;
using namespace std;

typedef enum _err {
    ERR_NONE = 0,
    
    TRAIN_CREATE_FEATURE_ERR = -10,
    PYRAMID_MATCH_NOT_FOUND_NEW_POINT = -20,
    PYRAMID_MATCH_FOUND_NEW_POINT_INTEGRITY_FAIL = -21,

    EXECUTE_CLIENT_EXCEPTION    = -30,
} ERR;


typedef struct _Pt
{    
    int x;
    int y;
    int z;
    _Pt(int a, int b ) {
        x = a; y = b; z = 0;
    };
    _Pt() { x = 0, y = 0; z = 0;};
} Pt;

typedef struct _FPt
{    
    float x;
    float y;
    float z;
    _FPt(float a, float b ) {
        x = a; y = b; z = 0.0;
    };
    _FPt() { x = 0.0, y = 0.0; z = 0.0;};
} FPt;

typedef struct _cc {
    Pt center;
    int radius;
} Cr;

typedef struct _matchpair {
    FPt train;
    FPt query;
    int distance;
    int pyramid_sacle;
    int kernel_size;

    _matchpair(FPt tpt, FPt qpt, int dist, int scl, int ksize ) {
        train = tpt;
        query = qpt;
        distance = dist;
        pyramid_sacle = scl;
        kernel_size = ksize;
    };
    
} MATCHPAIR;

typedef struct _maindata {
    int id = 0;
    char filename[100];
    Mat img;
#if defined GPU
    cv::cuda::GpuMat ori_img;
#else
    Mat ori_img;
#endif
    Mat mask_img;

    //for pyramid matching
    Mat pyramid[3];
    vector<KeyPoint>pyramid_ip[3];
    int pyramid_ip_per_pt[3];    
    Mat pyramid_desc[3];
    vector<MATCHPAIR>pyramid_pair[3];
    
    //Pt four_pt[4];
    FPt four_fpt[4];    
    FPt center;

    Mat rot_matrix;
    Mat trans_matrix;
    Point2f normal_vec[2];
    Mat projected_normal;

    double rod_norm;
    double rod_degree;
    Mat rod_rotation_matrix;

    vector<KeyPoint> ip;
    Mat desc;
    Mat matrix_fromimg;
    Mat matrix_scaledfromimg;

} SCENE;

typedef struct _adj { 
    double angle;
    double scale;
    double trans_x;
    double trans_y;
    double rotate_centerx;
    double rotate_centery;
    Rect rect;

} ADJST;

typedef struct _PARAM {
    int ground;     //Groud type
    int calibration_type;
    int roi_type;
    int roi_count;
    int masking_type;
    int circle_fixedpt_radius;
    int circle_fixedpt_radius_2nd;
    int match_type;
    int submatch_type;

    int pyramid_step;
    int pyramid_scale[3];
    int pyramid_patch[3];
    int stride[3];
    int base_kernel;
    int desc_kernel[3];
    float best_cut;
    float pixel_diff_cut;
    
    int blur_ksize;
    float blur_sigma;
    int desc_byte;
    bool use_ori;
    int nms_k;
    int fast_k;
    int minx;
    int p_scale;


    int pwidth;
    int pheight;
    float sensor_size;
    float focal;
    float* camera_matrix;
    float* skew_coeff;
    float normal[2][3];
    
    SCENE* world;    // World Coord 4 point
    Pt* region;     // Point array for polygon (ROI)
    Cr* circles;
    Pt* moved_region;

    void initialize() { 
        world = NULL;    // World Coord 4 point
        region = NULL;     // Point array for polygon (ROI)
        circles = NULL;
        moved_region = NULL;
    };

} PARAM;

enum _match {
    PLAIN_MATCH    = 1,
    PYRAMID_MATCH  = 4,
};

enum _submatch {
    SUBMATCH_NONE   = -1,
    BEST_MATCH     = 1,
    KNN_MATCH      = 2,
    SPLIT_MATCH    = 3,
};

enum _preset_calibration_type {
    PRESET_NONE_2D      = 0,
    PRESET_NONE_3D      = 1,
    RECALIBRATION_2D    = 2,
    RECALIBRATION_3D    = 3, 
};

enum roitype {
    POLYGON     = 1,
    CIRCLE      = 2,
    RECTANGLE   = 3,
};

enum roi_base_type {
    FOUR_POINT_BASE     = 1,
    USER_INPUT_CIRCLE   = 2,
    INNER_2POINT_BASE   = 3,
};

enum _coordwd { 
    FIRST_MATCH = 0,
    NORMAL_VECTOR_CAL = 1,
};

enum _keypoint_array {
    PLANE       = 0,
    CIRCULAR    = 1,
};

typedef enum _groundtype
{
    BaseballHome,
    BaseballGround,
    BasketballHalf,
    BasketballGround,
    Boxing,
    IceLinkHalf,
    IceLink,
    SoccerHalf,
    Soccer,
    Taekwondo,
    TennisHalf,
    Tennis,
    Ufc,
    VolleyballHalf,
    VolleyballGround,
    Football

} GROUNDTYPE;