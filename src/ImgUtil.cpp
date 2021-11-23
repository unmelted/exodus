
/*****************************************************************************
*                                                                            *
*                            ImgUtil          								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : ImgUtilr.cppp
    Author(S)       : Me Eunkyung
    Created         : 23 Oct 2021

    Description     : ImgUtil.cpp
    Notes           : Utility related image.
*/
#define _CRT_SECURE_NO_WARNINGS
#include "ImgUtil.hpp"

using namespace std;
using namespace cv;

ImgUtil::ImgUtil() {
    mtrx = MtrxUtil();
}

void ImgUtil::SaveImage(SCENE *sc, int type, SCENE* sc2, PARAM* p, int opt)
{
    if (type == 0)
        return;

    if (type == 1)
    {
        Mat img;
        char filename[50] = { 0, };
        drawKeypoints(sc->img, sc->ip, img);
#if defined _WIN_ || _WINDOWS       
        sprintf(filename, "recalibration\\saved\\%s_%d_keypoint.png", fname.c_str(), sc->id);
#else
        sprintf(filename, "recalibration/saved/%s_%d_keypoint.png", fname.c_str(), sc->id);
#endif
        imwrite(filename, img);        
    }
    else if (type == 2)
    {
        Mat img;
        char filename[50] = { 0, };
        img = sc->mask_img;
#if defined _WIN_ || _WINDOWS       
        sprintf(filename, "recalibration\\saved\\%s_%d_masking.png", fname.c_str(), sc->id);
#else
        sprintf(filename, "recalibration/saved//%s_%d_keypoint.png", fname.c_str(), sc->id);        
#endif

        imwrite(filename, img);                
    }
    else if (type == 3)
    {
        Mat img;
        char filename[50] = { 0, };
        drawKeypoints(sc->img, sc->ip, img);
        for (int i = 0; i < 4; i++)
        {
            circle(img,
                    Point((int)sc2->four_fpt[i].x/p->p_scale, (int)sc2->four_fpt[i].y/p->p_scale), 6, Scalar(255), -1);

        }        
#if defined _WIN_ || _WINDOWS               
        sprintf(filename, "recalibration\\saved\\%s_%d_keypoint.png", fname.c_str(), sc->id);
#else
        sprintf(filename, "recalibration/saved//%s_%d_keypoint.png", fname.c_str(), sc->id);
#endif        
        imwrite(filename, img);        
    }
    else if (type == 4) {
        Mat img;
        char filename[50] = { 0, };

        vector<Scalar> color;
        color.push_back(Scalar(0, 0, 255));
        color.push_back(Scalar(0, 255, 0));
        color.push_back(Scalar(255, 0, 0));
        color.push_back(Scalar(255, 255, 0));
        color.push_back(Scalar(255, 0, 255));
        color.push_back(Scalar(0, 255, 255));        

        drawKeypoints(sc->img, sc->ip, img);
        for (int i = 0; i < 4; i++)       {
            circle(img,
                    Point((int)sc2->four_fpt[i].x/p->p_scale, (int)sc2->four_fpt[i].y/p->p_scale), 6, Scalar(255), -1);
        }

        Point2f tp1 = Point2f(sc->projected_normal.at<Point2f>(0).x/p->p_scale, 
                sc->projected_normal.at<Point2f>(0).y/p->p_scale);

        for (int i = 1 ; i < sc->projected_normal.rows; i ++) {
            Logger("projected point %f %f ", sc->projected_normal.at<Point2f>(i).x,
                    sc->projected_normal.at<Point2f>(i).y);
            Point2f tp2 = Point2f(sc->projected_normal.at<Point2f>(i).x/p->p_scale, 
                    sc->projected_normal.at<Point2f>(i).y/p->p_scale);

            line(img, tp1, tp2, color[i], 2);
        }

#if defined _WIN_ || _WINDOWS       
        sprintf(filename, "recalibration\\saved/\\%s_%d_keypoint+normal.png", fname.c_str(), sc->id);
#else
        sprintf(filename, "recalibration/saved//%s_%d_masking.png", fname.c_str(), sc->id);
#endif

        imwrite(filename, img);        
    }
    else if (type == 5) {
        Mat img;
        char filename[50] = { 0, };
        img = sc->pyramid[0];
#if defined _WIN_ || _WINDOWS                   
        sprintf(filename, "recalibration\\saved\\%d_pyramid_check0.png", sc->id);        
#else
        sprintf(filename, "recalibration/saved//%d_pyramid_check0.png", sc->id);        
#endif        
        imwrite(filename, img);        
        img = sc->pyramid[1];
#if defined _WIN_ || _WINDOWS                           
        sprintf(filename, "recalibration\saved\\%d_pyramid_check1.png", sc->id);        
#else        
        sprintf(filename, "recalibration/saved//%d_pyramid_check1.png", sc->id);        
#endif        
        imwrite(filename, img);        
        img = sc->pyramid[2];
#if defined _WIN_ || _WINDOWS                           
        sprintf(filename, "recalibration\\saved\\%d_pyramid_check2.png", sc->id);
#else
        sprintf(filename, "recalibration/saved//%d_pyramid_check2.png", sc->id);
#endif                
        imwrite(filename, img);                

    }
    else if (type == 6) {
        Mat img;
        char filename[50] = { 0, };
        for(int step = p->pyramid_step - 1 ; step >= 0; step --){ 
            int scl = p->pyramid_scale[step];
            int cnt = sc->pyramid_ip[step].size(); 
            img = sc->pyramid[step].clone();            
            for (int i = 0; i < cnt; i++) {
                circle(img,
                    Point((int)sc->pyramid_ip[step][i].pt.x, (int)sc->pyramid_ip[step][i].pt.y), 
                    1, Scalar(125, 125,0), -1);

            }
#if defined _WIN_ || _WINDOWS
            sprintf(filename, "recalibration\\saved\\%s_%d_%d_pr_sel_point.png", fname.c_str(), sc->id, step);
#else            
            sprintf(filename, "recalibration/saved//%s_%d_%d_pr_sel_point.png", fname.c_str(), sc->id, step);
#endif            
            imwrite(filename, img);
        }

    } else if (type == 7) {
        Mat img;
        char filename[50] = { 0, };
        int step = opt;
        img = sc->pyramid[step].clone();
        int scl = p->pyramid_scale[step];
        for (int i = 0; i < 4; i++) {
/*             Logger(" pyramid pair %d %d " , (int)sc->pyramid_pair[step][i].query.x, 
                                    (int)sc->pyramid_pair[step][i].query.y);
 */            
            circle(img,
                Point((int)sc->pyramid_pair[step][i].query.x/scl, (int)sc->pyramid_pair[step][i].query.y/scl), 2, Scalar(255), -1);
#if defined _WIN_ || _WINDOWS
            sprintf(filename, "recalibration\\saved\\%s_%d_%d_pr2_best_point.png", fname.c_str(), sc->id, step);
#else            
            sprintf(filename, "recalibration/saved//%s_%d_%d_pr2_best_point.png", fname.c_str(), sc->id, step);
#endif            
            imwrite(filename, img);
        }
    }
    else if (type == 8) {
            Mat img;
        char filename[50] = { 0, };
        for (int step = 2; step >= 0; step--) {
            img = sc->pyramid[step].clone();
#if defined _WIN_ || _WINDOWS
            sprintf(filename, "recalibration\\saved\\%s_%d_train_create_err.png", fname.c_str(), step);
#else            
            sprintf(filename, "recalibration/saved//%s_%d_train_create_err.png", fname.c_str(), step);
#endif            
            imwrite(filename, img);
        }
    }
 
}
vector<Mat> ImgUtil::LoadImages(const string &path, vector<string>* dsc_id)
{
#if defined _MAC_    
    const int FK = 3500;
    const int FHD = 1900;
    vector<string>image_paths;

    namespace fs = std::__fs::filesystem;

    for (const auto &entry : fs::directory_iterator(path)) {
        if (fs::is_regular_file(entry) &&
           (entry.path().extension().string() == ".png" ||
            entry.path().extension().string() == ".jpg")) {
            image_paths.push_back(entry.path().string());
        }
    }

    sort(begin(image_paths), end(image_paths), less<string>());
    vector<Mat> images;
    int len = path.length();
    for (const string &ip : image_paths) {        
        images.push_back(imread(ip));
        string dsc = ip.substr(len, len + 6);
        dsc = dsc.substr(0, 7);
        dsc_id->push_back(dsc);
        Logger("Read image : %s , desc_id %s ", ip.c_str(), dsc.c_str());
    }

#else
    vector<Mat> images;
#endif    
    return images;
}

void ImgUtil::SimpleColorBalance(Mat& in, Mat& out, float percent){

    float half_percent = percent / 200.0f;

    vector<Mat> tmpsplit; split(in,tmpsplit);
    for(int i=0;i<3;i++) {
        //find the low and high precentile values (based on the input percentile)
        Mat flat; tmpsplit[i].reshape(1,1).copyTo(flat);
        cv::sort(flat,flat, SORT_EVERY_ROW + SORT_ASCENDING);
        int lowval = flat.at<uchar>(cvFloor(((float)flat.cols) * half_percent));
        int highval = flat.at<uchar>(cvCeil(((float)flat.cols) * (1.0 - half_percent)));
        cout << lowval << " " << highval << endl;
        
        //saturate below the low percentile and above the high percentile
        tmpsplit[i].setTo(lowval,tmpsplit[i] < lowval);
        tmpsplit[i].setTo(highval,tmpsplit[i] > highval);
        
        //scale the channel
        normalize(tmpsplit[i],tmpsplit[i],0,255,NORM_MINMAX);
    }
    merge(tmpsplit,out);      
}

int ImgUtil::AdjustImage(SCENE* sc, ADJST adj) {

    Size sz = Size(sc->ori_img.cols, sc->ori_img.rows);
    double angle = adj.angle + 90;
    double rad = angle * M_PI / 180.0;
    Logger("Adjust Image angle %f rad %f ", angle, rad);

    //Mat flipm = Mat::eye(3, 3, CV_32FC1);
    Mat mrot = mtrx.GetRotationMatrix(rad, adj.rotate_centerx, adj.rotate_centery);
    Mat mscale = mtrx.GetScaleMatrix(adj.scale, adj.scale, adj.rotate_centerx, adj.rotate_centery);
    Mat mtran = mtrx.GetTranslationMatrix(adj.trans_x, adj.trans_y);
    Mat mscaleout = mtrx.GetScaleMatrix(1, 1);

    Mat mfm = mscaleout * mtran * mscale * mrot;
    Logger("5 assembed mfm matrix ");
    for (int i = 0; i < mfm.rows; i++)
        for (int j = 0; j < mfm.cols; j++)
            Logger("[%d][%d] %f ", i, j, mfm.at<float>(i, j));

    Mat mret = mfm(Rect(0, 0, 3, 2));
    Logger("6 mret = submatrix of mfm");
    for (int i = 0; i < mret.rows; i++)
        for (int j = 0; j < mret.cols; j++)
            Logger("[%d][%d] %f ", i, j, mret.at<float>(i, j));

    Mat final;
    warpAffine(sc->ori_img, final, mret, Size(sc->ori_img.cols, sc->ori_img.rows));

    static int index = 0;
    char filename[30] = {
        0,
    };
    sprintf(filename, "recalibration/saved//%s_%2d_perspective.png", fname.c_str(), index);
    imwrite(filename, final);

    return ERR_NONE;
}

void ImgUtil::ColorCorrection(Mat& ref, Mat& src, Mat& out) {
    const float HISTMATCH = 0.000001;
    const int HIST_MAX = 256;    
    Mat ref_hist, src_hist;
    Mat ref_hist_acc, src_hist_acc;
    double min, max;
    int channel[] = { 0 };
    const int histSize[] = { HIST_MAX };
    float graylevel[] = {0, HIST_MAX};
    const float* range[] = { graylevel };

    calcHist(&ref, 1, 0, Mat(),
        ref_hist, 1, histSize, range, true, false);
    calcHist(&src, 1, 0, Mat(),
        src_hist, 1, histSize, range, true, false);

    minMaxLoc(ref_hist, &min, &max);
    normalize(ref_hist, ref_hist, min/max , 1.0, NORM_MINMAX);
    minMaxLoc(src_hist, &min, &max);
    normalize(src_hist, src_hist, min/max , 1.0, NORM_MINMAX);

    ref_hist.copyTo(ref_hist_acc);
    src_hist.copyTo(src_hist_acc);

    float* ref_cdf = ref_hist_acc.ptr<float>();
    float* src_cdf = src_hist_acc.ptr<float>();

    for(int i = 1 ; i < HIST_MAX; i ++) {
        ref_cdf[i] += ref_cdf[i -1];
        src_cdf[i] += src_cdf[i -1];        
    }

    minMaxLoc(ref_hist_acc, &min, &max);
    normalize(ref_hist_acc, ref_hist_acc, min/max , 1.0, NORM_MINMAX);
    minMaxLoc(src_hist_acc, &min, &max);
    normalize(src_hist_acc, src_hist_acc, min/max , 1.0, NORM_MINMAX);

    Mat lut(1, 256, CV_8UC1);
    uchar *M = lut.ptr<uchar>();
    uchar last = 0;

    for( int j = 0 ; j < src_hist_acc.rows; j ++) {
        float F1 = src_cdf[j];
        
        for (uchar k = last; k < ref_hist_acc.rows; k++) {
            float F2 = ref_cdf[k];
            if (abs(F2- F1) < HISTMATCH || F2 > F1) {
                M[j] = k;
                last = k;
                break;
            }
        }
    }
    
    Logger("color correction .. ");
    LUT(src, lut, out);
    // imwrite("recalibration/saved//cc_ref.png", ref);
    // imwrite("recalibration/saved//cc_src.png", src);    
    // imwrite("recalibration/saved//cc_out.png", out); 
}
