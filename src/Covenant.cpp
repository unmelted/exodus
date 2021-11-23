
/*****************************************************************************
*                                                                            *
*                            Covenant       								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : Covenant.Cpp
    Author(S)       : Me Eunkyung
    Created         : 24 Nov 2021

    Description     : Covenant.Cpp
    Notes           : Auto calibration concept main procedure.
*/
#define _CRT_SECURE_NO_WARNINGS
#include "Covenant.hpp"

using namespace std;
using namespace cv;

Covenant::Covenant(int width, bool _use_gpu) {
    if(_use_gpu) {
        use_gpu = true;
    }
    mtrx = MtrxUtil();
    genutil = ExpUtil();
    imgutil = ImgUtil();
    dl = Dlog();
    t = new TIMER();    
    InitializeData(width);
}

Covenant::Covenant(int width, int cnt, int *gline, int *iline)
{
    mtrx = MtrxUtil();
    genutil = ExpUtil();
    imgutil = ImgUtil();
    dl = Dlog();    
    t = new TIMER();    

    dl.SetLogFilename("TEST");
    InitializeData(3840, cnt, gline, iline);
    imgset = "image/new/";
    imgs = imgutil.LoadImages(imgset, &dsc_id);

}

Covenant::~Covenant()
{
    dl.Logger("finish 1 ");
    delete ground;
    delete base;

    if(p->region != NULL)
        g_os_free(p->region);
    dl.Logger("finish 2 ");        
    if(p->circles != NULL)        
        g_os_free(p->circles);
    dl.Logger("finish 3 ");        
    if(p->moved_region != NULL)        
        g_os_free(p->moved_region);
    delete p->world;
    dl.Logger("finish 4 ");    
    g_os_free(p->camera_matrix);
    g_os_free(p->skew_coeff);
    g_os_free(p);
    dl.Logger("finish 5 ");
    delete t;
    dl.Logger("finish 6 ");        
//  meminfo();    
}


void Covenant::InitializeData(int width, int cnt, int* gline, int* iline)
{
    p = (PARAM *)g_os_malloc(sizeof(PARAM));
    p->initialize();
    int index = 1;
    ground = new SCENE();
    base = new SCENE();
    dl.Logger(" %d %d %d %d ", gline[1], gline[2], gline[3], gline[4]);
    dl.Logger(" %d %d %d %d ", iline[1], iline[2], iline[3], iline[4]);        

    for(int i = 0 ; i < cnt ; i=i+4) {
        LINE l1(Pt(gline[i], gline[i+1]), Pt(gline[i+2], gline[i+3]));
        LINE l2(Pt(iline[i], iline[i+1]), Pt(iline[i+2], iline[i+3]));
        ground->in_line.push_back(l1);
        base->in_line.push_back(l2);
    }

}

int Covenant::Process(int step, int* buffers){

    dl.Logger("Process start step %d ", step);
    return 1;
}
    
int Covenant::Execute() {

    int index = 0;
    int ret = -1;

    TIMER* all;
    all = new TIMER();    
    StartTimer(all);
    
    for (int i = 0 ; i < imgs.size(); i ++)
    {
        StartTimer(t);
        SCENE sc;
        sc.id = index;
#if not defined GPU
        sc.ori_img = imgs[i];
#endif
        ProcessImages(&sc);

        if (index == 0)
        {
            sc.id = 0;

        }

        if(p->match_type == PLAIN_MATCH) {
            ImageMasking(&sc);
            ret = GetFeature(&sc);
            dl.Logger("[%d] feature extracting  %f ", index, LapTimer(t));
            if (sc.id == 0) {
                cal_group.push_back(sc);
                SetCurTrainScene(p->world);
                SetCurQueryScene(&cal_group[index]);
            }
            else {
                cal_group.push_back(sc);
                SetCurTrainScene(&cal_group[index - 1]);
                SetCurQueryScene(&cal_group[index]);
            }
        }
        else {
            if(sc.id == 0) {
                ret = 0; //CreateFeature(&sc, true, false);
                if( ret != ERR_NONE)
                    return ret;

                dl.Logger("[%d] feature extracting  %f ", index, LapTimer(t));
                cal_group.push_back(sc);
                index++;
                continue;
                //break;
            }
            else {
                ret = 0; //CreateFeature(&sc, false, true);    
                if( ret != ERR_NONE)
                    return ret;

                dl.Logger("[%d] feature extracting  %f ", index, LapTimer(t));
                cal_group.push_back(sc);                
                SetCurTrainScene(&cal_group[0]);                         
                //SetCurTrainScene(&cal_group[index - 1]);
                SetCurQueryScene(&cal_group[index]);
            }
        }
        ret = Match();
        dl.Logger("return from FindHomography------  %d", ret);
        dl.Logger("[%d] match consuming %f ", index, LapTimer(t));        

        if (ret > 0) {
            PostProcess();
#ifdef _IMGDEBUG
//            imgutil.SaveImage(&sc, 4, cur_train, p);
#endif            
        }
        else {
            dl.Logger("Match Pair Fail. Can't keep going.");
        }

        dl.Logger("------- [%d] end  consuming %f ", index, LapTimer(t));

        index++;
        // if (index == 2)
        //     break;
    }

    //Export result to josn file
    //genutil.Export(dsc_id, cal_group, p);
    //genutil.ExportforApp(dsc_id, cal_group, p);
    dl.Logger("All process time..  %f ", LapTimer(all));

    return ERR_NONE;
}

int Covenant::ProcessImages(SCENE* sc) {

}

int Covenant::ImageMasking(SCENE* sc)
{
    dl.Logger("Image masking function start ");
    Mat mask = Mat::zeros(sc->img.rows, sc->img.cols, CV_8UC1);

    if (p->masking_type == FOUR_POINT_BASE)
    {
        for (int i = 0; i < 4; i++)
        {
            if (sc->id == 0) {
                dl.Logger("masking point 1 %f %f ", sc->four_fpt[i].x, sc->four_fpt[i].y);
                circle(mask, Point((int)sc->four_fpt[i].x/p->p_scale, (int)sc->four_fpt[i].y/p->p_scale),
                       int(p->circle_fixedpt_radius), Scalar(255), -1);
            }
            else {
                dl.Logger("masking point 2 %f %f ", cal_group.back().four_fpt[i].x, cal_group.back().four_fpt[i].y);
                circle(mask,
                       Point((int)cal_group.back().four_fpt[i].x/p->p_scale, (int)cal_group.back().four_fpt[i].y/p->p_scale),
                       int(p->circle_fixedpt_radius), Scalar(255), -1);
            }
        }
    }
    else if (p->masking_type == USER_INPUT_CIRCLE) {
        for (int i = 0; i < p->roi_count; i++) {
            dl.Logger("masking point 3 %d %d %d ", p->circles[i].center.x, p->circles[i].center.y, p->circles[i].radius);
            circle(mask,
                   Point(p->circles[i].center.x, p->circles[i].center.y),
                   p->circles[i].radius, Scalar(255), -1);
        }
    }

    sc->mask_img = mask;

#if defined _IMGDEBUG
//    SaveImage(sc, 2);
#endif
    return ERR_NONE;    
}

int Covenant::GetFeature(SCENE *sc) {
    // FAST + BRIEF
    //auto feature_detector = FastFeatureDetector::create(p->fast_k, true, FastFeatureDetector::TYPE_9_16);

    //auto feature_detector = AgastFeatureDetector::create(); //AGAST
    //Ptr<xfeatures2d::MSDDetector>feature_detector = xfeatures2d::MSDDetector::create(9,11,15); //MSDETECTOR
    //Ptr<xfeatures2d::BriefDescriptorExtractor> dscr;
    //dscr = xfeatures2d::BriefDescriptorExtractor::create(p->desc_byte, p->use_ori);

    //AKAZE
    //auto feature_detector = AKAZE::create(false, false, 0.001f, 4, 1, KAZE::DIFF_PM_G1); //KAZE
    auto feature_detector = AKAZE::create();
    Ptr<AKAZE> dscr = feature_detector;

    Mat desc;
    vector<KeyPoint> kpt;
    vector<KeyPoint> f_kpt;

    feature_detector->detect(sc->img, kpt, sc->mask_img);
    dl.Logger("extracted keypoints count : %d", kpt.size());
    //f_kpt = KeypointMasking(&kpt);
    dscr->compute(sc->img, kpt, desc);

    sc->ip = kpt;
    sc->desc = desc;

#if defined _DEBUG
/*     for (int i = 0 ; i < sc->ip.size(); i ++) {
        dl.Logger("Keypoint index %2d id %3d, x %f y %f ",i, sc->ip[i].class_id, 
                            sc->ip[i].pt.x, sc->ip[i].pt.y);
    }
 */
#endif

    return ERR_NONE;
}

vector<KeyPoint> Covenant::KeypointMasking(vector<KeyPoint> *oip)
{
    vector<KeyPoint> ip;
    //dl.Logger("Before masking %d ", oip->size());
    int left = 0;
    int total = 0;
    int del = 0;

    for (int i = 0; i < p->roi_count; i++)
    {
        dl.Logger("roi check %d %d ", p->moved_region[i].x, p->moved_region[i].y);
    }

    for (auto it = oip->begin(); it != oip->end(); it++)
    {
        total++;
        Pt cp(int(it->pt.x), int(it->pt.y));
        int ret = isInside(p->moved_region, p->roi_count, cp);

        if (ret == 0 && it != oip->end())
        {
            del++;
        }
        else if (it->pt.x >= 1760 && it->pt.y >= 930 &&
                 it->pt.x <= 2070 && it->pt.y <= 1190)
            del++;
        else
        {
            ip.push_back(*it);
            left++;
        }
    }

    dl.Logger("new vector  %d. left %d  del %d / total ip %d ", ip.size(), left, del, total);
    return ip;
}

int Covenant::Match() {

    int ret = -1;
    if(p->match_type == PLAIN_MATCH) {
        if (cur_query->id == 0) {
            ret = FindBaseCoordfromWd();
        } else {
            ret = MatchPlain();
        }
    }
    
    return ret;
}

int Covenant::MatchPlain() {

    Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create(DescriptorMatcher::FLANNBASED);
    vector<DMatch> matches;
    vector<DMatch> good;

    if (cur_train->desc.type() != CV_32F || cur_query->desc.type() != CV_32F)
    {
        cur_train->desc.convertTo(cur_train->desc, CV_32F);
        cur_query->desc.convertTo(cur_query->desc, CV_32F);
    }
    dl.Logger("Match start %d %d ", cur_train->ip.size(), cur_query->ip.size());

    if(p->submatch_type == KNN_MATCH) {
        vector<vector<DMatch>> in;        
        const float ratio_thresh = 0.90f;        
        matcher->knnMatch(cur_query->desc, cur_train->desc, in, 2); //knn matcher
        for( int i = 0 ; i < in.size(); i++) {
            if(in[i][0].distance < ratio_thresh * in[i][1].distance) {
                good.push_back(in[i][0]);
            }
        }

    } else if (p->submatch_type == BEST_MATCH || p->submatch_type == SPLIT_MATCH) {
        matcher->match(cur_query->desc, cur_train->desc, good); //normal mathcer
    }

    sort(good.begin(), good.end());
    dl.Logger("First matche size %d ", good.size());

    matches = RefineMatch(good);
    
    if (matches.size() < 5 ) {
        return -1; 
    }

    vector<Point2f> train_pt, query_pt;
    vector<Point2f> scaled_train_pt, scaled_query_pt;

    for (vector<DMatch>::const_iterator it = matches.begin(); it != matches.end(); it++)
    {
        float tx = cur_train->ip[it->trainIdx].pt.x;
        float ty = cur_train->ip[it->trainIdx].pt.y;

        float qx = cur_query->ip[it->queryIdx].pt.x;
        float qy = cur_query->ip[it->queryIdx].pt.y;

        //        dl.Logger("_pt push %f %f %f %f ", tx, ty, qx, qy);

        if ((tx > 0 && ty > 0) && (tx < cur_train->img.cols && ty < cur_train->img.rows) &&
            (qx > 0 && qy > 0) && (qx < cur_query->img.cols && qy < cur_query->img.rows))
        {
            if (p->p_scale != 1)
            {
                scaled_train_pt.push_back(Point2f(tx, ty));
                scaled_query_pt.push_back(Point2f(qx, qy));
                train_pt.push_back(Point2f(tx * p->p_scale, ty * p->p_scale));
                query_pt.push_back(Point2f(qx * p->p_scale, qy * p->p_scale));
            }
            else
            {
                train_pt.push_back(Point2f(tx, ty));
                query_pt.push_back(Point2f(qx, qy));
            }
        }
    }

    if(p->submatch_type == KNN_MATCH  || p->submatch_type == BEST_MATCH) {
        Mat _h = findHomography(train_pt, query_pt, FM_RANSAC);
        //Mat _h = getAffineTransform(query_pt, train_pt);
        //Mat _h = estimateAffine2D(train_pt, query_pt);
        //Mat _h = estimateRigidTransform(query_pt, train_pt, false);
        cur_query->matrix_fromimg = _h;

        if (p->p_scale != 1)
        {
            Mat sc_h = findHomography(scaled_train_pt, scaled_query_pt, FM_RANSAC);
            cur_query->matrix_scaledfromimg = sc_h;
        }
#if defined _DEBUG
        static int fi = 0;
        Mat outputImg = cur_train->img.clone();
        drawMatches(cur_query->img, cur_query->ip, cur_train->img, cur_train->ip,
                    matches, outputImg, Scalar::all(-1), Scalar::all(-1), vector<char>(), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
        char filename[30] = {
            0,
        };
        sprintf(filename, "saved/%2d_match.png", fi);
        imwrite(filename, outputImg);
        fi++;
        
#endif        
    } else if (p->submatch_type == SPLIT_MATCH) {

#if defined _DEBUG
        static int fi2 = 0;
        Mat outputImg = cur_train->img.clone();
        drawMatches(cur_query->img, cur_query->ip, cur_train->img, cur_train->ip,
                    matches, outputImg, Scalar::all(-1), Scalar::all(-1), vector<char>(), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
        char filename[40] = {
            0,
        };
        sprintf(filename, "saved/%2d_match_splitbefore.png", fi2);
        imwrite(filename, outputImg);
        fi2++;
        
#endif        
        MatchSplit(train_pt, query_pt);
    }


    return matches.size();
}

int Covenant::MatchSplit(vector<Point2f> m_train, vector<Point2f>m_query) {

    dl.Logger(" Split Match Start! m_train size %d %d   ", m_train.size(), m_query.size());
    vector<Point2f> mr_train[4];
    vector<Point2f> mr_query[4];
    int count[4] = { 0, };    

    float range = (float)p->circle_fixedpt_radius * 1.2;

    for(int i = 0 ; i < m_train.size(); i ++) {
        int index = -1;
        for(int j = 0 ; j < 4 ; j ++) {
//            dl.Logger("value check %f %f ", abs(cur_train->four_fpt[j].x - m_train[i].x), 
//                abs(cur_train->four_fpt[j].y - m_train[i].y));

            if ( abs(cur_train->four_fpt[j].x - m_train[i].x) < range &&
                    abs(cur_train->four_fpt[j].y - m_train[i].y) < range) {
                    index = j;
                    break;
            }
        }

        if(index >= 0) {
            mr_train[index].push_back(Point2f(m_train[i].x, m_train[i].y));
            mr_query[index].push_back(Point2f(m_query[i].x, m_query[i].y));            
            count[index]++;            
//            dl.Logger("train %f %f query %f %f insert to index %d (%f %f) count %d", m_train[i].x, m_train[i].y,
//            m_query[i].x, m_query[i].y, index, cur_train->four_fpt[index].x, cur_train->four_fpt[index].y, //count[index]); 

        } /*else {
            dl.Logger("Didn't belong to any range..%f %f ", m_train[i].x, m_train[i].y);
        } */
    }
    dl.Logger("region match count %d %d %d %d ", count[0], count[1], count[2], count[3]);

    int max_index = -1;
    int max_val = 0;
    for(int i = 0; i < 4; i++) {
        if(count[i] > max_val) {
            max_index = i;
            max_val = count[i];
        }
    }

    if(max_val < 6) {
        dl.Logger("Can't find any available answer..");
        return -1;
    } else 
        dl.Logger("max value index %d, score %d", max_index, max_val);


    for(int k = 0 ; k < 4 ; k ++) {
        Mat _h;
        if(count[k] > 6) {
            //Mat _h = estimateRigidTransform(mr_train[k], mr_query[k], true);
            _h = estimateAffine2D(mr_train[k], mr_query[k]);            
            FPt newpt = mtrx.TransformPtbyAffine(cur_train->four_fpt[k], _h);
            dl.Logger("Split match point move[%d] %f %f -> %f %f ", k, cur_train->four_fpt[k].x, 
                    cur_train->four_fpt[k].y, newpt.x, newpt.y);
            cur_query->four_fpt[k].x = newpt.x;
            cur_query->four_fpt[k].y = newpt.y;        
        }
        else {
            dl.Logger("Point is not enough..- max_index value apply ");
            //Mat _h = estimateRigidTransform(mr_train[max_index], mr_query[max_index], true);                
            _h = estimateAffine2D(mr_train[max_index], mr_query[max_index]);
            FPt newpt = mtrx.TransformPtbyAffine(cur_train->four_fpt[k], _h);
            dl.Logger("Split match point move[%d] %f %f -> %f %f ", k, cur_train->four_fpt[k].x, 
                    cur_train->four_fpt[k].y, newpt.x, newpt.y);
            cur_query->four_fpt[k].x = newpt.x;
            cur_query->four_fpt[k].y = newpt.y;        

        }

//        float ncc_ = ncc(k, _h);
    }

    //Mat _h = estimateRigidTransform(mr_train[max_index], mr_query[max_index], true);
    /*
    Mat _h = estimateAffine2D(mr_train[max_index], mr_query[max_index]);    
    FPt newpt = mtrx.TransformPtbyAffine(cur_train->four_fpt[max_index], _h);
    dl.Logger("Split match point move[%d] %f %f -> %f %f ", max_index, cur_train->four_fpt[max_index].x, 
            cur_train->four_fpt[max_index].y, newpt.x, newpt.y);
    cur_query->four_fpt[max_index].x = newpt.x;
    cur_query->four_fpt[max_index].y = newpt.y;        

    if( ncc_ >= 0.7) {
        dl.Logger("Homography OK. apply another point.. ");
        for(int k = 0 ; k < 4 ; k ++) {
            if(k != max_index) {
                FPt newpt = mtrx.TransformPtbyAffine(cur_train->four_fpt[k], _h);
                cur_query->four_fpt[k].x = newpt.x;
                cur_query->four_fpt[k].y = newpt.y;        
            }
        }
    }
    else {
        dl.Logger("Homography fail. Try to searching again");
    }
 */
    return ERR_NONE;
}

vector<DMatch> Covenant::RefineMatch(vector<DMatch> good) {

    vector<DMatch>matches;
    vector<DMatch>last;

    int *t_hist = (int *)g_os_malloc(sizeof(int) * cur_train->ip.size());
    int *q_hist = (int *)g_os_malloc(sizeof(int) * cur_query->ip.size());
    for (int a = 0; a < cur_train->ip.size(); a++)
        t_hist[a] = 0;
    for (int b = 0; b < cur_query->ip.size(); b++)
        q_hist[b] = 0;

    for (int t = 0; t < good.size(); t++)
    {
        if (t_hist[good[t].trainIdx] == 0 && q_hist[good[t].queryIdx] == 0)
        {
            matches.push_back(good[t]);
            t_hist[good[t].trainIdx]++;
            q_hist[good[t].queryIdx]++;
        }
//        else
//            dl.Logger("double check %d %d ", good[t].trainIdx, good[t].queryIdx);
    }

    g_os_free(t_hist);
    g_os_free(q_hist);

    last = RemoveOutlier(matches);
    if(p->submatch_type != SPLIT_MATCH) {

        if (last.size() > 100) {
            while (last.size() >= 100) {
                last.pop_back();
            }
            dl.Logger("matches->pop_back size %d ", last.size());
        }
    }
    return last;

};

vector<DMatch> Covenant::RemoveOutlier(vector<DMatch> matches) {

    vector<DMatch> result;
    dl.Logger("Remove Outlier is called %d ", matches.size());
    double covar_deg = 0;
    double covar_dist = 0;
    double t_covar_deg = 0;        
    double t_covar_dist = 0;
    double degsum = 0;
    double dist_sum = 0;
    double degavg = 0;
    double distavg = 0;

    for( vector<DMatch>::const_iterator it = matches.begin(); it != matches.end(); it++) {
        float tx = cur_train->ip[it->trainIdx].pt.x;
        float ty = cur_train->ip[it->trainIdx].pt.y;
        float qx = cur_query->ip[it->queryIdx].pt.x;
        float qy = cur_query->ip[it->queryIdx].pt.y;
        float dx = qx - tx;
        float dy = qy - ty;

        dist_sum += sqrt( dx * dx + dy * dy );
        float orideg = fastAtan2( dx, dy);
        if (orideg > 180 )
            orideg -= 180;
        degsum += orideg;
//        dl.Logger(" diff %f %f ", sqrt( dx * dx  + dy * dy ), orideg);
    }


    distavg = dist_sum / (float)matches.size();
    degavg = degsum / (float)matches.size();

    for(vector<DMatch>::const_iterator it = matches.begin(); it != matches.end(); it++) {
        float tx = cur_train->ip[it->trainIdx].pt.x;
        float ty = cur_train->ip[it->trainIdx].pt.y;
        float qx = cur_query->ip[it->queryIdx].pt.x;
        float qy = cur_query->ip[it->queryIdx].pt.y;
        float dx = qx - tx;
        float dy = qy - ty;

        t_covar_dist += (sqrt( dx * dx + dy * dy)  - distavg) * (sqrt( dx * dx + dy * dy ) - distavg);

        float orideg = fastAtan2( dx, dy);
        if (orideg > 180 )
            orideg -= 180;

        t_covar_deg += (orideg - degavg)*(orideg - degavg);
    }

    covar_dist = sqrt(t_covar_dist/(float)matches.size());
    covar_deg = sqrt(t_covar_deg/(float)matches.size());    

    dl.Logger("covar %f %f threshold %f %f  ", covar_deg, covar_dist);

    for(vector<DMatch>::const_iterator it = matches.begin(); it != matches.end(); it++) {
        float tx = cur_train->ip[it->trainIdx].pt.x;
        float ty = cur_train->ip[it->trainIdx].pt.y;
        float qx = cur_query->ip[it->queryIdx].pt.x;
        float qy = cur_query->ip[it->queryIdx].pt.y;
        float dx = qx - tx;
        float dy = qy - ty;

        float mh_distance_dist = abs(sqrt(dx * dx + dy * dy) - distavg )/covar_dist;
        float orideg = fastAtan2( dx, dy);
        if (orideg > 180 )
            orideg -= 180;

        float mh_distance_deg = abs(orideg - degavg)/covar_deg;

        //dl.Logger(" mh distance %f -> %f  angle %f -> %f ", sqrt(dx * dx + dy * dy), mh_distance_dist, orideg, mh_distance_deg);

        if( mh_distance_deg >= 1.0 || mh_distance_dist >= 1.0)
            //dl.Logger("mh distance is over limit %f %f ", mh_distance_deg, mh_distance_dist);
            continue;
        else
           result.push_back(*it);
    }

    return result;
}

int Covenant::PostProcess() {

    if(p->match_type == PYRAMID_MATCH) {
        for(int i = 0 ; i < p->roi_count ; i ++)
            dl.Logger("DIFF %2.3f %2.3f ", cur_train->four_fpt[i].x - cur_query->four_fpt[i].x,
                cur_train->four_fpt[i].y - cur_query->four_fpt[i].y);
        return ERR_NONE;
    }

    if (cur_query->id == 0) {
        FindBaseCoordfromWd(NORMAL_VECTOR_CAL);
        return ERR_NONE;
    }
    dl.Logger(" Post Process start.. ");

    if (p->submatch_type == SPLIT_MATCH) {
        FindBaseCoordfromWd(NORMAL_VECTOR_CAL);    
        return ERR_NONE;
    }

    float err = 0;
    //move centerpoint
    FPt newcen = mtrx.TransformPtbyHomography(cur_train->center, cur_query->matrix_fromimg);
    dl.Logger("Query center answer(%f, %f) - (%f, %f)", cur_query->center.x, cur_query->center.y, newcen.x, newcen.y);
    cur_query->center = newcen;

    //move 4point
    for (int i = 0; i < p->roi_count; i++)
    {
        FPt newpt = mtrx.TransformPtbyHomography(cur_train->four_fpt[i], cur_query->matrix_fromimg);

        dl.Logger(" four pt move [%d] answer (%f, %f) - (%f, %f) ", i,
               cur_query->four_fpt[i].x, cur_query->four_fpt[i].y,
               newpt.x, newpt.y);
        cur_query->four_fpt[i].x = newpt.x;
        cur_query->four_fpt[i].y = newpt.y;
    }

    FindBaseCoordfromWd(NORMAL_VECTOR_CAL);

    //move user circle input
    if (p->roi_type == CIRCLE && p->masking_type == USER_INPUT_CIRCLE)
    {
        Mat apply_homo;
        if (p->p_scale != 1)
            apply_homo = cur_query->matrix_scaledfromimg;
        else
            apply_homo = cur_query->matrix_fromimg;

        for (int i = 0; i < p->roi_count; i++)  {
            Pt newpt = mtrx.TransformPtbyHomography(&p->circles[i].center, apply_homo);
            p->circles[i].center = newpt;
        }
    }
    return ERR_NONE;
}

int Covenant::FindBaseCoordfromWd(int mode)
{
    dl.Logger("FindBaseCoordfromW start ");
    Mat cm(3, 3, CV_32F, p->camera_matrix);
    Mat sc(4, 1, CV_32F, p->skew_coeff);
    Mat ppset1(4, 3, CV_32F);
    Mat ppset2(4, 2, CV_32F);
    for (int i = 0; i < 4; i++)
    {
        ppset1.at<float>(i, 0) = p->world->four_fpt[i].x;
        ppset1.at<float>(i, 1) = p->world->four_fpt[i].y;
        ppset1.at<float>(i, 2) = p->world->four_fpt[i].z;

        ppset2.at<float>(i, 0) = (float)cur_query->four_fpt[i].x;
        ppset2.at<float>(i, 1) = (float)cur_query->four_fpt[i].y;

    }

    Mat ret1(3, 1, CV_32F);
    Mat ret2(3, 1, CV_32F);

    bool result = solvePnP(ppset1, ppset2,
                        cm, sc,
                        ret1, ret2,
                        false, SOLVEPNP_ITERATIVE);

    cur_query->rot_matrix = ret1;
    cur_query->trans_matrix = ret2;

    dl.Logger("solve pnp is done result %d :  %d %d %d %d ", result,
        ret1.cols, ret1.rows, ret2.cols, ret2.rows);
    for (int i = 0; i < cur_query->rot_matrix.rows; i++)
        for (int j = 0; j < cur_query->rot_matrix.cols; j++)
            dl.Logger("[%d][%d] %f ", i, j, cur_query->rot_matrix.at<float>(i, j));

    dl.Logger(" ---- ");
    for (int i = 0; i < cur_query->trans_matrix.rows; i++)
        for (int j = 0; j < cur_query->trans_matrix.cols; j++)
            dl.Logger("[%d][%d] %f ", i, j, cur_query->trans_matrix.at<float>(i, j));

    Mat projectedNormal;
    float v_normal[6][3];
    float loca_x = 50.0;
    float loca_y = 50.0;        
    float nlen = 10;
    v_normal[0][0] = loca_x;
    v_normal[0][1] = loca_y;
    v_normal[0][2] = 0.0;

    v_normal[1][0] = loca_x + nlen;
    v_normal[1][1] = loca_y;
    v_normal[1][2] = 0;

    v_normal[2][0] = loca_x;
    v_normal[2][1] = loca_y + nlen;
    v_normal[2][2] = 0;

    v_normal[3][0] = loca_x;
    v_normal[3][1] = loca_y;
    v_normal[3][2] = -nlen * 3;

    v_normal[4][0] = loca_x - nlen;
    v_normal[4][1] = loca_y;
    v_normal[4][2] = 0;

    v_normal[5][0] = loca_x;
    v_normal[5][1] = loca_y - nlen;
    v_normal[5][2] = 0;

    Mat tvec(6, 3, CV_32F, v_normal);
    projectPoints(tvec, cur_query->rot_matrix, cur_query->trans_matrix,
        cm, sc, projectedNormal);

    cur_query->projected_normal = projectedNormal.clone();

    Point2f tp1 = Point2f(cur_query->projected_normal.at<Point2f>(0).x, 
                cur_query->projected_normal.at<Point2f>(0).y);

    Point2f tp2 = Point2f(cur_query->projected_normal.at<Point2f>(3).x, 
            cur_query->projected_normal.at<Point2f>(3).y);

    double distance = norm(tp2 - tp1);
    double degree = fastAtan2(tp1.y - tp2.y, tp1.x - tp2.x);
    degree = 90 - degree; // skew ratio of vector basd 90deg normal vector
    degree = 360 - degree; // inverse rotation for image
    dl.Logger(" tp1 %f %f tp2 %f %f dx %f dy %f ", tp1.x, tp1.y, tp2.x, tp2.y, tp1.x-tp2.x, tp1.y - tp2.y);
    cur_query->rod_norm = distance;
    cur_query->rod_degree = degree;
    double scale = cur_train->rod_norm / cur_query->rod_norm;    
    dl.Logger("normal vector norm %f degree %f scale %f ", distance, degree, scale);


    if (cur_query->id == 0)
        cur_query->rod_rotation_matrix = getRotationMatrix2D(Point2f(cur_query->center.x, cur_query->center.y), degree, 1);
    else
        cur_query->rod_rotation_matrix = getRotationMatrix2D(Point2f(cur_query->center.x, cur_query->center.y), degree, scale);

    return 1;
}


ADJST Covenant::CalAdjustData()
{
    ADJST newadj;

    double interval = cur_query->rod_norm - cur_train->rod_norm;
    /*     double agvx = (cur_train->center.x + cur_query->center.x)/2;
    double agvy = (cur_train->center.y + cur_query->center.y)/2;
 */
    double angle = -cur_query->rod_degree;
    if (angle >= -180)
        angle += -90;
    else
        angle += 270;

    double scale = cur_train->rod_norm / cur_query->rod_norm;
    double adjustx = cur_train->center.x - cur_query->center.x;
    double adjusty = cur_train->center.y - cur_query->center.y;
    double rotatex = cur_query->center.x;
    double rotatey = cur_query->center.y;

    double angleadjust = -1 * (angle + 90);
    double radian = angleadjust * M_PI / 180;
    double width = cur_query->ori_img.cols;
    double height = cur_query->ori_img.rows;

    Point2f pt1, pt2, pt3, pt4, ptR1, ptR2, ptR3, ptR4;
    pt1.x = (float)(cur_query->center.x * (1 - scale));
    pt1.y = (float)(cur_query->center.y * (1 - scale));

    pt2.x = (float)(pt1.x + width * scale);
    pt2.y = pt1.y;

    pt3.x = pt2.x;
    pt3.y = (float)(pt1.y + height * scale);

    pt4.x = pt1.x;
    pt4.y = pt3.y;

    Point2f ptcenter(rotatex, rotatey);
    ptR1 = mtrx.GetRotatePoint(ptcenter, pt1, radian);
    ptR2 = mtrx.GetRotatePoint(ptcenter, pt2, radian);
    ptR3 = mtrx.GetRotatePoint(ptcenter, pt3, radian);
    ptR4 = mtrx.GetRotatePoint(ptcenter, pt4, radian);

    int margin_l = 0;
    int margin_t = 0;
    int margin_r = width;
    int margin_b = height;

    if (ptR1.x + adjustx > margin_l)
        margin_l = (int)(ptR1.x + adjustx);
    if (ptR1.y + adjusty > margin_t)
        margin_t = (int)(ptR1.y + adjusty);

    if (ptR2.x + adjustx < margin_r)
        margin_r = (int)(ptR2.x + adjustx);
    if (ptR2.y + adjusty > margin_t)
        margin_t = (int)(ptR2.y + adjusty);

    if (ptR3.x + adjustx < margin_r)
        margin_r = (int)(ptR3.x + adjustx);
    if (ptR3.y + adjusty < margin_b)
        margin_b = (int)(ptR3.y + adjusty);

    if (ptR4.x + adjustx > margin_l)
        margin_l = (int)(ptR4.x + adjustx);
    if (ptR4.y + adjusty < margin_b)
        margin_b = (int)(ptR4.y + adjusty);

    if (margin_l > margin_t * width / height)
        margin_t = margin_l * height / width;
    else
        margin_l = margin_t * width / height;

    if (margin_r < margin_b * width / height)
        margin_b = margin_r * height / width;
    else
        margin_r = margin_b * width / height;

    // Margin 뒤집힘 현상 임시 방지
    if (margin_l > margin_r)
    {
        int nTemp = margin_l;
        margin_l = margin_r;
        margin_r = nTemp;
    }
    if (margin_t > margin_b)
    {
        int nTemp = margin_t;
        margin_t = margin_b;
        margin_b = nTemp;
    }

    newadj.angle = angle;
    newadj.rotate_centerx = rotatex;
    newadj.rotate_centery = rotatey;
    newadj.scale = scale;
    newadj.trans_x = adjustx;
    newadj.trans_x = adjusty;
    newadj.rect = Rect(margin_l, margin_t, margin_r - margin_l, margin_b - margin_t);

    dl.Logger("Adjust data ..");
    dl.Logger("angle  %f  centerx %f  centery %f scale %f ", newadj.angle, newadj.rotate_centerx, newadj.rotate_centery,
           newadj.scale);
    dl.Logger(" trans x %f y %f rect %f %f %f %f ", newadj.trans_x, newadj.trans_y,
           newadj.rect.x, newadj.rect.y, newadj.rect.width, newadj.rect.height);

    return newadj;
}

int Covenant::WarpingStep1()
{
    //Test
    //Image Warping
    vector<Point2f> t_pset;
    vector<Point2f> q_pset;

    for (int i = 0; i < 4; i++)
    {
        t_pset.push_back(Point2f(cur_train->four_fpt[i].x, cur_train->four_fpt[i].y));
        q_pset.push_back(Point2f(cur_query->four_fpt[i].x, cur_query->four_fpt[i].y));

        dl.Logger("t_pset %f %f  -- t_qset %f %f", t_pset[i].x, t_pset[i].y,
               q_pset[i].x, q_pset[i].y);
    }
    dl.Logger(" Start find homography ");

    Mat _h = findHomography(q_pset, t_pset, 0);
    dl.Logger(" Get homography ");

    for (int i = 0; i < _h.rows; i++)
        for (int j = 0; j < _h.cols; j++)
            dl.Logger("[%d][%d] %lf ", i, j, _h.at<double>(i, j));

    dl.Logger("ori image size %d %d ", cur_query->ori_img.cols, cur_query->ori_img.rows);

    Mat final;
    warpPerspective(cur_query->ori_img, final, _h, Size(cur_query->ori_img.cols, cur_query->ori_img.rows));

    static int index = 0;
    char filename[30] = {
        0,
    };
    sprintf(filename, "saved/%2d_perspective.png", index);
    imwrite(filename, final);

    Mat mcenter(3, 1, CV_64F);
    mcenter.at<double>(0) = cur_train->center.x;
    mcenter.at<double>(1) = cur_train->center.y;
    mcenter.at<double>(2) = 1;

    Mat mresult = _h * mcenter;

    dl.Logger("estimated center point %f %f ", mresult.at<double>(0), mresult.at<double>(1));
    dl.Logger("error  %f %f ", cur_query->center.x - mresult.at<double>(0), cur_query->center.y - mresult.at<double>(1));
    //ROI Warping

    return ERR_NONE;
}

int Covenant::DecomposeHomography()
{

    Mat _h = cur_query->matrix_fromimg;
    Mat cm(3, 3, CV_32F, p->camera_matrix);
    dl.Logger("_h from img is cols %d rows %d ", _h.cols, _h.rows);

    vector<Mat> Rs_decomp, ts_decomp, normals_decomp;

    int solutions = decomposeHomographyMat(_h, cm, Rs_decomp, ts_decomp, normals_decomp);

    for (int i = 0; i < solutions; i++)
    {
        Mat rvec_decomp;
        Rodrigues(Rs_decomp[i], rvec_decomp);
        dl.Logger("Solution %d : ", i);
        Mat rvec_t = rvec_decomp.t();
        dl.Logger("rvec from homography decomposition: %f %f %f ", rvec_t.at<double>(0), rvec_t.at<double>(1), rvec_t.at<double>(2));
        Mat ts_decom_t = ts_decomp[i].t();
        Mat normal_decom_t = normals_decomp[i].t();
        dl.Logger("tvec from homography decomposition: %f %f %f ", ts_decom_t.at<double>(0), ts_decom_t.at<double>(1),
               ts_decom_t.at<double>(2));
        dl.Logger("plane normal from homography decomposition: %f %f %f ", normal_decom_t.at<double>(0),
               normal_decom_t.at<double>(1), normal_decom_t.at<double>(2));
    }

    return ERR_NONE;
}

int Covenant::FindHomographyP2P()
{

    Mat ppset1(4, 2, CV_32F);
    Mat ppset2(4, 2, CV_32F);
    for (int i = 0; i < 4; i++)
    {
        ppset1.at<float>(i, 0) = cur_train->four_fpt[i].x;
        ppset1.at<float>(i, 1) = cur_train->four_fpt[i].y;

        ppset2.at<float>(i, 0) = cur_query->four_fpt[i].x;
        ppset2.at<float>(i, 1) = cur_query->four_fpt[i].y;

        dl.Logger("ppset %f %f %f -- %f %f", ppset1.at<float>(i, 0), ppset1.at<float>(i, 1),
               ppset2.at<float>(i, 0), ppset2.at<float>(i, 1));
    }

    Mat h = findHomography(ppset1, ppset2);
    for (int i = 0; i < h.rows; i++)
        for (int j = 0; j < h.cols; j++)
            dl.Logger("[%d][%d] %lf ", i, j, h.at<double>(i, j));

    Mat mcenter(3, 1, CV_64F);
    mcenter.at<double>(0) = cur_train->center.x;
    mcenter.at<double>(1) = cur_train->center.y;
    mcenter.at<double>(2) = 1;
    Mat mret = h * mcenter;

    double newx = mret.at<double>(0) / mret.at<double>(2);
    double newy = mret.at<double>(1) / mret.at<double>(2);
    dl.Logger("transformed cetner by P2P : %f %f ", newx, newy);

    return ERR_NONE;
}