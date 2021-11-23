  
/*****************************************************************************
*                                                                            *
*                            Util         								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : Util.cpp
    Author(S)       : Me Eunkyung
    Created         : 13 Oct 2021

    Description     : Util.Cpp
    Notes           : Uility - logging, json export
*/

#include "ExpUtil.hpp"
using json = nlohmann::json;

void ExpUtil::Export(vector<string>image_paths, vector<SCENE>cal_group, PARAM* p) {

    //world coord 
    //world point1 x = p->world->four_fpt[0].x;
    //world point1 y = p->world->four_fpt[0].x;
    //world center x = p->world->center.x;
    //world center y = p->world->center.y;    

    for (vector<SCENE>::const_iterator it = cal_group.begin(); it != cal_group.end(); it++)
    {
        //id = it-> image_paths[it->id] 
        //point1 x = it->four_fpt[0].x
        //point1 y = it->four_fpt[0].y
        //center x = it->center.x;
        //center y = it->center.y;
    }

    //file write. path = "/saved/xxx.json"

    json jObj = json::object();
    //World
    json world;
    world["X1"] = p->world->four_fpt[0].x;
    world["Y1"] = p->world->four_fpt[0].y;
    world["X2"] = p->world->four_fpt[1].x;
    world["Y2"] = p->world->four_fpt[1].y;
    world["X3"] = p->world->four_fpt[2].x;
    world["Y3"] = p->world->four_fpt[2].y;
    world["X4"] = p->world->four_fpt[3].x;
    world["Y4"] = p->world->four_fpt[3].y;
   
    //jObj["stadium"] = GROUNDTYPE.get();
    jObj["stadium"] = "SoccerHalf";
    jObj["world_coords"] = world;

    //2dPoint
    json point2d = json::object();
    point2d["UpperPosX"] = -1.0;
    point2d["UpperPosY"] = -1.0;
    point2d["MiddlePosX"] = -1.0;
    point2d["MiddlePosX"] = -1.0;
    point2d["LowerPosX"] = -1.0;
    point2d["LowerPosX"] = -1.0;
    //swipe
    json swipe;
    swipe["X1"] = -1.0;
    swipe["Y1"] = -1.0;
    swipe["X2"] = -1.0;
    swipe["Y2"] = -1.0;

    auto arr = json::array();

    for (vector<SCENE>::const_iterator it = cal_group.begin(); it != cal_group.end(); it++)
    {
  
        json jDsc = json::object();
  
        jDsc["dsc_id"] = image_paths[it->id];
        jDsc["pts_2d"] = point2d;

        //3dPoint
        json point3d = json::object();
        point3d["X1"] = it->four_fpt[0].x;
        point3d["Y1"] = it->four_fpt[0].y;
        point3d["X2"] = it->four_fpt[1].x;
        point3d["Y2"] = it->four_fpt[1].y;
        point3d["X3"] = it->four_fpt[2].x;
        point3d["Y3"] = it->four_fpt[2].y;
        point3d["X4"] = it->four_fpt[3].x;
        point3d["Y4"] = it->four_fpt[3].y;
        point3d["CenterX"] = it->center.x;
        point3d["CenterY"] = it->center.y;
   
        jDsc["pts_3d"] = point3d;
        jDsc["pts_swipe"] = swipe;


        arr.push_back(jDsc);

    }
    jObj["points"] = arr;

    //file write
    std::ofstream file("saved/UserPointData_" + getCurrentDateTime("date")  +"_" + getCurrentDateTime("now")+ ".pts");

    file << std::setw(4) << jObj << '\n';
    //std::cout << std::setw(4) << jObj << '\n';
}

void ExpUtil::ExportforApp(vector<string>image_paths, vector<SCENE>cal_group, PARAM* p) {


    json jObj = json::object();
    //World
    json world; json world1;
    world["X1"] = p->world->four_fpt[0].x;
    world["Y1"] = p->world->four_fpt[0].y;
    world["X2"] = p->world->four_fpt[1].x;
    world["Y2"] = p->world->four_fpt[1].y;
    world["X3"] = p->world->four_fpt[2].x;
    world["Y3"] = p->world->four_fpt[2].y;
    world["X4"] = p->world->four_fpt[3].x;
    world["Y4"] = p->world->four_fpt[3].y;

    
    world1["group"] = "Group1";                   //group
    world1["stadium"] = "SoccerHalf";             //stadium
    world1["world_coords"] = world;

    //2dPoint
    json point2d = json::object();
    point2d["UpperPosX"] = -1.0;
    point2d["UpperPosY"] = -1.0;
    point2d["MiddlePosX"] = -1.0;
    point2d["MiddlePosX"] = -1.0;
    point2d["LowerPosX"] = -1.0;
    point2d["LowerPosX"] = -1.0;

    json pt2d = json::object();
    pt2d["IsEmpty"] = false;
    pt2d["X"] = -1.0;
    pt2d["Y"] = -1.0;

    point2d["Upper"] = pt2d;
    point2d["Middle"] = pt2d;
    point2d["Lower"] = pt2d;

    auto world_arr = json::array();
    world_arr.push_back(world1);
    jObj["worlds"] = world_arr;



    auto arr = json::array();

    for (vector<SCENE>::const_iterator it = cal_group.begin(); it != cal_group.end(); it++)
    {

        json jDsc = json::object();

        jDsc["dsc_id"] = image_paths[it->id];
        jDsc["flip"] = 0;
        jDsc["Group"] = "Group1";
        jDsc["Width"] = "w";
        jDsc["Height"] = "h";
        jDsc["infection_point"] = 0;
        jDsc["swipe_base_length"] = -1.0;
        jDsc["ManualOffsetY"] = 0;
        jDsc["FocalLength"] = 0.0;
        jDsc["pts_2d"] = point2d;

        //3dPoint
        json point3d = json::object();
        point3d["X1"] = it->four_fpt[0].x;
        point3d["Y1"] = it->four_fpt[0].y;
        point3d["X2"] = it->four_fpt[1].x;
        point3d["Y2"] = it->four_fpt[1].y;
        point3d["X3"] = it->four_fpt[2].x;
        point3d["Y3"] = it->four_fpt[2].y;
        point3d["X4"] = it->four_fpt[3].x;
        point3d["Y4"] = it->four_fpt[3].y;
        point3d["CenterX"] = it->center.x;
        point3d["CenterY"] = it->center.y;
        json pt3d_1;
        pt3d_1["IsEmpty"] = false;
        pt3d_1["X"] = it->four_fpt[0].x;
        pt3d_1["Y"] = it->four_fpt[0].y;
        point3d["Point1"] = pt3d_1;
        json pt3d_2;
        pt3d_2["IsEmpty"] = false;
        pt3d_2["X"] = it->four_fpt[1].x;
        pt3d_2["Y"] = it->four_fpt[1].y;
        point3d["Point2"] = pt3d_2;
        json pt3d_3;
        pt3d_3["IsEmpty"] = false;
        pt3d_3["X"] = it->four_fpt[2].x;
        pt3d_3["Y"] = it->four_fpt[2].y;
        point3d["Point3"] = pt3d_3;
        json pt3d_4;
        pt3d_4["IsEmpty"] = false;
        pt3d_4["X"] = it->four_fpt[3].x;
        pt3d_4["Y"] = it->four_fpt[3].y;
        point3d["Point4"] = pt3d_4;
        json cent;
        cent["IsEmpty"] = false;
        cent["X"] = it->center.x;
        cent["Y"] = it->center.y;
        point3d["Center"] = cent;

        jDsc["pts_3d"] = point3d;

        arr.push_back(jDsc);

    }
    jObj["points"] = arr;

    //file write
    std::ofstream file("saved/app_UserPointData_" + getCurrentDateTime("date") +"_" + getCurrentDateTime("now") + ".pts");

    file << std::setw(4) << jObj << '\n';
    //cout << std::setw(4) << jObj << '\n';
}
