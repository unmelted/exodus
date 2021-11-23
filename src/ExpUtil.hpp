  
/*****************************************************************************
*                                                                            *
*                            Util         								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : Util.hpp
    Author(S)       : Me Eunkyung
    Created         : 13 Oct 2021

    Description     : Util.hpp
    Notes           : Uility - logging, json export
*/


#include <cmath>
#include "DefData.hpp"
#include "common/json.hpp"

class ExpUtil {

public:
    void Export(vector<string>image_path, vector<SCENE>cal_group, PARAM* p);
    void ExportforApp(vector<string>image_path, vector<SCENE>cal_group, PARAM* p);
    
private:

};