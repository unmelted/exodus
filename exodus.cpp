
/*****************************************************************************
*                                                                            *
*                            Genesis          								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : genesis.pp
    Author(S)       : Me Eunkyung
    Created         : 21 Sep 2021

    Description     : genesis.cpp
    Notes           : Python - C library connector.
*/

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include "src/DefData.hpp"
#include "src/Extractor.hpp"

using namespace std;

#define VER "0.1.0"

/* int TestFeature(unsigned char* framedata)
{
    printf("Enter! \n");
    static int index = 0;    
    printf("Enter! 2\n");
    IMG* bframe = CreateImage(1024, 688, 0, framedata);
    //IMG* bframe = CreateImage(1024, 688, 0);    
    printf(" 2 \n");
    char fname[25] = {0, };
    sprintf(fname, "test/saveimage_%d.png", index);
    index++;
    printf(" 3 %s \n", fname);    
    SaveImagePNG(bframe, fname);
    printf(" 4 \n");
    DestroyImage(bframe);
}
 */
void Process(int cnt, int* region, char* img_path);
void Finish();

extern "C" {
/*     void Feature(unsigned char* buffers) {
        TestFeature(buffers);
    }
 */    
    void GetVerion() {
        cout<< "Cur Version : " << VER << endl;        
    }

    int Extract(int* buffers, char* img_path) {
        int cnt = buffers[0];

        Logger( "received count %d path %s", cnt, img_path);
        Process(cnt, buffers, img_path);
        return 1;
    }

    void Exit() {
        Finish();
    }
}

Extractor* ext;
void Process(int cnt, int* region, char* img_path) {
        string imgset(img_path);    
        ext = new Extractor(imgset, cnt, region);
        ext->Execute();

        //ext->DrawInfo();
        
}

void Finish() {    
    ext->~Extractor();
}