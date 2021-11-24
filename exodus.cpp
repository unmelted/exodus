
/*****************************************************************************
*                                                                            *
*                            Exodus          								 *
*                                                                            *
*   Copyright (C) 2021 By 4dreplay, Incoporated. All Rights Reserved.        *
******************************************************************************

    File Name       : exodus.cpp
    Author(S)       : Me Eunkyung
    Created         : 23 Nov 2021

    Description     : exodus.cpp
    Notes           : Python - C library connector in c side.
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
#include "src/Covenant.hpp"

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
void Start(int cnt, int* gline, int* iline);
int ProcessStep(int step, int* buffers);
void Finish();

extern "C" {
/*     void Feature(unsigned char* buffers) {
        TestFeature(buffers);
    }
 */    
    void GetVerion() {
        cout<< "Cur Version : " << VER << endl;        
    }

    int Extract(int* gr_line, int* im_line) {
        int gcnt = gr_line[0];
        int icnt = im_line[0];
        if( gcnt != icnt) {
            Logger( "something wrong .. %d %d", gcnt, icnt);
        }
        Logger( "received count %d %d ", gcnt, icnt);
        Start(gcnt, gr_line, im_line);
        return 1;
    }

    int Process(int step, int* buffers) {
        int result = ProcessStep(step, buffers);
        return result;
    }

    void Exit() {
        Finish();
    }
}

Covenant* cov;
void Start(int cnt, int* gline, int* iline)  {
        cov = new Covenant(3840, cnt, gline, iline);
        //cov->Execute();

        //ext->DrawInfo();
        
}
int ProcessStep(int step, int* buffers)  {

    Logger("recived step %d ", step);
    return 1;
}

void Finish() {    
    cov->~Covenant();
}