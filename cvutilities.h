#ifndef _CVUTILITIES_H
#define _CVUTILITIES_H

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "imageprocess.h"

void cvMat2IMG_RGB(const cv::Mat& m, IMG_RGB& i);
void IMG_RGB2cvMat(const IMG_RGB& i, cv::Mat& m);

std::string type2str(int type);

#endif
