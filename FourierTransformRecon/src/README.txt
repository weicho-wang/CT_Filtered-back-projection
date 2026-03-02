CT二维图像重建的几种方法
代码在MATLAB R2010b调试成功，有问题联系 2014550763@qq.com

DFR.m   直接傅里叶反变换法

BPF.m 	先反投影后滤波
BPF_conv2.m	原始图像与PSF（点扩展函数）的二维卷积比较

FBP.m	滤波反投影

CBP	卷积反投影
    CBP.m 	卷积反投影主程序
    DirectBackProject.m		累加法反投影的MATLAB代码
    Backproject.c, Backproject.mexw64	mex 编译的C代码，供CBP.m调用。
