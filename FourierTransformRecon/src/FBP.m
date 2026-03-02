clear all;
close all;

%% 获得投影数据
p_N = 256;                  % 图像默认大小256
theta_N = 180;              % 180度平行投影

P = phantom(p_N);           % 256*256的头骨幻影数据
theta = 0:(theta_N-1);      % 0~180度平行投影
[R,xp] = radon(P,theta);    % radon变换后，R就是投影数据，即正弦图

%% 求取直接反投影后的断层图像，以及滤波反投影后的图像
[I H] = iradon(R,theta,'linear','none');         % iradon变换，不滤波，则I就是直接反投影后的断层图像
[I1 H1] = iradon(R,theta,'linear','Ram-Lak');    % iradon变换，斜坡滤波，则I1就是滤波反投影后的图像

figure, 
subplot(1,3,1), imshow(P,[]),  title('原始图像');
subplot(1,3,2), imshow(I,[]),  title('不滤波直接反投影后的图像');
subplot(1,3,3), imshow(I1,[]), title('先斜坡滤波后反投影重建的图像');
