clear all;
close all;

%% 获得原始图像、投影数据、以及直接反投影后的断层图像
p_N = 256;                                  % 图像默认大小256
theta_N = 180;                              % 180度平行投影
pad_N = 1024;                               % 计算二维卷积时对数据扩展，变为1024

P = phantom(p_N);                           % 256*256的头骨幻影数据
theta = 0:(theta_N-1);                      % 0~180度平行投影
[R,xp] = radon(P,theta);                    % radon变换后，R就是投影数据，即正弦图
[I H] = iradon(R,theta,'linear','none');    % iradon变换，不滤波，则I就是直接反投影后的断层图像

%% 获得PSF点扩展函数 h(x,y), h(x,y) = 1/sqrt(x^2+y^2)
N = pad_N;
xP = zeros(N,N);
yP = zeros(N,N);
hP = zeros(N,N);
for i = 1:1:N
    for j = 1:1:N
        xP(i,j) = ((-(N-1)/2)+(j-1)*1)*2*pi/N;
        yP(i,j) = ((-(N-1)/2)+(i-1)*1)*2*pi/N;
        hP(i,j) = sqrt(xP(i,j)^2+yP(i,j)^2);
        if(hP(i,j)>0)
            hP(i,j) = 1/hP(i,j);
        end
    end
end

%% 计算二维卷积，Pc = p(x,y) ** h(x,y)
row_pad = N;
row_P = size(P,1);
row_start = (row_pad-row_P)/2;
P_pad = zeros(row_pad,row_pad);                     % 扩充成1024*1024
for i = 1:1:row_P
    for j = 1:1:row_P
        P_pad(i+row_start,j+row_start) = P(i,j);    % 里面256*256是原始头骨幻影数据，外面的都是0
    end
end

Pc = conv2(hP,P_pad,'same');                        % 二维离散卷积，结果Pc也是1024*1024
crop = row_start;
Pc = Pc(crop+1:end-crop,crop+1:end-crop);           % 原始数据256*256，因此需要对数据截断

figure, 
subplot(1,3,1), imshow(P,[]), title('原始图像');
subplot(1,3,2), imshow(Pc,[]), title('原始图像与PSF二维卷积后的图像');
subplot(1,3,3), imshow(I,[]), title('直接反投影后的断层图像');






