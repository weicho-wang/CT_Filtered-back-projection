clear all;
close all;

%% 获得原始图像、投影数据、以及滤波反投影后的图像
p_N = 256;  
theta_N = 180;  
pad_N = 1024;        
       
P = phantom(p_N);                                   % 256*256的头骨幻影数据
theta = 0:(theta_N-1);                              % 0~180度平行投影
[R,xp] = radon(P,theta);                            % radon变换后，R就是投影数据，即正弦图
[I H] = iradon(R,theta,'linear','Ram-Lak');         % iradon变换，斜坡滤波，则I就是滤波反投影后的图像

%% 求取Ram-Lak斜坡滤波器的时域表达式 hRL
N = pad_N;  d = 1.0;                                % d是滤波因子，0~1之间
N_div2 = round(N/2);    
hRL = zeros(N,1);
hRL(N_div2+1) = 0.25*d*d;                           % 卷积核的中心点最大，为 0.25d^2
for i = 1:1:(N_div2-1)
    if(mod(i,2)==0)
        hRL(i+N_div2+1) = 0;                        % 偶数坐标为0
    else
        hRL(i+N_div2+1) = -1.0/(i*i*d*d*pi*pi);     % 奇数坐标的值，逐渐减小
    end
end

for i = 2:1:N_div2
    hRL(i) = hRL(N-i+2);                            % 左右对称
end

%% 对投影数据做填充，然后每个角度下的投影数据与hRL做一维卷积
proj_sino = R;  [rowR colR] = size(R);
if mod(length(proj_sino(:,1)),2)==1                       % 如果是奇数个接收栅格，则需要扩展为偶数，填充0
   proj_sino = [proj_sino;zeros(1,size(proj_sino,2))];    % 填充一行，367*180变为了368*180
end
pad_row = (pad_N-size(proj_sino,1))/2;                    % （1024-368）/2 = 328    
proj_sino = padarray( proj_sino,[pad_row 0],0,'both');    % 368*180，上下填充328行，变为了1024*180
row_start = ((pad_N-rowR)/2);   
row_start = floor(row_start);

filter_proj_sino = zeros(rowR,colR);
for j = 1:1:theta_N                                       % 每个投影角度下的数据，都要做卷积
    x = proj_sino(:,j);                                   % 一列数据，就是一个投影角度下的投影值
    y = conv(x,hRL,'same');                               % 一维卷积，获取卷积结果的中间部分
    for i = 1:1:rowR
        filter_proj_sino(i,j) = y(i+row_start);           % 丢弃填充的数据
    end
end

%% 卷积后的数据反投影
theta = pi*theta/180;   
costheta = cos(theta);  
sintheta = sin(theta);
imgC = Backproject(filter_proj_sino, costheta, sintheta, p_N, 1);                     % 反投影函数，调用C实现
imgM = DirectBackproject(filter_proj_sino, costheta, sintheta, p_N, 'linear');        % 反投影函数，直接调用MATLAB
imgC = imgC*pi/(2*length(theta));   
imgM = imgM*pi/(2*length(theta));

figure, 
subplot(1,4,1), imshow(P,[]), title('原始图像');
subplot(1,4,2), imshow(I,[]), title('先滤波后反投影重建的图像');
subplot(1,4,3), imshow(imgC,[]), title('先卷积后反投影重建的图像，调用C函数');
subplot(1,4,4), imshow(imgM,[]), title('先卷积后反投影重建的图像，直接MATLAB');





