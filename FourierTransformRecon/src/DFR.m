clear all;
close all;

%% 获得投影数据
p_N = 256;                  % 图像默认大小256
theta_N = 180;              % 180度平行投影
pad_N = 1024;               % 投影后变为367*180,每列数据就是当前角度下的投影值。367<512，考虑到傅里叶变换时需要基2对齐，故扩展为1024*180

P = phantom(p_N);           % 256*256的头骨幻影数据
theta = 0:(theta_N-1);      % 0~180度平行投影
[R,xp] = radon(P,theta);    % radon变换后，R就是投影数据，即正弦图

%% 根据每个投影角度下的投影数据，求取其一维傅里叶变换
proj_sino = R;
if mod(length(proj_sino(:,1)),2)==1                       % 如果是奇数个接收栅格，则需要扩展为偶数，填充0
   proj_sino = [proj_sino;zeros(1,size(proj_sino,2))];    % 填充一行，367*180变为了368*180
end
pad_row = (pad_N-size(proj_sino,1))/2;                    % （1024-368）/2 = 328    
proj_sino = padarray( proj_sino,[pad_row 0],0,'both');    % 368*180，上下填充328行，变为了1024*180
L_pad = pad_row + ceil(((p_N.*sqrt(2)+2)-p_N)/2)+1;       % 原始图像大小为p_N，傅里叶逆变换后为pad_N，则需要截断数据

proj_sino = ifftshift(proj_sino,1);                       % 将实际数据移动到两端，中间的是填充的0
f_p = fft(proj_sino,[ ],1);                               % 求取正弦图的一维傅里叶变换，1024*180
f_p = fftshift(f_p,1);                                    % 反向处理，中间的数据移动到两端，两端的数据移动到中间

%% 极坐标栅格化为笛卡尔坐标
nfp = length(f_p(:,1));
omega_sino = (-(nfp-1)/2:(nfp-1)/2).*(2*pi/size(f_p,1));  % 极半径范围 -pi~pi，分成1024等份
theta = theta*pi/180;                                     % 角度转换为弧度，范围0~pi，180等份
[theta_grid, omega_grid] = meshgrid(theta,omega_sino);    % 网格后，theta_grid和omega_grid都是1024*180,theta_grid每一列都是一样的角度，omega_grid每一行都是一样的位置

omega_image = omega_sino;                                               % 根据极半径大小，建立笛卡尔坐标系
[omega_grid_x, omega_grid_y] = meshgrid(omega_image, omega_image);      % 网格后，omega_grid_x和omega_grid_y都是1024*1024，omega_grid_x每一列的x坐标一样，omega_grid_y每一行的y坐标都一样
[coo_th_fft2, coo_r_fft2] = cart2pol(omega_grid_x,omega_grid_y);        % coo_th_fft2 = atan(omega_grid_y/omega_grid_x), coo_r_fft2 = sqrt(omega_grid_y^2+omega_grid_x^2)     
coo_r_fft2 = coo_r_fft2.*sign(coo_th_fft2);                             % 第一象限和第二象限内，笛卡尔半径为负
coo_th_fft2(coo_th_fft2<0) = coo_th_fft2(coo_th_fft2<0) + pi;           % 第二象限和第四象限，笛卡尔角度为0~pi/2；第一象限和第三象限，笛卡尔角度为pi/2~pi

%% 根据极坐标处的值，二维插值得出笛卡尔坐标处对应的值，插值后也是1024*1024的复数矩阵
Fourier2_radial = interp2(theta_grid,omega_grid,f_p,coo_th_fft2,coo_r_fft2,'cubic',(0+1i.*0));

%% 根据插值结果，进行二维傅里叶逆变换，得出目标图像 
crop = L_pad;
target = fftshift(ifft2(ifftshift(Fourier2_radial)));                   % 二维傅里叶逆变换，得到的仍然是1024*1024的矩阵
target = target(crop+1:end-crop,crop+1:end-crop);                       % 原始数据大小256*256，因此需要对上面的数据截断
I_a = abs(target);                                                      % 复数的模值 
I_a = (I_a-min(I_a(:)))./(max(I_a(:))-min(I_a(:)));                     % 归一化为 0~1
Lg_I_a = log(1+I_a);                                                    % 为了好的显示效果，取对数

figure, 
subplot(1,2,1), imshow(flipud(Lg_I_a),[ ]); title('DFR重建的图像')
subplot(1,2,2), imshow(P,[ ]); title('原始的图像')

