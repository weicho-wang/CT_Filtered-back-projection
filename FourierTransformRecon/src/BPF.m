clear all;
close all;

%% 获得原始图像、投影数据、以及直接反投影后的断层图像
p_N = 256; 
theta_N = 180; 
pad_N = 1024;

P = phantom(p_N);
theta = 0:(theta_N-1);
[R,xp] = radon(P,theta);
[I H] = iradon(R,theta,'linear','none');                % iradon变换，不滤波，则I就是直接反投影后的断层图像

row_pad = pad_N; 
row_I = size(I,1); 
row_start = (row_pad-row_I)/2;
fb = zeros(row_pad,row_pad);                            % 将I扩充成1024*1024
for i = 1:1:row_I
    for j = 1:1:row_I
        fb(i+row_start,j+row_start) = I(i,j);           % 里面256*256是直接反投影后的数据，外面的都是0
    end
end

%% 获取二维锥形滤波器Frow， Frow = sqrt(u^2+v^2)
N = pad_N; 
uP = zeros(N,N); 
vP = zeros(N,N); 
Frow = zeros(N,N);
for i = 1:1:N
    for j = 1:1:N
        uP(i,j) = ((-(N)/2)+(j)*1)*2*pi/N;
        vP(i,j) = ((-(N)/2)+(i)*1)*2*pi/N;
        Frow(i,j) = sqrt(uP(i,j)^2+vP(i,j)^2);
    end
end

%% 求取断层图像fb的二维傅里叶变换Fb，根据Frow,二维滤波后得出Fvu，Fvu = Fb*Frow
Fb = fft2(fb);                                          % 二维傅里叶变换
Fb = ifftshift(Fb);                                     % 直流分量移动到中央，中心是低频，向外是高频
Fvu = zeros(row_pad,row_pad);
for i=1:1:row_pad
		for j=1:1:row_pad
			if(i==(row_pad/2+1) || j==(row_pad/2+1) )	
				Fvu(i,j) = Fb(i,j);                     % 中央的直流分量，仍然保持未滤波前的值
			else
				Fvu(i,j) = Fb(i,j)*Frow(i,j);			% 二维滤波
			end
		end
end

%% 对滤波后的数据Fvu，求取二维傅里叶逆变换，得到fxy
fxy = (ifft2((Fvu)));                                   % 求取二维傅里叶逆变换
target = fxy;
crop = row_start; 
target = target(crop+1:end-crop,crop+1:end-crop);
I_a = abs(target);                                      % 复数的模值 
I_a = (I_a-min(I_a(:)))./(max(I_a(:))-min(I_a(:)));     % 归一化为 0~1
Lg_I_a = log(1+I_a);                                    % 为了好的显示效果，取对数

figure, 
subplot(1,3,1), imshow(P,[]), title('原始图像');
subplot(1,3,2), imshow(I,[]), title('直接反投影后的断层图像');
subplot(1,3,3), imshow(Lg_I_a,[]), title('先反投影后滤波恢复的图像');
