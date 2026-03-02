function img = DirectBackProject(F, costheta, sintheta, num, interp)
% 累加法，实现反投影

n_proj = size(costheta,2);          % 投影角度
N = num;                            % 重建图像大小
ctr = floor((N-1) / 2);             
xleft = double(-ctr);               % X坐标最小值
ytop = double(ctr);                 % Y坐标最大值

len = size(F,1);
ctr_idx = floor(len/2);             % 接收栅格的中点
ctr_idx_1 = ctr_idx + 1;            % 向前一个点
switch interp
    case 'nearest'
        InterpType = 0;
    case 'linear'
        InterpType = 1;
    otherwise
        InterpType = 1;
end
%% 每个角度下得出一个矩阵，所有角度下的矩阵叠加
img = zeros(N,N);                   % 存储重建的图像
for theta_idx = 1:1:n_proj
   cos_theta = costheta(theta_idx);
   sin_theta = sintheta(theta_idx);

   for x_idx = 1:1:N
       x = xleft + x_idx-1;
       x_cos_theta = x*cos_theta;
       t = x_cos_theta + ytop*sin_theta;    % 根据正弦图，r = x*costheta + y*sintheta
       for y_idx = 1:1:N
            if(InterpType == 0)             % 最近插值，直接选取最邻近的点
                    if(t>0)
                        a = floor(t);
                    else
                        a = ceil(t);
                    end
                    img(y_idx,x_idx) = img(y_idx,x_idx) + F(a+ctr_idx,theta_idx);
            else                            % 线性插值，根据直线的两点式方程
                    a = floor(t);
                    img(y_idx,x_idx) = img(y_idx,x_idx) + ...
                        (t-a)*(F(a+ctr_idx_1,theta_idx)-F(a+ctr_idx,theta_idx)) + F(a+ctr_idx,theta_idx);
            end
            t = t-sin_theta;
       end
   end
end

