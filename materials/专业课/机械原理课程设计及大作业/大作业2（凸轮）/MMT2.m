% authored by Wang Rui
%受限于代码水平，看个乐呵就行，估计也没人看得懂这坨屎（写完我自己也看不懂了^_^）

r0 = 22;  % 基圆半径 
h = 35;   % 推杆行程 
e = 14; % 偏心距 
s0 = sqrt(r0 * r0 + e * e); 
rr = 18; % 滚子半径
N = 36; % 数据点个数
pi = 3.1415926; % 圆周率
ddelta = 2 * pi / N;
delta1 = 0:ddelta:pi; % 推程 
delta2 = pi:ddelta:25*pi/18;  % 远休程 
delta3 = (25*pi/18):ddelta:(33*pi/18); % 回程角 
delta4 = (33*pi/18):ddelta:2*pi; % 近休止角 
rdeg = [delta1, delta2(2:end), delta3(2:end), delta4(2:end-1)];
s1 = h * (delta1/pi - 1/(2*pi)*sin(delta1)); % 推程推杆位移
s2 = ones(1, numel(delta2)-1) * h; % 远休程推杆位移
s3 = h -  h / pi / pi / 16 * 81 * (delta3-25*pi/18) .* (delta3-25*pi/18); % 回程推杆位移
s4 = ones(1, numel(delta4)-1) * 0; % 近休程推杆位移
s = [s1, s2(2:end), s3(2:end), s4]; % 推杆位移
disp(numel(s));
x = (s + s0) .* sin(rdeg) + e .* cos(rdeg);
y = (s + s0) .* cos(rdeg) - e .* sin(rdeg);
dx = diff(x);
dx = [0,dx];
dx = dx / ddelta;
dy = diff(y);
dy = [0,dy];
dy = dy / ddelta;
%disp(cosbeta);
cosbeta = -dy ./ sqrt(dx.*dx+dy.*dy);
sinbeta =  dx ./ sqrt(dx.*dx+dy.*dy);
xr = x - rr * cosbeta;
yr = y - rr * sinbeta;
ds1 = diff(s1);
ds3 = diff(s3);
ds1 = ds1 / ddelta;
ds3 = ds3 / ddelta;
e1 = ones(1, numel(s1)-1);
e3 = ones(1, numel(s3)-1);
alpha1 = atan2(ds1 + e1 * e , s0 * e1 + s1(2:end));
alpha3 = atan2(ds3 - e3 * e , s0 * e3 + s3(2:end));
[alpha1max,deltaalpha1max] = max(alpha1); % 推程最大压力角
[alpha3max,deltaalpha3max] = max(abs(alpha3)); % 回程最大压力角
alpha1max = alpha1max * 180 / pi;
alpha3max = alpha3max * 180 / pi;
deltaalpha1max = delta1(deltaalpha1max);
deltaalpha3max = delta3(deltaalpha3max);

disp(alpha1max);
disp(alpha3max);
disp(deltaalpha1max);
disp(deltaalpha3max);

ddx = diff(dx);
ddx = [0,ddx];
ddx = ddx / ddelta;
ddy = diff(dy);
ddy = [0,ddy];
ddy = ddy / ddelta;

rho = power(dx.*dx + dy.*dy, 1.5) ./ abs(dx.*ddy - dy.*ddx);

disp(rho);
% x1 = sqrt(s0*s0+e*e)*cos(delta);
% y1 = sqrt(s0*s0+e*e)*sin(delta);
% x2 = sqrt((s0+h)*(s0+h)+e*e)*cos(delta);
% y2 = sqrt((s0+h)*(s0+h)+e*e)*sin(delta);

figure;
plot(x, y, 'r-', 'LineWidth', 2);
hold on;
plot([x(N), x(1)], [y(N), y(1)], 'r-', 'LineWidth', 2);
hold on;
plot(xr, yr, 'b-', 'LineWidth', 2);
hold on;
plot([xr(N), xr(1)], [yr(N), yr(1)], 'b-', 'LineWidth', 2);
% hold on;
% plot(x1, y1, 'r-', 'LineWidth', 2);
% hold on;
% plot(x2, y2, 'g-', 'LineWidth', 2);
axis equal;
title('凸轮轮廓');
xlabel('X (mm)');
ylabel('Y (mm)');
grid on;

alpha1allow = 35; % 推程许用压力角
alpha2allow = 65; % 回程许用压力角

[rhomin,deltarhomin] = min(rho); % 理论轮廓最小曲率半径
rhoaminallow = 0.35 * rr; % 实际轮廓最小曲率许应半径
deltarhomin = rdeg(deltarhomin);






