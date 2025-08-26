clc
clear
% 读取输入数据
time=xlsread('时间.xlsx');
I=xlsread('电流实测值.xlsx');
% 假设电流实测值保存在变量measured_currents中
% 假设点堆动力学参数保存在变量params中
params=struct('beta',[2.16287E-04 1.46220E-03 1.35047E-03 2.81505E-03 9.54088E-04 3.22744E-04
],'lamda',[1.24988E-02	3.08168E-02	1.15258E-01	3.11078E-01	1.24124E+00	3.33321E+00], ...
'rho_0',0,'Lambda',2.66315E-05);


%% 
% 初始化参数
beta = params.beta; % 延迟中子比例
lamda = params.lamda; % 
Lambda = params.Lambda;
rho_0 = params.rho_0; % 初始反应性
sumbeta=sum(beta);

N = length(time); % 数据点数量
dt=zeros(N,1); % 无拟合的非均匀时间步长
% 初始化反应性数组
rho = zeros(N,1);
rho(1) = rho_0;
for i = 2:N 
    dt(i)=time(i)-time(i-1);
end
xi=zeros(N,length(beta));
d1=zeros(N,length(beta));
d2=zeros(N,length(beta));
d3=zeros(N,length(beta));
% for i= 1:6
%     xi(i,1)=I(1)/Lambda(i)*(1-exp(-Lambda(i)*dt(2)));
% end
xi(1,:)=I(1)./lamda(:).*(1-exp(-lamda(:).*dt(2)));
d1(1,:)=exp(-lamda(:).*dt(2));
d2(1,:)=1./lamda(:).*(1-1./(lamda(:).*dt(2)).*(1-exp(-lamda(:).*dt(2))));
d3(1,:)=1./lamda(:).*(exp(-lamda(:).*dt(2))-1./(lamda(:).*dt(2)).*(1-exp(-lamda(:).*dt(2))));
%% 

% 计算反应性随时间的变化
for i = 2:N
    
    % 根据离散式逆点堆方程计算反应性
    %drho_dt = (measured_currents(i) - measured_currents(i-1)) / (beta * Lambda * dt);
    first  = Lambda*(I(i)-I(i-1))/(sumbeta*dt(i)*I(i));
    second = I(1)*( sum(beta(:).*exp(-lamda(:).*dt(i))) )/(I(i)*sumbeta);
    third_first = xi(i-1,:).*d1(i-1,:) + I(i).*d2(i-1,:)+I(i-1).*d3(i-1,:);
    third  = 1/I(i)*( sum(lamda(:).*beta(:).*(third_first(:))) )/sumbeta;
    rho(i) = sumbeta*(1 + first - second - third);
    %更新下一步
    if i==N
        break
    end
    xi(i,:)=xi(i-1,:).*d1(i-1,:)+I(i).*d2(i-1,:)+I(i-1).*d3(i-1,:);
    d1(i,:)=exp(-lamda(:).*dt(i+1));
    d2(i,:)=1./lamda(:).*(1-1./(lamda(:).*dt(i+1)).*(1-exp(-lamda(:).*dt(i+1))));
    d3(i,:)=1./lamda(:).*(exp(-lamda(:).*dt(i+1))-1./(lamda(:).*dt(i+1)).*(1-exp(-lamda(:).*dt(i+1))));
end

% 绘制反应性随时间的变化曲线
%time = (0:N-1) * dt;
time=time-100;
plot(time(100:205), rho(100:205),'-r','LineWidth',2);
hold on
plot(time(205:341), rho(205:341),'-b','LineWidth',2);
xlabel('时间 (s)');
ylabel('反应性');
title('反应性随时间的变化');
grid on;
xticks(0:25:200)
yticks(-0.3:0.05:0.05)
legend('正反应性','负反应性');  % 设置图例
axis([0 205 -0.3 0.05])
xlswrite('反应性.xlsx',rho,'sheet1','A1') 