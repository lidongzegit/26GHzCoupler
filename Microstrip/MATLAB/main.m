close all
clear
clc

%% 微带线仿真
% 常数赋值
eps0 = 8.854e-12;
U = 1;
% 读入从COMSOL中导出的数据文件
no2xy = importdata('node.txt',' '); %每个点的xy坐标
el2no = importdata('elem.txt',' '); %每个单元的三个顶点编号
noIn = []; %微带线节点编号
noEx = []; %地节点编号

noNum = size(no2xy,1); %节点数
elNum = size(el2no,1); %单元数

K = zeros(noNum);
P = zeros(noNum,1);

for i = 1:noNum
    x = no2xy(i,1);
    y = no2xy(i,2);
    if (y == 0.25 && x<=0.36)||(x == 0.36 && y >= 0.25 && y <= 0.285)||(y == 0.285 && x<=0.36)
        noIn = [noIn;i];
    end
    if (y == 0)||(x == 2.5)||(y == 1)
        noEx = [noEx;i];
    end
end

% 计算总体系数矩阵K
for elIdx = 1:elNum
    no = el2no(elIdx,:);
    xy = no2xy(no,:);
    K_el = CmpElMtx(xy);
    K(no,no) = K(no,no)+K_el;
end

% 处理强加边界条件
no_known = union(noIn,noEx);
no_all = 1:noNum;
no_ess = setdiff(no_all,no_known);
K_right = K(no_ess,no_known);
K_ess = K(no_ess,no_ess);
P = P(no_ess);
phi = zeros(length(no_all),1);
phi(noIn) = U*ones(length(noIn),1);
phi_known = phi(no_known);
phi_ess = K_ess\(P-K_right*phi_known);

phi = zeros(length(no_all),1);
phi(no_known) = phi_known;
phi(no_ess) = phi_ess;

figure(1)
trisurf(el2no,no2xy(:,1),no2xy(:,2),phi,'FaceColor','interp','LineWidth',0.1);
view([0,90]);
patch([0 2 2 0],[0 0 0.25 0.25],[1 1 1 1],'b','FaceAlpha',0,'EdgeColor','w');
patch([0 0.36 0.36 0],[0.25 0.25 0.285 0.285],[1 1 1 1],'k','FaceAlpha',1);

% 求解场量
EST = zeros(noNum,1);
ST = zeros(noNum,1);
for elIdx = 1:elNum
    no = el2no(elIdx,:);
    xy = no2xy(no,:);
    z = phi(no);
    G(elIdx,:) = CmpGElMtx(xy,z);
    EST(no,1) = EST(no,1) + G(elIdx,1)*G(elIdx,2)*ones(3,1);
    ST(no,1) = ST(no,1) + G(elIdx,2)*ones(3,1);
end
E = EST./ST;

figure(2)
trisurf(el2no,no2xy(:,1),no2xy(:,2),E,'FaceColor','interp','LineWidth',0.1);
view([0,90]);





