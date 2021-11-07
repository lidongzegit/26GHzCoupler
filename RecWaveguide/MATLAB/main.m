close all
clear
clc
%% 矩形波导仿真
a = 5;
b = 0.25;
Hmax = 0.1;

%创建2D矩形模型
Model = createpde(1);
R = [3,4,0,a,a,0,-b/2,-b/2,b/2,b/2]';
g = decsg(R);
geometryFromEdges(Model,g);

%线性网格剖分
mesh = generateMesh(Model,'GeometricOrder','linear','Hmax',Hmax);
figure(1);
pdeplot(Model);
xlim([0 5]);
ylim([-0.2 0.2]);

%读取节点信息
[p,e,t] = meshToPet(mesh);
Ele = t(1:3,:)';
No = p';

%获取节点编号
N_num = size(No,1);

%获取系数矩阵的值
[K,T] = KTMat(Ele,No);

%求解特征值和特征向量
[V,D] = eigs(T\K, 2, 1e-5);

% 读取特征值与特征向量
kc2 = diag(D);
if ~issorted(kc2)
       [kc2,idx] = sort(kc2);
       V = V(:, idx);
end
kc = sqrt(kc2);

%去掉特征值约为零的模
if abs(kc(1)) < 0.0001
       kc2 = kc2(2:end);
       kc = kc(2:end);
       V = V(:, 2:end);
end

%画出主模
figure(2);
Hz_all = V(:,1);
trisurf(Ele,No(:,1),No(:,2),Hz_all);
view(2);
xlabel('x'); ylabel('y'); zlabel('H_z');axis('equal');
title(sprintf('k_c=%0.5f',kc(1)));
colorbar('Location','westoutside');
subtitle('主模的纵向磁场分布');