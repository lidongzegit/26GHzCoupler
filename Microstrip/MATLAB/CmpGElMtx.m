function Ge = CmpGElMtx(xy,z)

% 三角形的三边矢量
s1 = xy(2,:)-xy(1,:);
s2 = xy(3,:)-xy(2,:);
s = [s1;s2];

Z = [z(2,:)-z(1,:);z(3,:)-z(2,:)];

bc = s\Z;

% 三角形的三边矢量
s1 = xy(3,:)-xy(2,:);
s2 = xy(1,:)-xy(3,:);
s3 = xy(2,:)-xy(1,:);

Atot = 0.5*(s2(1)*s3(2)-s2(2)*s3(1)); %计算三角形单元面积，矢量叉乘

G = sqrt(bc(1)^2+bc(2)^2);

Ge = [G Atot];

