function Z1=NT(X,Y,Z,X1,Y1)%X,Y,Z为插值点坐标，T为待插值点坐标横纵坐标，P为T的插值结果
M=Z';
[n,m]=size(M);%取得插值矩阵的大小
%首先对x方向newton递推
x=Y(:,1)';%因为转置后，X,Y也翻转了。
M1=zeros(n,m);
for i1=1:n
    L=M(i1,:);
    S=netwon(x,L);%对每一行进行牛顿递推
    M1(i1,:)=S;
end
%然后对y方向进行thiele递推
y=X(1,:);
M2=zeros(n,m);
P=zeros(1,m); 
for j1=1:m
    L1=M1(:,j1);
    L1=L1';
        [S1,p]=thiele1(y,L1);
    P(1,j1)=p;
    M2(:,j1)=S1;
end
%利用M2表中的值来计算y待插值点值
y1=X1(1,:);
n1=size(y1,2);
%获得y的坐标
A2=zeros(m,n1);
for j1=1:m
    L1=M2(:,j1);
    L1=L1';
    L2=thieledis1(y1,L1,y,P(1,j1));
    A2(j1,:)=L2;
end
%计算（x,y）待插值点;
x1=Y1(:,1)';
m1=size(x1,2);
Z1=zeros(n1,m1);
%获取x坐标
for i1=1:n1
    l1=A2(:,i1);
    l1=l1';
    Z=netwondis(x1,l1,x);
    Z=Z';
    Z1(i1,:)=Z;
end
Z1=Z1';
end
function C=thieledis1(X,Y,Z,p)%X为待插值点，Y为逆差商，Z为原插值节点横坐标，p为奇异点坐标
XH=size(X,2);
YH=size(Y,2);
if(p~=0)
    c=netwondis(X,Y(1,p+1:YH),Z(1,p+1:YH));
    c=c';
    c1=(X-Z(1,p)).*c+Y(1,p);
    for i=1:XH
        if(c1(1,i)>0.0001)
    C=thieledis(X,[Y(1,1:p-1),c1(1,i)],Z(1,1:p));
        else
         C=thieledis(X,Y(1,1:p-1),Z(1,1:(p-1)));
        end
    end
else
    C=thieledis(X,Y,Z);
end
end
function c=thieledis(X,Y,Z)%x为待求点值，y为thiele多项式系数，z为原插值点值
y=size(Y,2);
y1=size(X,2);
c=ones(1,y1);
for i=1:y1
A=ones(y+1,1); 
A(1,1)=1;A(2,1)=Y(1,1);
B=ones(y+1,1);
B(1,1)=0;B(2,1)=1;
for k=3:y+1
    A(k,1)=A(k-1,1)*Y(1,k-1)+(X(1,i)-Z(k-2))*A(k-2,1);
    B(k,1)=B(k-1,1)*Y(1,k-1)+(X(1,i)-Z(k-2))*B(k-2,1);
end
c(1,i)=A(y+1,1)/B(y+1,1);
end
end  