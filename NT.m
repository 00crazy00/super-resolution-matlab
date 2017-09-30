function Z1=NT(X,Y,Z,X1,Y1)%X,Y,ZΪ��ֵ�����꣬TΪ����ֵ������������꣬PΪT�Ĳ�ֵ���
M=Z';
[n,m]=size(M);%ȡ�ò�ֵ����Ĵ�С
%���ȶ�x����newton����
x=Y(:,1)';%��Ϊת�ú�X,YҲ��ת�ˡ�
M1=zeros(n,m);
for i1=1:n
    L=M(i1,:);
    S=netwon(x,L);%��ÿһ�н���ţ�ٵ���
    M1(i1,:)=S;
end
%Ȼ���y�������thiele����
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
%����M2���е�ֵ������y����ֵ��ֵ
y1=X1(1,:);
n1=size(y1,2);
%���y������
A2=zeros(m,n1);
for j1=1:m
    L1=M2(:,j1);
    L1=L1';
    L2=thieledis1(y1,L1,y,P(1,j1));
    A2(j1,:)=L2;
end
%���㣨x,y������ֵ��;
x1=Y1(:,1)';
m1=size(x1,2);
Z1=zeros(n1,m1);
%��ȡx����
for i1=1:n1
    l1=A2(:,i1);
    l1=l1';
    Z=netwondis(x1,l1,x);
    Z=Z';
    Z1(i1,:)=Z;
end
Z1=Z1';
end
function C=thieledis1(X,Y,Z,p)%XΪ����ֵ�㣬YΪ����̣�ZΪԭ��ֵ�ڵ�����꣬pΪ���������
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
function c=thieledis(X,Y,Z)%xΪ�����ֵ��yΪthiele����ʽϵ����zΪԭ��ֵ��ֵ
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