 function [B1,X]=tuxiangcj(X,apha,p)%XΪ������ͼƬ��aphaΪ�Ŵ������ȼ���Ϊ2,����ƽ̹��ϵ��
%Ԥ�����Ƚ�X�淶�����ܽ���4*4�Ŀ�����
X=rgb2gray(X);
X=im2double(X);
[x,y]=size(X);
A=X;
if(mod(x,2)~=0)
    A=[A;A(x,:)];
end
if(mod(y,2)~=0)
    A=[A,A(:,y)];
end
%��ͼ����зŴ�
[x1,y1]=size(A);
x2=x1*apha;y2=y1*apha;
B=zeros(x2,y2);%���ɲ�ֵ�����ģ��
for i=apha:apha:x2
    for j=apha:apha:y2
        B(i,j)=A(i/apha,j/apha);
    end
end
%��ԭͼ���ֵ��Ӧ������ͼ���λ��
B1=B;
%��ͼ�����������ƽ̹���ж�;
[xc,yc]=meshgrid(linspace(0.1,0.4,4),linspace(0.1,0.4,4));
[xxc,yyc]=meshgrid(linspace(0.05,0.4,8),linspace(0.05,0.4,8));
for i=1:2:x1-2
    for j=1:2:y1-2
        L=[A(i,j),A(i,j+1),A(i,j+2),A(i,j+3);
            A(i+1,j),A(i+1,j+1),A(i+1,j+2),A(i+1,j+3);
            A(i+2,j),A(i+2,j+1),A(i+2,j+2),A(i+2,j+3);
            A(i+3,j),A(i+3,j+1),A(i+3,j+2),A(i+3,j+3)];
        x=panduan(L,p);
        if(x==1)%��ƽ̹��ʹ��ţ��ţ�ٷ�
            
           L1=netnetwon(xc,yc,L,xxc,yyc);%X,Y,ZΪ��ֵ�����꣬TΪ����ֵ������������꣬PΪT�Ĳ�ֵ���
        else
           L1=NT(xc,yc,L,xxc,yyc);%X,Y,ZΪ��ֵ�����꣬TΪ����ֵ������������꣬PΪT�Ĳ�ֵ���  
        end
        if(i==1&&j==1)
        B1(2*i-1:2*i+6,2*j-1:2*j+6)=L1;
        end
        if(i==1&&j~=1)
        L2=(L1(1:8,1:4)+B1(2*i-1:2*i+6,2*j-1:2*j+2))./2;
        L1=[L2,L1(1:8,5:8)];
        B1(2*i-1:2*i+6,2*j-1:2*j+6)=L1;
        end 
        if(i~=1&&j==1)
        L2=(L1(1:4,1:8)+B1(2*i-1:2*i+2,2*j-1:2*j+6))./2;
        L1=[L2;L1(5:8,1:8)];
        B1(2*i-1:2*i+6,2*j-1:2*j+6)=L1;
        end
        if(i~=1&&j~=1)
        L2=(L1(1:4,1:8)+B1(2*i-1:2*i+2,2*j-1:2*j+6))./2;
        L3=(L1(1:8,1:4)+B1(2*i-1:2*i+6,2*j-1:2*j+2))./2;
        L4=(L2(1:4,1:4)+L3(1:4,1:4))./2;
        L1=[L4,L2(1:4,5:8);L3(5:8,1:4),L1(5:8,5:8)];
        B1(2*i-1:2*i+6,2*j-1:2*j+6)=L1;
        end
    end
end
    B1=im2uint8(B1);
  imwrite(B1,'E:\test\new002.bmp');
  imwrite(X,'E:\test\001.bmp');
end