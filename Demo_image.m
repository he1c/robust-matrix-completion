clear

Image=imread('data\lenna.bmp');
CrossMark=imread('data\CrossMark.bmp');

I=double(rgb2gray(Image))./255;
[n1,n2]=size(I);
r=50;
p=0.5;
v1=0.0001;
v2=1;
c=0.1;

LRM=Imagerankapprox(I,r);
LRM(LRM<0)=0;

Mask=1-CrossMark;

MissM_clean=LRM.*Mask;
sizeX = size(MissM_clean);

Mask_v=find(Mask(:)==1);
IDX = Mask_v; 

y = reshape(MissM_clean,n1*n2,1);
y = y(IDX);
y = y+noisemix(length(IDX),1,c,v1,v2,'gaussian');

MissM_v=zeros(size(MissM_clean));
MissM_v(IDX)=y;
MissM=reshape(MissM_v,n1,n2);

option.U=abs(0.001*rand(n1,r));
option.V=abs(0.001*rand(r,n2));
option.yita=2;
option.sigmamin=0.01;
option.maxitr=500;
    
tic
option.stop_1=1e6;
option.stop_2=1e-4;
M_HQ_PF=HQ_PF(MissM,Mask,option);
error3=norm(LRM-M_HQ_PF,'fro');
disp(['HQPF_acc error:' num2str(error3) '. Time ' num2str(toc) 's']);

tic
option.stop_1=1e6;
option.stop_2=1e-7;
M_HQ_ASD=HQ_ASD(MissM,Mask,option);
error4=norm(LRM-M_HQ_ASD,'fro');
disp(['HQASD error:' num2str(error4) '. Time ' num2str(toc) 's']);
    

