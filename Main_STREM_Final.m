clc
clear 
close all

 load  Brain_ms.mat
 load BW3.mat
 [t1,t2,pd]=NORM_BRAIN(t1,t2,pd,1);
 figure, imagesc(t1.*BW3);
title (['Removed Skull '])
% %%%%%%%%%%%%%%%%%%%%%%%%% 2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  x = zeros([size(t1) 3]); % input image that includes all t1 t2 pd data

x(:,:,1)=t1;
x(:,:,2)=t2;
x(:,:,3)=pd;

p_value=0.025;
[m,n,d]=size(x); % d= dimension 
N=m*n;
K=3; % the number of clusters
x=reshape(x,[],K);
 bw3=reshape(BW3,[],1);
 
 rw_bw3=find(bw3==1);
 
 

%     mu_old=randn(d,K);
  mu_old=[0.4504,0.2661,0.5815;0.5504,0.7803,0.4307;0.8348,0.8790,0.7341];
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& %
%  Cov_old=repmat(eye(d),1,1,K);
Cov_old(:,:,1)= [ 0.0020   -0.0012   -0.0010;   -0.0012    0.0037    0.0013;   -0.0010    0.0013    0.0016];
Cov_old(:,:,2)=[0.0092   -0.0125   -0.0021;   -0.0125    0.0180    0.0030;   -0.0021    0.0030    0.0012];
Cov_old(:,:,3)=   1.0e-03 * [.4112   -0.0947   -0.0770;   -0.0947    0.8402    0.0985;   -0.0770    0.0985    0.7836];

%   Alpha_old=[1/3,1/3,1/3]';
Alpha_old=[0.3262,0.1703,0.5035]';
iter_EM=3;

w_ml=zeros(m*n,K);
% 
 [mu_new_ml,Cov_new_ml,Alpha_new_ml,w_ml(rw_bw3,:,:)] = EM_image(x(rw_bw3,:,:),mu_old,Cov_old,Alpha_old,iter_EM);
 mu_new_init=mu_new_ml;
 Cov_new_init=Cov_new_ml;
 Alpha_new_init=Alpha_new_ml;
% 
% 
[val_ml,cluster_idx_ml]=max(w_ml,[],2);
Clustring=reshape(cluster_idx_ml,m,n).*BW3;
figure,imagesc(Clustring)

w_final=zeros(m*n,K);
iter=3;


for it=1:iter
     P_new=zeros(N,K);
     for i=1:length(rw_bw3)  
         for k=1:K
            P_new(i,k)=Pk(x(rw_bw3(i),:)',mu_new_ml(:,k),Cov_new_ml(:,:,k))*Alpha_new_ml(k);
         end
     end
     r=-log(sum(P_new,2));
     [val_s,idx_s]=sort(r);

     idx_trimmed=idx_s(1:fix(length(rw_bw3)*0.96));
     x_trimmed=x(rw_bw3(idx_trimmed),:);

    [mu_new_ml,Cov_new_ml,Alpha_new_ml,~] = EM_image(x_trimmed,mu_new_ml,Cov_new_ml,Alpha_new_ml,iter_EM);   

end
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& %
%                 Classification                     %
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& %

% 1-1 Computing mahalanobis distance

D=zeros(length(rw_bw3),K);
for i=1:length(rw_bw3)
    for k=1:K
        D(i,k)=(x(rw_bw3(i),:)'-mu_old(:,k))'*inv(Cov_new_ml(:,:,k))*(x(rw_bw3(i),:)'-mu_old(:,k));
    end
end
% 1-2- For a vector y i , if the Mahalanobis distance between each
% class j is greater than the critical value of χ 2 m distribution for a given p value,
% then the vector is considered as an outlier and belongs to c k+1 .

threshold = chi2inv(1-p_value,d);
C_k=(D>threshold);

C_k1=sum(C_k,2)>2;
C_k1_full=zeros(N,1);
C_k1_full(rw_bw3)=C_k1;
figure,imshow((reshape(C_k1_full,m,n)+1).*t1);
title (['Oultiers '])


idx_wm=3;


WM_avg=mu_new_ml(:,idx_wm);
HyperInt_val=WM_avg+(WM_avg)/5;
HyperInt_lower=WM_avg-(WM_avg)/5;
CK_F=x(rw_bw3,:).*repmat(C_k1,1,d)>repmat(HyperInt_val',length(rw_bw3),1);
CK_F3=x(rw_bw3,:).*repmat(C_k1,1,d)<repmat(HyperInt_lower',length(rw_bw3),1);


MS_lesion=sum(CK_F(:,2:3),2)>1; % it means if in T2 and PD is hyperintense then it is MS_lesion
Necrosis=CK_F3(:,1).*MS_lesion;


MS_lesion_full=zeros(size(x,1),1);

Necrosis_full=zeros(size(x,1),1);

MS_lesion_full(rw_bw3)=MS_lesion;
Necrosis_full(rw_bw3)=Necrosis;

figure,imagesc((reshape(MS_lesion_full,m,n)+1).*t1);
title (['Others Lesions'])
figure,imagesc((reshape(Necrosis_full,m,n)+1).*t1);
title (['Necrosis Lesions'])
% figure,imagesc(pd),colormap gray


