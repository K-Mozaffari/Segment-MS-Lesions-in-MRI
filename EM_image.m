function [mu_new,Cov_new,Alpha_new,w] = EM_image(x,mu_old,Cov_old,Alpha_old,iter)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
K=length(Alpha_old);
N=size(x,1);
d=size(x,2);
mu_new=zeros(d,K);
Cov_new=zeros(d,d,K);
Alpha_new=zeros(1,K);
w=zeros(N,K);

for it=1:iter
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& %
%        Computing Membership Weights                %
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& %    
  for k=1:K
    for i=1:N
         sum_Pm=0;
         for l=1:K
           sum_Pm=Pk(x(i,:)',mu_old(:,l),Cov_old(:,:,l))*Alpha_old(l)+sum_Pm;
         end
          w(i,k)=Pk(x(i,:)',mu_old(:,k),Cov_old(:,:,k))*Alpha_old(k)/sum_Pm;
    end
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& %
%        Computing new mean and Alpha                %
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& %
    Alpha_new(k)=sum(w(:,k))/N;
    mu_new(:,k)=sum(repmat(w(:,k),1,3).*x)'/sum(w(:,k))';
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& %
%        Computing new Covariance                    %
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& %
    s=0;
    for i=1:N
       s= (1/sum(w(:,k)))*w(i,k)*(x(i,:)'-mu_new(:,k))*(x(i,:)'-mu_new(:,k))'+s; 
    end
    Cov_new(:,:,k)=s;
  end
  Alpha_old=Alpha_new;
  mu_old =mu_new;
  Cov_old=Cov_new;
end



