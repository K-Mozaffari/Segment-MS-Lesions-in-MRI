 function Result=Pk(x,mu,Cov)
   Dim=length(mu);
   Result=(1/(((2*pi)^(Dim/2))*sqrt(det(Cov))))*exp(-0.5*(x-mu)'*Cov^-1*(x-mu));
 end