function [t1,t2,pd]=NORM_BRAIN(t1,t2,pd,a)
maxt1=max(t1,[],"all");
mint1=min(t1,[],"all");

maxt2=max(t2,[],"all");
mint2=min(t2,[],"all");

maxpd=max(pd,[],"all");
minpd=min(pd,[],"all");
t1 = a.* (((t1)-mint1)) ./ (maxt1-mint1);
t2 =  a.*(((t2)-mint2)) ./ (maxt2-mint2);
pd =  a.*(((pd)-minpd)) ./ (maxpd-minpd);

end
