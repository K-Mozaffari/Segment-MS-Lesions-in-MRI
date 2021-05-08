% levelset section
% This Matlab code demonstrates an edge-based active contour model as an application of
%  the Distance Regularized Level Set Evolution (DRLSE) formulation in the following paper:
%
%  C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation",
%     IEEE Trans. Image Processing, vol. 19 (12), pp. 3243-3254, 2010.
%
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com
%         li_chunming@hotmail.com
% URL:  http://www.imagecomputing.org/~cmli//
%end of levelset section
close all;
clear all
clc
load  Brain_ms.mat
Img=t1;
Img=double(255*(Img(:,:,1)/max(max(Img(:,:,1)))));
imshow((Img),[])
[xx yy]=size(Img);

%% parameter setting
timestep=5;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
iter_inner=45;
iter_outer=35;
lambda=5; % coefficient of the weighted length term L(phi)
alfa=1.5;  % coefficient of the weighted area term A(phi)
epsilon=1.5; % papramater that specifies the width of the DiracDelta function
sigma=1.5;     % scale parameter in Gaussian kernel
Gamma=fspecial('gaussian',15,sigma);
Img_smooth=conv2(Img,Gamma,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.
% initialize LSF as binary step function
c0=2;
initialLSF=c0*ones(size(Img));
% generate the initial region R0 as a rectangle
initialLSF(1:xx-5, 1:yy-5)=-c0;
phi=initialLSF;
% figure(1);
mesh(-phi);   % for a better view, the LSF is displayed upside down
hold on;  contour(phi, [0,0], 'r','LineWidth',2);
title('Initial level set function');
% view([-80 35]);
figure(2);
imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
title('Initial zero level contour');
pause(0.5);
potential=2;
if potential ==1
    potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model
elseif potential == 2
    potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
else
    potentialFunction = 'double-well';  % default choice of potential function
end
% start level set evolution
for n=1:iter_outer
    phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
    if mod(n,2)==0
        figure(2);
        imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
    end
end
% refine the zero level contour by further level set evolution with alfa=0
alfa=0;
iter_refine = 10;
phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
finalLSF=phi;
figure(2);
imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
hold on;  contour(phi, [0,0], 'r');
str=['Final zero level contour, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
title(str);
pause(1);
figure;
mesh(-finalLSF); % for a better view, the LSF is displayed upside down
hold on;  contour(phi, [0,0], 'r','LineWidth',2);
str=['Final level set function, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
title(str);
axis on;
figure;
% imshow(adf,[])
% hold on
kk=contour(phi, [0,0], 'r');
saveas(gcf,'tmpim.jpg')

lim=imread('tmpim.jpg');

figure;imshow(lim);
rec =[118.5100   52.5100  670.9800  528.9800];
[clim ]=imcrop(lim,rec);
slim=clim(:,:,1)-clim(:,:,3);
figure;imshow(slim)
tslim=slim>20;
BW2 = imfill(tslim,'holes');
BW4=imerode(BW2,strel('disk',62,4));
BW3=imresize(BW4,[xx yy]);
save('BW3','BW3');
figure(11);imshow(BW3)
figure(77777);imshow(Img,[])

Orginalimage = Img;
Mask = BW3;
NC = 3;
minimum = min(Orginalimage(:));
Orginalimage = Orginalimage - minimum;
maximum = max(Orginalimage(:));

X = zeros(size(Orginalimage, 1), size(Orginalimage, 2));

positions = and(Orginalimage <= maximum, Mask == 1);
X(positions) = 3;

positions = and(Orginalimage <= 2 * maximum / 3, Mask == 1);
X(positions) = 2;

positions = and(Orginalimage <= maximum / 3, Mask == 1);
X(positions) = 1;
u = zeros(1, NC);
seg = zeros(1, NC);
beta = 0;