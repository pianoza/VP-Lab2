%% Part 1
% step 1
au=557.0943; av=712.9824; u0=326.3819; v0=298.6679;
f=80; % mm
tx=100; ty=0; tz=1500; %mm
Phix=0.8*pi/2; Phiy=-1.8*pi/2; Phix1=pi/5;
pN = 6; % number of points;
%intrinsic transformation matrix
intMat = [au 0 u0 0;
          0 av v0 0;
          0  0  1 0];

%extrinsic transformation matrix; Use Euler XYX1
rx = [1 0 0;
      0 cos(Phix) -sin(Phix);
      0 sin(Phix) cos(Phix)];

ry = [cos(Phiy) 0 sin(Phiy);
      0 1 0;
      -sin(Phiy) 0 cos(Phiy)];
rx1 = [1 0 0;
       0 cos(Phix1) -sin(Phix1);
       0 sin(Phix1) cos(Phix1)];
R = rx*ry*rx1;
extMat = [R(1,:), tx;
          R(2,:), ty;
          R(3,:), tz;
          0 0 0 1];
%Step 3
% 6 random points in 3d space in range [-480:480;-480:480;-480:480].
points3d = [randi([-480, 480], [pN, 1]), randi([-480, 480], [pN, 1]), randi([-480, 480], [pN, 1])];
% homogeneous points3d
hgPoints = points3d';
hgPoints(4,:) = 1;
%Step 4
proj = intMat*extMat*hgPoints;
%Step 5
%Plot the projected points
projNorm = zeros(2, pN);
% normalize the projected points
projNorm(1,:) = proj(1,:)./proj(3,:);
projNorm(2,:) = proj(2,:)./proj(3,:);
figure, scatter3(points3d(:,1), points3d(:,2), points3d(:,3), 'bo'), title('3D points');
figure, scatter(projNorm(1,:), projNorm(2,:), 'bo', 'MarkerFaceColor', 'b'), title('projected 2D points');

%Step 6
% Get Hall transformation matrix
A = HallMatrix(points3d, projNorm);
%Step 7
tr = intMat * extMat;
tr = tr(:,:)/tr(3, 4);
% A and tr are the same;

%Step 8-9
std = 0.5; % standard deviation
noisePoints2d = zeros(2, pN);
% function randn gives normal (Gaussian) distributed random numbers;
noisePoints2d(1,:) = projNorm(1,:)+std*randn(1, pN);
noisePoints2d(2,:) = projNorm(2,:)+std*randn(1, pN);
noiseA = HallMatrix(points3d, noisePoints2d);

%get new 2d points using noiseA
noiseProj = noiseA*hgPoints;
noiseProj(1,:) = noiseProj(1,:)./noiseProj(3,:);
noiseProj(2,:) = noiseProj(2,:)./noiseProj(3,:);
noiseProj(3,:) = [];

figure;
projFig = scatter(projNorm(1,:), projNorm(2,:), 'bo');
title('Step 8'), hold on;
noiseFig = scatter(noisePoints2d(1,:), noisePoints2d(2,:), 'r+');
noiseProjFig = scatter(noiseProj(1,:), noiseProj(2,:), 'g.'); 
hold off;
legend([projFig,noiseFig,noiseProjFig],'Initial points without noise','Points with noise',...
    'Projected points','NorthEastOutside');

error = mean(sqrt((noiseProj(1,:)-projNorm(1,:)).^2+((noiseProj(2,:)-projNorm(2,:)).^2)));
fprintf('Step 9: Number of points: %d Mean error: %.16f\n',pN,error);

%Part 2
% Step 10
[fougerasTransMatr, intMatFoug, extMatFoug] = FougerasMatrix(points3d, projNorm);
% Step 11
std = [0.5, 1.0, 1.5];
for i = 1:length(std)
    noisePoints2d = zeros(2, pN);
    noisePoints2d(1,:) = projNorm(1,:)+std(i)*randn(1, pN);
    noisePoints2d(2,:) = projNorm(2,:)+std(i)*randn(1, pN);
    [f1] = FougerasMatrix(points3d, noisePoints2d);
    noiseProj = f1*hgPoints;
    noiseProj(1,:) = noiseProj(1,:)./noiseProj(3,:);
    noiseProj(2,:) = noiseProj(2,:)./noiseProj(3,:);
    noiseProj(3,:) = [];
    error = sqrt((noiseProj(1,:)-projNorm(1,:)).^2+((noiseProj(2,:)-projNorm(2,:)).^2));
    fprintf('Step 11 (Fougeras): Number of points: %d Mean error %.16f (sigma=%f)\n',pN,mean(error),std(i));
    [f2] = HallMatrix(points3d, noisePoints2d);
    noiseProj = f2*hgPoints;
    noiseProj(1,:) = noiseProj(1,:)./noiseProj(3,:);
    noiseProj(2,:) = noiseProj(2,:)./noiseProj(3,:);
    noiseProj(3,:) = [];
    error = sqrt((noiseProj(1,:)-projNorm(1,:)).^2+((noiseProj(2,:)-projNorm(2,:)).^2));
    fprintf('Step 11 (Hall): Number of points: %d Mean error %.16f (sigma=%f)\n',pN,mean(error),std(i));
end;

%Part 3
% Step 12
% camera coordinate system origin at 0 0 0;
figure, scatter3(0,0,0,'b.'), hold on, title('World scene simulation'); 
%axis vis3d; % uncomment this to keep aspect ratio
text(0, 0, 0, 'Camera');
% draw x and y axis of camera coord system;
line([0, 150], [0 0], [0 0]); text(155, 0, 0, 'x');
line([0, 0], [0 150], [0 0]); text(0, 155, 0, 'y');
% define world coordinate system axis definition;
WorldCoord = [100 0 0; 0 100 0; 0 0 100; 1 1 1];
% World coordinate origin with respect to camera coord system;
Wcenter= extMat*[0;0;0;1];
% world coord system with respect to camera
WorldCoord2Cam = extMat*WorldCoord;
% draw world coordinate system axis and give labels;
line([Wcenter(1) WorldCoord2Cam(1, 1)],[Wcenter(2) WorldCoord2Cam(2, 1)],[Wcenter(3) WorldCoord2Cam(3, 1)]);
line([Wcenter(1) WorldCoord2Cam(1, 2)],[Wcenter(2) WorldCoord2Cam(2, 2)],[Wcenter(3) WorldCoord2Cam(3, 2)]);
line([Wcenter(1) WorldCoord2Cam(1, 3)],[Wcenter(2) WorldCoord2Cam(2, 3)],[Wcenter(3) WorldCoord2Cam(3, 3)]);
text(Wcenter(1), Wcenter(2), Wcenter(3), 'World Coordinate System');
text(WorldCoord2Cam(1, 1), WorldCoord2Cam(2, 1), WorldCoord2Cam(3, 1), 'X');
text(WorldCoord2Cam(1, 2), WorldCoord2Cam(2, 2), WorldCoord2Cam(3, 2), 'Y');
text(WorldCoord2Cam(1, 3), WorldCoord2Cam(2, 3), WorldCoord2Cam(3, 3), 'Z');

% 3D point coordinates with respect to camera;
World2Cam = extMat*hgPoints;
% plot them;
scatter3(World2Cam(1,:), World2Cam(2,:), World2Cam(3,:), 'ro');

% define image plane with size 640x480
imagePlane = [
        -320, 320, 320, -320, -320;
         240, 240,-240, -240, 240;
         f, f, f, f, f
    ];
% convert from pixels to mm.
imagePlane(1,:) = imagePlane(1,:)/ku;
imagePlane(2,:) = imagePlane(2,:)/kv;
% plot the image plane.
plot3(imagePlane(1,:), imagePlane(2,:), imagePlane(3,:), 'b-');

% projNorm - are the pixels on the image plane.
RD = zeros(2, pN);
RD(1,:) = u0;
RD(2,:) = v0;
% convert 2d point coordinates on image plane into focal point
RD(1,:) = RD(1,:) - projNorm(1,:);
RD(2,:) = RD(2,:) - projNorm(2,:);
ku = -au/f;
kv = -av/f;
CU = zeros(3, pN);
% convert 2d coordinates into 3d points with respect to the camera
CU(1,:) = RD(1,:)/ku;
CU(2,:) = RD(2,:)/kv;
CU(3,:) = f;
% plot them
scatter3(CU(1,:), CU(2,:), CU(3,:), 'r.');

% rays from 3d points to camera center
for i = 1:pN
    line([0 World2Cam(1,i)], [0 World2Cam(2,i)], [0 World2Cam(3,i)]);
end;


