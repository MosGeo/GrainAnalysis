%% Read the Folders for STL files
clear all
mainDirectory = pwd;
mainDirectory = '';
grainDirectoryName = '06 Grain_.5';

grainDirectoryName = 'Real_.2';

stlDirectoryName = ''
resultsFileName = 'Results.dat';

grainDirectory = fullfile (mainDirectory,grainDirectoryName);
grainsListing = dir(fullfile(grainDirectory, stlDirectoryName));

grainsListing = grainsListing(3:5002);
nGrains =  numel(grainsListing);

% load Results File
nameColumn = 1;
positionColumns = [2 3 4];
rotationColumns = [5 6 7];
scaleColumns = [8 9 10];
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s';
resultsFileName = fullfile(mainDirectory, grainDirectoryName, resultsFileName);
fileID = fopen(resultsFileName,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);
cell2num = @(inputCell) arrayfun(@(x) str2double(x), inputCell);
grainNames = dataArray{nameColumn};
grainPosition = [cell2num(dataArray{positionColumns(1)}) cell2num(dataArray{positionColumns(2)}) cell2num(dataArray{positionColumns(3)})];
grainRotation = [cell2num(dataArray{rotationColumns(1)}) cell2num(dataArray{rotationColumns(2)}) cell2num(dataArray{rotationColumns(3)})];
grainScale = [cell2num(dataArray{scaleColumns(1)}) cell2num(dataArray{scaleColumns(2)}) cell2num(dataArray{scaleColumns(3)})];
clear nameColumn positionColumns rotationColumns scaleColumns formatSpec fileID resultsFileName dataArray
% LOAD GRAINS
grains(nGrains) = Grain(nGrains);
tic
for grainNumber = 1:nGrains
    disp(['Loading Grain ' num2str(grainNumber) ' Out of ' num2str(nGrains)])
    selectedGrainName = grainsListing(grainNumber).name;
    grainNameExtract = strsplit(selectedGrainName, ' - ');
    grainNameExtract = grainNameExtract{2};
    [~, ind] = ismember(grainNameExtract(1:end-4), grainNames);
    grainFileName = fullfile(grainDirectory, stlDirectoryName, selectedGrainName);
    [vertices,faces,normals,~] = stlRead(grainFileName);
    
    grains(grainNumber).setMeshData(vertices, faces, normals);
    grains(grainNumber).setName(selectedGrainName);
    grains(grainNumber).setPRSdata(grainPosition(ind,:), grainRotation(ind,:), grainScale(ind,:))
    grains(grainNumber).translateGrain(grainPosition(ind,:));
end
clear selectedGrainName ind vertices faces normals grainPosition grainRotation grainScale 
toc

% CREATE GRAIN PACK
grainPack = GrainPack(grains);
grainPack.createBinaryGrainPack(.05);
im = grainPack.extractSubVolume(.05,.05,[.05 .1]);
%grainPack.resetSubVolume();
%grainPack.visulize3D();
%stackFileName = [grainDirectory '\Image Stack\GrainPack'];
%grainPack.saveBinaryImageStack(stackFileName, 1);
grainPackFileName = fullfile(grainDirectory, [grainDirectoryName, '_grainPack.mat']);
save (grainPackFileName, 'grainPack')

%porosity = grainPack.calculatePorosity()
%grainPack.visulizeSlices( .5, .5, .5);

%%
im = grainPack.extractSubVolume(.1,.1,[.1 .5]);
im(im==0) = 2;
[x, y, z] = size(im);
clear opt
opt(1).radbound=1; % head surface element size bound
opt(2).radbound=1; % brain surface element size bound
opt(1).side='lower'; %
opt(2).side='lower'; %
maxvol = 100;
dofix = 1;
isovalues =  [1];
method = 'cgalsurf';
[node,elem,face]=v2m(uint8(im),isovalues,opt,maxvol,'cgalmesh');

 plotmesh(node,face);
 %savestl(node(:,1:3),elem,'comsol.stl','Grains')
 saveComsol(node, face, elem, 'comsol.mphtxt');
 %mesh2comsol(node(:,1:3),face,elem,'comsol2.mphtxt');
 
 %%
 
clear opt
opt.keepratio=0.1; % this option is only useful when vol2mesh uses 'simplify' method
opt.radbound=3;    % set the target surface mesh element bounding sphere be <3 pixels in radius.
tic
[node,elem,face]=vol2mesh(im,1:size(im,1),1:size(im,2),1:size(im,3),opt,100,1);
toc

%%
bw = grainPack.extractSubVolume(.1,.1,[.1 .1]);
n= 0;
bw = bw==n;
%grainPack.visulizeSlices( .5, .5, .5, double(bw));
D = -bwdist(~(bw));
%D = imgaussfilt3(D,1);
%grainPack.visulizeSlices( .5, .5, .5, D);
connectivity  = 26;
mask = imextendedmin(D,2,connectivity);
D2 = imimposemin(D,mask);
%grainPack.visulizeSlices( .5, .5, .5, D2);
Ld2 = watershed(D2,connectivity);
Ld3 = (Ld2 == 0);
for (i = 1:size(Ld3,3))
   Ld3(:,:,i) = bwmorph(Ld3(:,:,i),'skel',Inf);
   Ld3(:,:,i) = bwmorph(Ld3(:,:,i),'thin'); 
end
% % 
% for (i = 1:size(Ld3,2))
% %    Ld3(:,i,:) = bwmorph(squeeze(Ld3(:,i,:)),'skel',Inf);
%     Ld3(:,i,:) = bwmorph(squeeze(Ld3(:,i,:)),'thin'); 
%  end
% Ld3 = Skeleton3D(Ld2==0);
% grainPack.visulizeSlices( .5, .5, .5, double(Ld3));
bw3 = bw;
bw3(Ld3) = n;
grainPack.visulizeSlices( .5, .5, .5, double(bw3));



%
% find(fullMask == 1)
% [x,y,z] = ind2sub(size(fullMask),find(fullMask == 1));
% [gm, samc] = mcurvature_vec(x,y,z)




%%
[grainImage, v2smap]=surf2vol(vertices,faces,boxRaster,boxRaster,boxRaster);
Volume=polygon2voxel(FV,[901 901 901],'auto');
