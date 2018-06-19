

grainDirectoryNames = {'RealSH_0','RealSH_.1','RealSH_.2','RealSH_.3', 'RealSH_.4',...
    'RealSH_.5', 'RealSH_.6', 'RealSH_.7', 'RealSH_.8', 'RealSH_.9', 'RealSH_1'}

grainDirectoryNames = {'RealSHMONO_0','RealSHMONO_.1','RealSHMONO_.2','RealSHMONO_.3', 'RealSHMONO_.4',...
    'RealSHMONO_.5', 'RealSHMONO_.6', 'RealSHMONO_.7', 'RealSHMONO_.8', 'RealSHMONO_.9', 'RealSHMONO_1'} %, 

nGrainDirectories = numel(grainDirectoryNames);

for i=1:nGrainDirectories
   grainDirectoryName = grainDirectoryNames{i};
   doGrainStuff(grainDirectoryName);
end

%%

grainDirectoryNames = {'RealSH_0Sph','RealSH_.1','RealSH_.2','RealSH_.3', 'RealSH_.4',...
    'RealSH_.5', 'RealSH_.6', 'RealSH_.7', 'RealSH_.8', 'RealSH_.9', 'RealSH_1'}

%grainDirectoryNames = {'RealSHMONO_0Sph','RealSHMONO_.1','RealSHMONO_.2','RealSHMONO_.3', 'RealSHMONO_.4','RealSHMONO_.5', 'RealSHMONO_.6', 'RealSHMONO_.7', 'RealSHMONO_.8', 'RealSHMONO_.9', 'RealSHMONO_1'}

amp = 0:.1:1

nGrainDirectories = numel(grainDirectoryNames);

por = [];
for i=1:nGrainDirectories
    currentFileName = fullfile(grainDirectoryNames{i}, [grainDirectoryNames{i}, '_grainPack.mat']);
    currentFileName = fullfile(pwd, grainDirectoryNames{i}, [grainDirectoryNames{i}, '_grainPack.mat']);
    load(currentFileName);
    im = grainPack.extractSubVolume(.05,.05,[.05 .15]);
    por(i) = grainPack.calculatePorosity()
end

%% Figure Plotting
figure('Color', 'White', 'Units','inches', 'Position',[0 0 10 4],'PaperPositionMode','auto');
subplot(1,2,2)
hold on
plot(amp,por,'--k', 'LineWidth', 1)
scatter(amp, por, 100, 'w', 'filled', 'o','MarkerEdgeColor','k')
box on

ylabel({'Porosity (fraction)'},... 
    'FontUnits','points',... 
    'FontWeight','normal',... 
    'FontSize',16,... 
    'FontName','Times')

xlabel({'Perlin noise amplitude'},... 
    'FontUnits','points',... 
    'FontWeight','normal',... 
    'FontSize',16,... 
    'FontName','Times')

set(gca,...
'Units','normalized',... 
'FontUnits','points',... 
'FontWeight','normal',... 
'FontSize',14,... 
'FontName','Times')

%%
subplot(1,2,1)

minX  = 0
maxX  = 3
deltaX = .1;
meanX = 1;
stdX  = 0.5;

% Calculate the varience
x = (minX:deltaX:maxX)';
varX = stdX^2;

% Calculate mu and sigma 
mu    = log(meanX/sqrt(1+varX/meanX^2));
sigma = (log(1+varX/meanX^2))^.5;
pdfX = lognpdf(x,mu, sigma);
pdfX = round(pdfX,4)/sum(round(pdfX,4));

x = x/sum(x.*pdfX)
xlim([0,maxX])
% Figure plotting
hold on
plot(x,pdfX,'--k', 'LineWidth', 1)
scatter(x, pdfX, 100, 'w', 'filled', 'o','MarkerEdgeColor','k')
box on

xlabel({'Normalized equivalent spherical diameter'},... 
    'FontUnits','points',... 
    'FontWeight','normal',... 
    'FontSize',16,... 
    'FontName','Times')

ylabel({'Propobility density'},... 
    'FontUnits','points',... 
    'FontWeight','normal',... 
    'FontSize',16,... 
    'FontName','Times')

set(gca,...
'Units','normalized',... 
'FontUnits','points',... 
'FontWeight','normal',... 
'FontSize',14,... 
'FontName','Times')
%% Sphericity Plot

grainDirectoryNames = {'RealSHMONO_0Sph','RealSHMONO_.1','RealSHMONO_.2','RealSHMONO_.3', 'RealSHMONO_.4','RealSHMONO_.5', 'RealSHMONO_.6', 'RealSHMONO_.7', 'RealSHMONO_.8', 'RealSHMONO_.9', 'RealSHMONO_1'}

nGrainPacks = numel(grainDirectoryNames);

meanSphericity = zeros(nGrainPacks,1);
meanVolume = zeros(nGrainPacks,1);
meanSurfaceArea =  zeros(nGrainPacks,1);
meanCmean =  zeros(nGrainPacks,1);
meanCgaussian =  zeros(nGrainPacks,1);

sphericity = {}; surfaceArea={}; volume = {};
Cmean={} ; Cgaussian={};

for i = 1:nGrainPacks
i
grainPack = doGrainStuff(grainDirectoryNames{i});
[sphericity{i}, volume{i}, surfaceArea{i}] = grainPack.calculateSphericity();
[Cmean{i} , Cgaussian{i}]  =  grainPack.calculateMeanGrainCurvature();

meanSphericity(i) = mean(sphericity{i});
meanVolume(i) = mean(volume{i});
meanSurfaceArea(i) = mean(surfaceArea{i});
meanCmean(i) = mean(Cmean{i});
meanCgaussian(i) = mean(Cgaussian{i});
end

plot(amp, meanCmean)
