function fixationSalience_and_significanceSMT(FIXATIONFILE,imageX,imageY)
% created by Seth Koenig 6/21/2012 modified for SMT task on 10/25/2013

% function determines normalized salience and image intensity values at 
% fixation locations. The function also calculates if fixations occur at these values 
% at a probability greater than expected by chance using a z-test. 

% Inputs:
%   FIXATIONFILE: Fixations extracted from cortex e.g. MP120606_1-fixation.mat
%   imset: image set number
%   imageX & imageY: X and Y dimensions of the set's images,respstivley

%Outputs:
%   A .mat file named [FIXATIONFILE(1:end-14) '-FixationStatistics']
%   containg fixations statistics. See variable statvariablenames for
%   detailed explanation of variables in .mat file.

if nargin < 2
    error(['Not enough inputs: function requires ViewingBehaviorFile,'...
        'image directory, imageX, and imageY.'])
elseif nargin < 3
    imageX = 800;
    imageY = 600;
end

load(FIXATIONFILE,'fixationstats','images');

maxfixs = 0;
numtrials = length(fixationstats);
shuffunshuff = cell(1,2); %1 for shuffled 2 for unshuffled data
for s = 1:length(shuffunshuff)
    corrs = cell(numtrials,2);
    for cndlop=1:length(fixationstats)
        
        
        if numtrials > 90
            if images(cndlop) < 10
                imgfile = ['00' num2str(images(cndlop)) '.bmp'];
                load(['00' num2str(images(cndlop)) '-saliencemap.mat'],'fullmap');
            elseif images(cndlop) < 100
                imgfile = ['0' num2str(images(cndlop)) '.bmp'];
                load(['0' num2str(images(cndlop)) '-saliencemap.mat'],'fullmap');
            else
                imgfile = [num2str(images(cndlop)) '.bmp'];
                load([num2str(images(cndlop)) '-saliencemap.mat'],'fullmap');
            end
        else
            if images(cndlop) < 10
                imgfile = ['0' num2str(images(cndlop)) '.bmp'];
                load(['0' num2str(images(cndlop)) '-saliencemap.mat'],'fullmap');
            else
                imgfile = [num2str(images(cndlop)) '.bmp'];
                load([num2str(images(cndlop)) '-saliencemap.mat'],'fullmap');
            end
        end
        img = imread(imgfile);
        if size(img,3) ~= 1;
             img = double(rgb2gray(img));
        else
            img = double(img);
        end
        img= img/max(max(img));
    
        saliencemap = fullmap;
        imageX = size(fullmap,2);
        imageY = size(fullmap,1);
        fixations = fixationstats{cndlop}.fixations;
        fixations(:,1) = []; %remove fixation often still on central crosshiar
        if ~isempty(fixations)

            numfixs = size(fixations,2);
            maxfixs = max(maxfixs,numfixs);
            corrs{cndlop,1} = NaN(1,numfixs);
            corrs{cndlop,2} = NaN(1,numfixs);
            
            for i = 1:numfixs
                if s == 1; %shuffled points
                    spot = [ceil(imageX*rand) ceil(imageY*rand)]; %fake x,y data
                else %unshuffled
                    spot = ceil(fixations(:,i));
                    %spot(2) = imageY-spot(2); %import code already flips
                    spot(spot < 1) = 1;
                    spot(1,spot(1) > imageX) = imageX;
                    spot(2,spot(2) > imageY) = imageY;
                end
                corrs{cndlop,1}(i) = saliencemap(spot(2),spot(1));
                corrs{cndlop,2}(i) = img(spot(2),spot(1));
            end
        end
    end
    shuffunshuff{s} = corrs;
end

for s = 1:length(shuffunshuff)
    combineddata = NaN(numtrials,2,maxfixs);
    for i = 1:size(combineddata,1)
        for ii = 1:size(combineddata,2)
            if ~isempty(shuffunshuff{s}{i,ii})
                combineddata(i,ii,1:length(shuffunshuff{s}{i,ii}))= ...
                    shuffunshuff{s}{i,ii};
            end
        end
    end
    
    meanvals = NaN(2,maxfixs);
    stdvals = NaN(2,maxfixs);
    numvals = NaN(2,maxfixs);
    for i = 1:size(meanvals,1)
        for ii = 1:size(meanvals,2);
            if  sum(~isnan(combineddata(:,i,ii))) > 5;
                meanvals(i,ii) = nanmean(combineddata(:,i,ii));
                stdvals(i,ii) = nanstd(combineddata(:,i,ii));
                numvals(i,ii) = sum(~isnan(combineddata(:,i,ii)));
            end
        end
    end
    
    shuffunshuffdata{s} = {meanvals stdvals combineddata};
end

% z-test of means agains random distributions assuming mean is larger
zp = NaN(size(meanvals)); %p-values
cI = NaN(size(meanvals)); %top confidence interval value, lowest is typticall 0/-Inf
for i = 1:size(meanvals,1)
    for ii = 1:size(meanvals,2)
        shuffleddata = shuffunshuffdata{1}{3}(:,i,ii);
        shuffleddata(isnan(shuffleddata)) = [];
        [~,p,ci] = ztest(shuffleddata,meanvals(i,ii),std(shuffleddata),...
            0.05,'left');
        zp(i,ii) = p;
        cI(i,ii) = ci(2);
    end
end

meanvals(:,isnan(meanvals(1,:)))= [];
stdvals(:,isnan(stdvals(1,:))) = [];
numvals(:,isnan(numvals(1,:))) = [];
zp(:,isnan(zp(1,:))) = [];
cI(:,isnan(cI(1,:))) = [];

statistics.meanvalues = meanvals;
statistics.stdvalues = stdvals;
statistics.numbervalues = numvals;
statistics.pvalues = zp;
statistics.confidenceintervals = cI;


statvariablenames = {
    'shuffunshuffdata: cell arrray containing shuffled {1} and unshuffled {2} values.';
    'The cell arrays contain meanvalues, stdvalues and all combined values for';
    'salience and image intensity at fixations, respectively.';
    'Each of these cells are arranged by row and column in the following';
    'manner:  Rows are arranged as Salience then image intensity; and columns by fixation number';
    '';
    'statistics: structure with results from a z-test to determine if fixations';
    'occur at salience and image intesntisy valeus at rates';...
    'higher than what would be expected by chance. Random distributions from each';
    'parameter is compared to the mean value at that parameter by fixation';
    'statistics contains means, std, number of fixations, p-values,and confidence intervals';
    'These variables are arranged by row and column. Rows are arranged as';
    'Salience and image intensity. Columns are organized by fixation number';
    };

save([FIXATIONFILE(1:end-13) '-FixationStatistics.mat'],...
    'shuffunshuffdata','statistics','statvariablenames')
end