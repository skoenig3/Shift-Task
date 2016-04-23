%% Runs whole Shift Task (SFT) analysis
% Similar to RunWholeSalienceModel
% sft set 7 doesn't have enough good calibration trials to calibrate properly
%
% search [x] to get to the desired section
% 1. Get the salience maps for multiple task sets
% 1.5 Get Average Saliency Map
% 2. Get Fixations and Saccades from Behavioral Files
% 3. Calculate Salience at Each Fixation
% 4. Calculate Average Salience at Each Fixation Across mutliple data sets
% 5. Extract Viewing Behavior
% 6. Combine Viewing Behavior by Monkey
% 7. N/A
% 8. Run BCRW
% 9. BCRW Centroids
% 10. Data Centroids
%%
%---[1] Get the salience maps for multiple task sets---%
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'}; 
for imset = 1:length(image_sets);
    dirName = [sft_image_dir image_sets{imset}];
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    imageindexes = [];
    for i = 1:length(fileList)
        bmps = strfind(fileList{i},'bmp');
        if ~isempty(bmps)
            if double(fileList{i}(bmps-2)) <= 57 %ascii for number
                imageindexes = [imageindexes i];
            end
        end
    end
    for i = 1:length(imageindexes)
        imagefile = fileList{imageindexes(i)};
        getSalienceMap(imagefile)
    end
end
%%
%---[1.5] Get Average Saliency Map---%%
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'}; 
imageX = 756;
imageY = 378;
avgmap = zeros(imageY,imageX);
gray = NaN(2,1);
for imset = 1:length(image_sets);
    dirName = [sft_image_dir image_sets{imset}];
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    imageindexes = [];
    for i = 1:90
        if i < 10
            load(['0' num2str(i) '-saliencemap'],'fullmap')
%         elseif i < 100
%             load(['0' num2str(i) '-saliencemap'],'fullmap');
        else
            load([num2str(i) '-saliencemap'],'fullmap');
        end
        if size(fullmap,2) == 755 %1 set has images that are 755 instead of 756
            fullmap = [zeros(378,1) fullmap];
            avgmap = avgmap+fullmap;
        else
            avgmap = avgmap+fullmap;
        end
    end
end
imagesc(avgmap)


%%
%---[2] Get Fixations and Saccades from Behavioral Files---%
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'}; 

imageX = 756;
imageY = 378;
CNDFile = [sft_image_dir 'sft_cnd.cnd'];
for imset = 1:length(image_sets);
    dirName = [sft_image_dir image_sets{imset}];
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    eyeindexes = [];
    for i = 1:length(fileList)
        period = strfind(fileList{i},'.');
        if (length(fileList{i}) == period(end)+1) && ...
                (double(fileList{i}(period+1)) <=57)
            eyeindexes = [eyeindexes i];
        end
    end
    for i = 1:length(eyeindexes)
        cortexfile = fileList{eyeindexes(i)};
        ITMFile = [image_sets{imset} '.itm'];
        getSFTeyedat(cortexfile,ITMFile,CNDFile,imageX,imageY)
    end
end
clearvars -except sft_image_dir image_sets
%%
%---[3] Calculate Salience at Each Fixation---%
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'}; 
imageX = 756; imageY = 378;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            eyedatafiles = [eyedatafiles i];
        end
    end
    for eyefile = eyedatafiles;
        fixationSalience_and_significanceSFT(matfiles.mat{eyefile},imageX,imageY)
    end
end
clearvars -except sft_image_dir image_sets
%%
%%---[4] Calculate Average Salience at Each Fixation Across mutliple data sets---%
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets ={'sft22','sft24','sft27','sft5'}; 
tags = {'WR'};
minlen = 100;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    matfiles = what;
    statfiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'FixationStatistics.mat');
        if ~isempty(str)
            statfiles = [statfiles i];
        end
    end
    for stat = statfiles;
        load(matfiles.mat{stat},'statistics')
        minlen = min(minlen,size(statistics.numbervalues,2));
    end
end

alldata = NaN(120*length(image_sets),2,minlen);
allshuffled = NaN(120*length(image_sets),2,minlen);
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    matfiles = what;
    statfiles = zeros(1,length(tags));
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'FixationStatistics.mat'));
            for ii = 1:length(tags);
                if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                    statfiles(ii) = i;
                end
            end
        end
    end
    
    for stat = 1:length(statfiles);
        if statfiles(stat) ~= 0
            i = length(tags)*(SET-1)+stat;
            load(matfiles.mat{statfiles(stat)})
            combineddata = shuffunshuffdata{2}{3};
            combinedshuffled = shuffunshuffdata{1}{3};
            alldata(120*(i-1)+1:(120*(i-1)+1+size(combineddata,1)-1),:,:) = combineddata(:,:,1:minlen);
            allshuffled(120*(i-1)+1:(120*(i-1)+1+size(combineddata,1)-1),:,:) = combinedshuffled(:,:,1:minlen);
        end
    end
end

allmeanvals = NaN(2,minlen);
allstdvals = NaN(2,minlen);
allnumvals = NaN(2,minlen);
for i = 1:size(allmeanvals,1)
    for ii = 1:minlen
        allmeanvals(i,ii) = nanmean(alldata(:,i,ii));
        allstdvals(i,ii) = nanstd(alldata(:,i,ii));
        allnumvals(i,ii) = sum(~isnan(alldata(:,i,ii)));
    end
end

% z-test of means agains random distributions assuming mean is larger
allzp = NaN(size(allmeanvals)); %p-values
allcI = NaN(size(allmeanvals)); %top confidence interval value, lowest is typticall 0/-Inf
for i = 1:size(allmeanvals,1)
    shuffledvals = allshuffled(:,i,:);
    shuffledvals(isnan(shuffledvals)) = [];
    for ii = 1:size(allmeanvals,2)
        [~,p,ci] = ztest(shuffledvals,allmeanvals(i,ii),std(shuffledvals),...
            0.05);
        allzp(i,ii) = p;
        allcI(i,ii) = ci(2);
    end
end

allstatistics.meanvalues = allmeanvals;
allstatistics.stdvalues = allstdvals;
allstatistics.numbervalues = allnumvals;
allstatistics.pvalues = allzp;
allstatistics.confidenceintervals = allcI;


figure
hold on
plot(allcI(1,1)*ones(1,size(allcI,2)),'--k');
errorbar(allmeanvals(1,:),allstdvals(1,:)./sqrt(allnumvals(1,:)),'r');
hold off
legend('Chance Levels','Salience @ Fixations')
xlabel('Fixation Number')
ylabel('Normalized Salience')

figure
hold on
plot(allcI(2,1)*ones(1,size(allcI,2)),'--k');
errorbar(allmeanvals(2,:),allstdvals(2,:)./sqrt(allnumvals(2,:)),'r');
hold off
legend('Chance Levels','Image Intensity @ Fixations')
xlabel('Fixation Number')
ylabel('Normalized Salience')

allstatvariablenames = {
    'alldata: cell arrray containing combined values for salience, salience';
    'contrast, and image intensity at fixations across all monkeys and data files.';
    'Each of these cells are arranged by row,column,and z-column in the following';
    'manner:  Rows are arranged as Salience, salience contrast, and intensity;';
    'Columns indicate if data is these parmaters at the average fixation cooridante';
    '(col 1) or the average of the parameters during a fixatoin; Z-column is';
    'organized by fixation number';
    '';
    'allshuffled: same as alldata but shuffled data';
    '';
    'allstatistics: structure with results from a z-test to determine if fixations';
    'occur at salience, salience contrasts, and image intesntisy valeus at rates';...
    'higher than what would be expected by chance. Random distributions from each';
    'parameter is compared to the mean value at that parameter by fixation';
    'statistics contains means, std, number of fixations, p-values,and confidence intervals';
    'These variables are arranged by row,column,z-column. Rows are arranged as';
    'Salience, salience contrast, and intensity. Columns indicate if data is these';
    'paramaters at the average fixation cooridante (col 1) or the average of the';
    'parameters during a fixatoin. Z-column is organized by fixation number';
    };

save('C:\Users\seth.koenig\Documents\MATLAB\Shift Task\Combined-Wilbur-Salience',...
    'alldata','allshuffled','allstatistics','allstatvariablenames');
clearvars -except sft_image_dir image_sets
%%
%---[5] Extract Viewing Behavior---%
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'}; 
tags = {'WR'};
PLOTOPTIONS = 'none';
imageX = 756;
imageY = 378;
SAMPRATE = 5;
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            eyedatafiles = [eyedatafiles i];
        end
    end
    for eyefile = eyedatafiles;
        getViewingBehaviorSFT(matfiles.mat{eyefile},SAMPRATE,imageX,imageY,PLOTOPTIONS)
    end
end
clearvars -except sft_image_dir image_sets
%%
%---[6] Combine Viewing Behavior by Monkey---%
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'}; 
tags = {'WR'};
imageX = 756; imageY = 378;
medianfix = NaN(length(image_sets),length(tags));
mediansac = NaN(length(image_sets),length(tags));
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    matfiles = what;
    statfiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'ViewingBehavior');
        if ~isempty(str)
            for ii = 1:length(tags);
                strt = strfind(matfiles.mat{i},tags{ii});
                if ~isempty(strt)
                    load(matfiles.mat{i},'avgfixprofile','avgsacprofile');
                    medianfix(SET,ii) = size(avgfixprofile,2);
                    mediansac(SET,ii) = size(avgsacprofile,2);
                end
            end
        end
    end
end
medianfix = round(nanmedian(medianfix));
mediansac = round(nanmedian(mediansac));

allview = cell(1,length(tags));
for i = 1:length(tags)
    allview{i}.densitymap = zeros(imageY,imageX);
    allview{i}.allfixations = [];
    allview{i}.allsaccades = [];
    allview{i}.persistence =[];
    allview{i}.anglebtwfix = [];
    allview{i}.sacangle_2fix = [];
    allview{i}.distanceprofile = [];
    allview{i}.distbtwnfix = [];
    allview{i}.fixduration = [];
    allview{i}.sacangle = [];
    allview{i}.sacdist = [];
    allview{i}.sacamplitude = [];
    allview{i}.sacduration = [];
    allview{i}.timebtwfix = [];
end
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    matfiles = what;
    statfiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'ViewingBehavior');
        if ~isempty(str)
            for ii = 1:length(tags);
                strt = strfind(matfiles.mat{i},tags{ii});
                if ~isempty(strt)
                    load(matfiles.mat{i});
                    if size(allfixations,2) == medianfix(ii);
                        timewarp = 1:size(allfixations,2);
                    else
                        timewarp = round(linspace(1,size(allfixations,2),medianfix(ii)));
                    end
                    distanceprofile.fix = distanceprofile.fix(:,timewarp);
                    persistence.fix = persistence.fix(:,timewarp);
                    persistence.fix = persistence.fix(:,6:end-5);
                    distanceprofile.fix = distanceprofile.fix(:,6:end-5);
                    allview{ii}.allfixations = [allview{ii}.allfixations;...
                        allfixations(:,timewarp,:)];
                    if size(allsaccades,2) == mediansac(ii);
                        timewarp = 1:size(allsaccades,2);
                    else
                        timewarp = round(linspace(1,size(allsaccades,2),mediansac(ii)));
                    end
                    allview{ii}.allsaccades = [allview{ii}.allsaccades;...
                        allsaccades(:,timewarp,:)];
                    distanceprofile.sac = distanceprofile.sac(:,timewarp);
                    distanceprofile.sac = distanceprofile.sac(:,6:end-5);
                    persistence.sac = persistence.sac(:,timewarp);
                    persistence.sac = persistence.sac(:,6:end-5);
                    allview{ii}.persistence = [ allview{ii}.persistence;
                        [persistence.sac persistence.fix]];
                    allview{ii}.anglebtwfix = [allview{ii}.anglebtwfix;anglebtwfix];
                    allview{ii}.sacangle_2fix = [allview{ii}.sacangle_2fix;...
                        sacangle_2fix];
                    allview{ii}.densitymap = allview{ii}.densitymap+densitymap;
                    allview{ii}.distanceprofile = [allview{ii}.distanceprofile;
                        [distanceprofile.sac distanceprofile.fix]];
                    allview{ii}.distbtwnfix = [allview{ii}.distbtwnfix;distbtwnfix];
                    allview{ii}.fixduration = [allview{ii}.fixduration;...
                        fixduration];
                    allview{ii}.sacangle = [allview{ii}.sacangle;sacangle];
                    allview{ii}.sacdist = [allview{ii}.sacdist;sacdist];
                    allview{ii}.sacamplitude = [allview{ii}.sacamplitude;sacamplitude];
                    allview{ii}.sacduration = [allview{ii}.sacduration;sacduration];
                    allview{ii}.timebtwfix = [allview{ii}.timebtwfix;...
                        timebtwfix];
                    allview{ii}.mediansac = mediansac(ii)-10;
                    allview{ii}.medianfix = medianfix(ii)-10;
                end
            end
        end
    end
end

clearvars -except sft_image_dir image_sets allview tags imageX imageY

SAMPRATE = 5;
n = (-180:180)*pi/180;
variables = {'Dist','vel','accel','rot'};
f = fspecial('gaussian',[256,256],24);
graphnum = gcf;
if graphnum == 1
    graphnum = 0;
end
for i = 1:length(tags)
    [allprobanglebtwfix] = hist(allview{i}.anglebtwfix(~isnan(allview{i}.anglebtwfix)),360);
    allprobanglebtwfix = [allprobanglebtwfix(36:-1:1) allprobanglebtwfix allprobanglebtwfix(end:-1:end-36)];
    allprobanglebtwfix = filtfilt(1/6*ones(1,6),1,allprobanglebtwfix);
    allprobanglebtwfix = allprobanglebtwfix(37:end-37);
    allprobanglebtwfix = allprobanglebtwfix/sum(allprobanglebtwfix);
    allprobanglebtwfix = [allprobanglebtwfix allprobanglebtwfix(1)];
    
    hax(1) = figure(graphnum+1);
    polar(n,allprobanglebtwfix)
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
    title(['Distribution of angles between fixations ' tags{i}])
    
    [allprobsacangle] = hist(allview{i}.sacangle(~isnan(allview{i}.sacangle)),360);
    allprobsacangle = [allprobsacangle(36:-1:1) allprobsacangle allprobsacangle(end:-1:end-36)];
    allprobsacangle = filtfilt(1/6*ones(1,6),1,allprobsacangle);
    allprobsacangle = allprobsacangle(37:end-37);
    allprobsacangle = allprobsacangle/sum(allprobsacangle);
    allprobsacangle = [allprobsacangle allprobsacangle(1)];
    
    hax(2) = figure(graphnum+2);
    polar(n,allprobsacangle)
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
    title(['Distribution of angles leaving a fixation ' tags{i}])
    
    %---Stats by fixation---%
    allStatsbyfixation{i}.fixatoinspertrial = sum(~isnan(allview{i}.fixduration),2);
    allStatsbyfixation{i}.meanfixationduration = SAMPRATE*nanmean(allview{i}.fixduration);
    allStatsbyfixation{i}.stdfixationduration = SAMPRATE*nanstd(allview{i}.fixduration);
    allStatsbyfixation{i}.numfix = sum(~isnan(allview{i}.fixduration));
    allStatsbyfixation{i}.meansacdistance = nanmean(allview{i}.sacdist);
    allStatsbyfixation{i}.stdsacdistance = nanstd(allview{i}.sacdist);
    allStatsbyfixation{i}.numsacs = sum(~isnan(allview{i}.sacdist));
    
    hax(3) = figure(graphnum+3);
    fixduration = allview{i}.fixduration(~isnan(allview{i}.fixduration))*SAMPRATE;
    fixduration(fixduration > 500) = [];
    hist(fixduration,95)
    xlabel('Time (ms)')
    title(['Distribution of Fixation Durations ' tags{i}])
    
    hax(4) = figure(graphnum+4);
    hist(allview{i}.distbtwnfix(~isnan(allview{i}.distbtwnfix)),100)
    xlabel('Distance (Pixels)')
    title(['Distribution of Distances between Fixations ' tags{i}])
    
    hax(5) = figure(graphnum+5);
    hist(1000./(allview{i}.timebtwfix(~isnan(allview{i}.timebtwfix))*SAMPRATE),100)
    xlabel('Hz')
    title(['Fixation Rate ' tags{i}])
    
    hax(6) = figure(graphnum+6);
    sacduration = allview{i}.sacduration(~isnan(allview{i}.sacduration))*SAMPRATE;
    sacduration(sacduration > 150) = [];
    hist(sacduration,25)
    xlabel('Time (ms)')
    title(['Distribution of Saccade Durations ' tags{i}])
    
    hax(7) = figure(graphnum+7);
    sacdist = allview{i}.sacdist(~isnan(allview{i}.sacdist));
    sacdist(sacdist > imageY) = [];
    hist(sacdist,100)
    xlabel('Distance (Pixels)')
    subtitle(['Distribution of Saccade Distances (arc length) ' tags{i}])
    
    hax(8) = figure(graphnum+8);
    plot(allview{i}.sacdist(~isnan(allview{i}.sacdist)),...
        allview{i}.sacduration(~isnan(allview{i}.sacduration))*SAMPRATE,...
        '*','markersize',2)
    xlabel('Distance (pixels)')
    ylabel('Duration (ms)')
    xlim([0 1000])
    ylim([0 200])
    title(['Correlation between saccade duration and distance ' tags{i}])
    
    hax(9) = figure(graphnum+9);
    title(['Probability Distribution of Fixations ' tags{i}])
    densitymap = allview{i}.densitymap;
    densitymap = imfilter(densitymap,f);
    densitymap = densitymap./sum(sum(densitymap));
    imagesc(densitymap)
    axis off
    
    hax(10) = figure(graphnum+10);
    hold all
    plot(allStatsbyfixation{i}.fixatoinspertrial)
    errorbar(allStatsbyfixation{i}.meanfixationduration,...
        allStatsbyfixation{i}.stdfixationduration./sqrt(allStatsbyfixation{i}.numfix))
    xl = find(allStatsbyfixation{i}.numfix < 50);
    xlim([0 xl(1)])
    ylim([0 400])
    errorbar(allStatsbyfixation{i}.meansacdistance,...
        allStatsbyfixation{i}.stdsacdistance./sqrt(allStatsbyfixation{i}.numsacs))
    if i == 1;
        ylabel('Number of Fixations, Fixation Duration (ms),Saccade Distance (Pixels)')
        set(get(gca,'YLabel'),'Position',[-7 -0.3 0])
    end
    if i == 2
        legend('# Fixations','Fixation Duration','Saccade Distance','Location',...
            'NorthEastOutside')
    end
    if i > 2
        xlabel('Trial Number, Fixation Number, or Saccade Number')
    end
    title(['Fixation Statistics by Trial Number, Fixation Number, or Saccade Number ' tags{i}])
    
    avgfixation= mean(allview{i}.allfixations,1);
    fixlen = size(avgfixation,2);
    avgfixprofile = zeros(size(avgfixation));
    for ii = 1:size(avgfixation,3);
        avgfixprofile(:,:,ii) = filtfilt(1/3*ones(1,3),1,avgfixation(:,:,ii));
        avgfixprofile(:,:,ii) = avgfixprofile(:,:,ii) - min(avgfixprofile(:,:,ii));
        avgfixprofile(:,:,ii) = avgfixprofile(:,:,ii)/max(avgfixprofile(:,:,ii));
    end
    avgsaccade= mean(allview{i}.allsaccades,1);
    saclen = size(avgsaccade,2);
    avgsacprofile = zeros(size(avgsaccade));
    for ii = 1:size(avgsaccade,3);
        avgsacprofile(:,:,ii) = filtfilt(1/3*ones(1,3),1,avgsaccade(:,:,ii));
        avgsacprofile(:,:,ii) =  avgsacprofile(:,:,ii) - min(avgsacprofile(:,:,ii));
        avgsacprofile(:,:,ii) = avgsacprofile(:,:,ii)/max(avgsacprofile(:,:,ii));
    end
    
    hax(11) = figure(graphnum+11);
    hold all
    h = area(5:fixlen-5,ones(1,fixlen-9));
    set(h,'FaceColor',[.75 .75 .75])
    set(h,'EdgeColor','none')
    for ii =  1:size(avgfixprofile,3);
        plot(avgfixprofile(:,:,ii),'linewidth',2)
    end
    hold off
    xlim([1 fixlen])
    set(gca,'XTick',[])
    set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
    if i == 1
        ylabel('Normalized Value')
        set(get(gca,'YLabel'),'Position',[-2 -0.2 0])
    end
    if i == 2
        legend([{'fixation'} variables],'Location','NorthEastOutside');
    end
    if i > 2
        xlabel('Warped Time')
    end
    title(['Average-Smoothed Fixation Profile by Parameter ' tags{i}])
    
    hax(12) = figure(graphnum+12);
    hold all
    h1 = area(1:5,ones(1,5));
    set(h1,'FaceColor',[.75 .75 .75])
    set(h1,'EdgeColor','none')
    h2 = area(saclen-4:saclen,ones(1,5));
    set(h2,'FaceColor',[.75 .75 .75])
    set(h2,'EdgeColor','none')
    for ii = 1:size(avgsacprofile,3)
        p(ii) = plot(avgsacprofile(:,:,ii),'linewidth',2);
    end
    hold off
    xlim([1 saclen])
    set(gca,'XTick',[])
    set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
    if i == 1
        ylabel('Normalized Value')
        set(get(gca,'YLabel'),'Position',[-2 -0.2 0])
    end
    if i == 2
        legend([h1 p],[{'fixation'} variables],'Location','NorthEastOutside');
    end
    if i > 2
        xlabel('Warped Time')
    end
    title(['Average-Smoothed Saccade Profile by Parameter ' tags{i}])
    
    hax(13) = figure(graphnum+13);
    hold on
    h  = area(allview{i}.mediansac+1:size(allview{i}.persistence,2),...
        ones(1,size(allview{i}.persistence,2)-allview{i}.mediansac));
    set(h,'FaceColor',[.75 .75 .75])
    set(h,'EdgeColor','none')
    p = plot(mean(allview{i}.persistence));
    hold off
    xlim([1 size(allview{i}.persistence,2)])
    if i == 1
        ylabel('Probability of Saccade Angle Changeing > 45 Degrees')
        set(get(gca,'YLabel'),'Position',[-3 -0.3 0])
    end
    if i == 2
        legend([h p],{'fixation','persistence'},'Location','NorthEastOutside');
    end
    if i > 2
        xlabel('Warped Time')
    end
    title(['Probability of Saccade Angle Changeing > 45 Degrees ' tags{i}])
    
    hax(14) = figure(graphnum+14);
    plot(nanmean(allview{i}.distanceprofile))
    xlabel('Warped Time')
    ylabel('Distance (pixels)')
    title(['Saccade and Fixation Distance ' tags{i}])
    
    [allprobangle2fix] = hist(allview{i}.sacangle_2fix(~isnan(allview{i}.sacangle_2fix)),360);
    allprobangle2fix = [allprobangle2fix(36:-1:1) allprobangle2fix allprobangle2fix(end:-1:end-36)];
    allprobangle2fix = filtfilt(1/6*ones(1,6),1,allprobangle2fix);
    allprobangle2fix = allprobangle2fix(37:end-37);
    allprobangle2fix = allprobangle2fix/sum(allprobangle2fix);
    allprobangle2fix = [allprobangle2fix allprobangle2fix(1)];
    
    hax(15) = figure(graphnum+15);
    polar(n,allprobangle2fix)
    ph=findall(gca,'type','text');
    set(ph,'fontweight','bold');
    title(['Distribution of saccade angles entering fixations ' tags{i}])
    
    hax(16) = figure(graphnum+16);
    sacdist = allview{i}.sacamplitude(~isnan(allview{i}.sacamplitude));
    sacdist(sacdist > imageY) = [];
    hist(sacdist,100)
    xlabel('Distance (Pixels)')
    title(['Distribution of Saccade Amplitudes ' tags{i}])
    
end

screen_size = get(0, 'ScreenSize');
figdir = ['C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'...
    'Plots\Viewing Behavior\'];
figuretitles = {
    'Distribution of Angles between fixations';
    'Distribution of Angles Leaving a fixation';
    'Distribution of Fixation Durations';
    'Distribution of Distances Between Fixations';
    'Fixation (Saccade) Rate';
    'Distribution of Saccade Durations-Arc Length';
    'Distribution of Saccade Distances';
    'Correlation between Saccade Duration and Saccade Distance';
    'Fixation and Saccade Statistics by Fixation or Saccade Number';
    '2-D Fixation PDF'
    'Average-Smoothed Fixation Profile';
    'Average-Smoothed Saccade Profile';
    'Persistence Profile';
    'Saccade and Fixation Distance Profiles';
    'Distribution of angles entering a fixation';
    'Distribution of Saccade Amplitudes';
    };

for ff = 1:graphnum+16;
    figure(graphnum+ff)
    set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    saveas(gcf,[figdir figuretitles{ff}])
    print(gcf,'-r300','-djpeg',[figdir figuretitles{ff}])
end

allviewvariables = {
    'tags: subject names';
    '';
    'allview: all combined data by subject';
    'allview.densitymap: positions of fixations';
    'allview.allfixations: fixation profile by parameters warped twice to median length';
    'allview.allsaccades: saccade profile by parameters warped twice to median length';
    'allview.persistence: persistence/probability of eye movement changing >45 angle double warped';
    'allview.anglebtwfix: angles between fixations by fixation';
    'allview.sacangle_2fix: saccade angles entering a fixation';
    'allview.distanceprofile: velocity of eye movements for saccade + subsequent fixation';
    'allview.distbtwnfix: distance between fixations by fixation number';
    'allview.fixduration: fixation duration by fixation number';
    'allview.sacangle: angle of saccade leaving a fixation by saccade number';
    'allview.sacdist: distance of saccade by saccade number (arc length)';
    'allview.sacamplitude: saccade amplitude from end to end';
    'allview.timebtwfix: time between fixations by fixation number';
    };

save(['C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'...
    'CombinedViewingBehavior-Wilbur.mat'],'tags','allview',...
    'allStatsbyfixation','allviewvariables');
%% ---[8] Run BCRW ---%
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'}; 
tags = {'WR'};
Combinedbehaviorfile = ['C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'...
    'CombinedViewingBehavior-Wilbur.mat'];
load(Combinedbehaviorfile,'allview')
imageX = 756;
imageY = 378;
plotoptions.runs = 'none';
plotoptions.probdens = 'none';
plotoptions.type = 'sal';
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    
    matfiles = what;
    saliencemapfiles = [NaN;NaN];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'saliencemap');
        if ~isempty(str)
            dash = strfind(matfiles.mat{i},'-');
            saliencemapfiles = [saliencemapfiles [i;str2num(matfiles.mat{i}(1:dash(1)-1))]];
        end
    end
    saliencemapfiles(:,1) = [];
    [~,si] = sort(saliencemapfiles(2,:));
    saliencemapfiles = si;
    
    for i = 1:length(saliencemapfiles)
        for t = 1:length(tags)
            disp(['Running ' tags{t} ' on image #' num2str(i) ' from ' image_sets{SET}])
            run_BCRWCF_CentralBias(allview{t},matfiles.mat{saliencemapfiles(i)},tags{t},imageX,imageY,plotoptions)
        end
    end
end
%% [9] BCRW Centroids
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'}; 
tags = {'WR'};
imageX = 756;
imageY = 378;
fixationpdfs = cell(1,9);
centroids = cell(1,9);
centerpoint = [];
for i = 1:9;
    fixationpdfs{i} = zeros(imageY,imageX);
    centerpoint = [centerpoint i*ones(1,100)];
    centroids{i} = NaN(90*length(image_sets),2);
end
img = NaN(9,90*length(image_sets));
center = NaN(9,90*length(image_sets));

imgcount = ones(1,9);
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    disp(SETNUM)
    for i = 1:90;
        for t = 1:length(tags);
            if i < 10
                load([tags{t} '-0' num2str(i) '-CB.mat'],'fixationtimes');
            else
                load([tags{t} '-' num2str(i) '-CB.mat'],'fixationtimes');
            end
            tempfixationpdf = cell(1,9);
            for ti = 1:9
                tempfixationpdf{ti} = zeros(imageY,imageX);
            end
            for run = 1:size(fixationtimes,1);
                tind = find(fixationtimes(run,:,1) > 0);
                for ii = 1:length(tind)
                    x = fixationtimes(run,tind(ii),1);
                    y = fixationtimes(run,tind(ii),2);
                    tempfixationpdf{centerpoint(run)}(y,x) = tempfixationpdf{centerpoint(run)}(y,x) +1;
                end
            end
        end
        for ti = 1:9
            fixationpdfs{ti} = fixationpdfs{ti}+tempfixationpdf{ti};
            [xcm ycm] = centroid(tempfixationpdf{ti});
            centroids{ti}(imgcount(ti),:) = [xcm ycm];
            img(ti,imgcount(ti)) = SET;
            center(ti,imgcount(ti)) = ti;
            imgcount(ti) = imgcount(ti)+1;
        end
    end
end

f = fspecial('gaussian',[256,256],24);
startpos =  [10     189;
    10      368;
    378     10;
    378     189;
    378     368;
    746     10;
    746     189;
    746     368;
    10      10];


order = [2 5 8 1 4 7 9 3 6]; %rearrange for plotting fixation pdfs
figure
for a = 1:9
    i = order(a);
    subplot(3,3,a)
    hold on
    imagesc(imfilter(fixationpdfs{i},f))
    plot(startpos(i,1),startpos(i,2),'+m','markersize',10)
    hold off
    axis off
    axis tight
end

runorder = [4 5 8 7 6 3 9 1 2]; %rearrange for plotting bar graphs
centroids = centroids(runorder);
stds = NaN(9,2);
means = NaN(9,2);
allcentroids = [];
centers = [];
images =[]; %image set
subjects = [];
for i = 1:9
    centroids{i}(:,1) = (centroids{i}(:,1)-imageX/2)/24;
    centroids{i}(:,2) = (centroids{i}(:,2)-imageY/2)/24;
    means(i,:) = nanmean(centroids{i});
    stds(i,:) = nanstd(centroids{i});
    allcentroids = [allcentroids; centroids{i}];
    images = [images;img(i,:)'];
    centers = [centers;center(i,:)'];
end

direction = {'Center','N','NE','E','SE','S','SW','W','NW'};
figure
bar(means);
set(gca,'Xtick',1:9)
set(gca,'Xticklabel',direction)
legend('Horizontal','Vertical')
%%
[P,ANOVATAB,STATS] = anovan(allcentroids(:,1),images)
multcompare(STATS);

[P,ANOVATAB,STATS] = anova1(allcentroids(:,1),centers)
multcompare(STATS);
set(gca,'Yticklabel',direction(end:-1:1))

[P,ANOVATAB,STATS] = anovan(allcentroids(:,2),images)
multcompare(STATS);

[P,ANOVATAB,STATS] = anova1(allcentroids(:,2),centers)
multcompare(STATS);
set(gca,'Yticklabel',direction(end:-1:1))
%% [10] Data Centroids
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'};  
tags = {'WR'};
imageX = 756;
imageY = 378;

startpos = [
    378   189;
    378   368;
    746   368;
    746   189;
    746    10;
    378    10;
    10    10;
    10   189;
    10   368];

fixationpdfs = cell(1,9);
centroids = cell(1,9);
image_position = cell(1,9);
centerpoint = [];
center = [];
for i = 1:9;
    fixationpdfs{i} = zeros(imageY,imageX);
    centroids{i} = [];
    center{i} = [];
end

for SET =1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    disp(SETNUM)
    
    ITMFile = [SETNUM '.itm'];
    itmfil=[];
    h =fopen(ITMFile);
    tline = 1;
    while tline > 0
        tline = fgetl(h);
        if ~isempty(itmfil)
            if length(tline)>size(itmfil,2)
                tline=tline(1:size(itmfil,2));
            end
        end
        tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
        if ischar(tline)
            itmfil=[itmfil; tline];
        else
            break
        end
    end
    fclose(h);
    
    imagepos = [];
    for img = 227:size(itmfil,1);
        temp = textscan(itmfil(img,:),'%d');
        imagepos = [imagepos temp{1}(3)];
    end
    
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            eyedatafiles = [eyedatafiles i];
        end
    end
    for eyefile = eyedatafiles;
        load(matfiles.mat{eyefile})
        for i = 1:length(fixationstats);
            crossposs = round(initialcross(i,:));
            startind = find(crossposs(1) == startpos(:,1) & crossposs(2) == startpos(:,2));
            
            tempfixationpdf = zeros(imageY,imageX);
            fixations = fixationstats{i}.fixations;
            for ii = 2:size(fixations,2)
                x = round(fixations(1,ii));
                y = round(fixations(2,ii));
                x(x < 1) = 1; x(x > imageX) = imageX;
                y(y < 1) = 1; y(y > imageY) = imageY;
                tempfixationpdf(y,x) =  tempfixationpdf(y,x)+1;
            end
            fixationpdfs{startind} = fixationpdfs{startind}+tempfixationpdf;
            [xcm ycm] = centroid(tempfixationpdf);
            centroids{startind} = [centroids{startind}; [xcm ycm]];
            center{startind} = [center{startind}; startind];
            image_position{startind} = [image_position{startind} imagepos(images(i))];
        end
    end
end

f = fspecial('gaussian',[256,256],24);
order = [9 2 3 8 1 4 7 6 5];
runorder = [4 5 8 7 6 3 9 1 2];
figure
for a = 1:9
    i = order(a);
    subplot(3,3,a)
    hold on
    imagesc(imfilter(fixationpdfs{i},f))
    plot(startpos(i,1),startpos(i,2),'+m','markersize',10)
    hold off
    axis off
    axis tight
end

centroids = centroids(runorder);
stds = NaN(9,2);
means = NaN(9,2);
allcentroids = [];
centers = [];
images =[];
subjects = [];
all_image_pos = [];
for i = 1:9
    centroids{i}(:,1) = (centroids{i}(:,1)-imageX/2)/24;
    centroids{i}(:,2) = (centroids{i}(:,2)-imageY/2)/24;
    means(i,:) = nanmean(centroids{i});
    stds(i,:) = nanstd(centroids{i});
    allcentroids = [allcentroids; centroids{i}];
    centers = [centers;center{i}];
    all_image_pos = [all_image_pos image_position{i}];
end

direction = {'Center','N','NE','E','SE','S','SW','W','NW'};
figure
bar(means);
set(gca,'Xtick',1:9)
set(gca,'Xticklabel',direction)
legend('Horizontal','Vertical')
ylabel('Dva')
%%
% Test if initial fixation cross hair position change centroids
[P,ANOVATAB,STATS] = anova1(allcentroids(:,1),centers);
multcompare(STATS);
set(gca,'Yticklabel',direction(end:-1:1))

[P,ANOVATAB,STATS] = anova1(allcentroids(:,2),centers);
multcompare(STATS);
set(gca,'Yticklabel',direction(end:-1:1))

% Test if image position change centroids
[P,ANOVATAB,STATS] = anova1(allcentroids(:,1),all_image_pos);
multcompare(STATS);

[P,ANOVATAB,STATS] = anova1(allcentroids(:,2),all_image_pos);
multcompare(STATS);

%Test if image position and crosshair position effect data
[P,ANOVATAB,STATS] = anovan(allcentroids(:,1),{all_image_pos',centers'},'model','interaction');

[P,ANOVATAB,STATS] = anovan(allcentroids(:,2),{all_image_pos',centers'},'model','interaction');