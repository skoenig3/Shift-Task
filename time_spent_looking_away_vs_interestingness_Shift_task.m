clar
Matlab_dir = 'C:\Users\seth.koenig\Documents\MATLAB\';
IOR_dir = [Matlab_dir 'IOR\'];

% sobel filters detect edges!
sobelx = [1     2   1;
    0     0   0;
    -1    -2  -1;];

sobely = [1     2   1;
    0     0   0;
    -1    -2  -1;];

tasks = {'Shift'};

image_dirs = {'\Shift Task\'};

sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft6','sft8','sft9'};

sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\'; 
image_sets = {'sft22','sft24','sft27','sft5'}; 

entropyvalues = cell(1,length(tasks));
edgevalues = cell(1,length(tasks));
lookingtime = cell(1,length(tasks));
rcount = ones(2,length(tasks));
for task = 1:length(image_dirs);
    
    entropyvalues{task} = NaN(2,1750);
    edgevalues{task} = NaN(2,1750);
    lookingtime{task} = NaN(2,1750);
    
    cd([IOR_dir 'Eye Data\' tasks{task}]);
    disp(tasks{task})
    a = what;
    for m = 1:length(a.mat)
        if isempty(strfind(a.mat{m},'SalienceIOR'))
            load(a.mat{m});
            if strcmpi(tasks{task},'SCM')
                novelconditions = 1:2:length(fixationstats);
                repeatconditions = 2:2:length(fixationstats);
                images = 1:36;
                validconditions = 1010;
            elseif strcmpi(tasks{task},'LIST')
                novelconditions = [];
                for im = 1:length(images)
                    ind = find(images == images(im));
                    if length(ind) == 1 || ind(1) == im
                        novelconditions = [novelconditions im];
                    end
                end
                repeatconditions = 1:length(images);
                [~,ia,~] = intersect(repeatconditions,novelconditions);
                repeatconditions(ia) = [];
                validconditions = 1010;
            elseif strcmpi(tasks{task},'SMT')
                novelconditions = [];
                for im = 1:ceil(length(images)/2);
                    ind = find(images == im);
                    if ~isempty(ind)
                        novelconditions = [novelconditions ind(1)];
                    end
                end
                repeatconditions = 1:length(images);
                [~,ia,~] = intersect(repeatconditions,novelconditions);
                repeatconditions(ia) = [];
                validconditions = 1080;
            elseif strcmpi(tasks{task},'VPLT')
                novelconditions = 1:length(fixationstats)/2;
                repeatconditions = novelconditions(end)+1:length(fixationstats);
                validconditions = 1010;
            end
            
            cfile = a.mat{m};
            cfile = [cfile(1:8) '.' cfile(10)];
            
            if strcmpi(cfile(1:2),'PW')
                cortexfile = ['R:\Buffalo Lab\Cortex Data\Vivian\' cfile];
            elseif strcmpi(cfile(1:2),'TT')
                cortexfile = ['R:\Buffalo Lab\Cortex Data\Timmy\' cfile];
            elseif strcmpi(cfile(1:2),'MP')
                cortexfile = ['R:\Buffalo Lab\Cortex Data\Peepers\' cfile];
            elseif strcmpi(cfile(1:2),'JN')
                cortexfile = ['R:\Buffalo Lab\Cortex Data\Guiseppe\' cfile];
            elseif strcmpi(cfile(1:2),'WR')
                cortexfile = ['R:\Buffalo Lab\Cortex Data\Wilbur\' cfile];
            elseif strcmpi(cfile(1:2),'RR')
                cortexfile = ['R:\Buffalo Lab\Cortex Data\Red\' cfile];
            elseif strcmpi(cfile(1:2),'BZ')
                cortexfile = ['R:\Buffalo Lab\Cortex Data\Bizzie\' cfile];
            elseif strcmpi(cfile(1:2),'IW')
                cortexfile = ['R:\Buffalo Lab\Cortex Data\Irwin\' cfile];
            else
                disp('could not find cortexfile. Check source folder and file name')
            end
            
            [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(cortexfile);
            
            cnd = [];
            new_eog_arr=[];
            per = [];
            numrpt = size(event_arr,2);
            valrptcnt = 0;
            for rptlop = 1:numrpt
                if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) >= validconditions)) ~=0
                    if size(find(event_arr(:,rptlop) == 200)) ~=0
                        perbegind = find(event_arr(:,rptlop) == 23,1,'first');
                        perendind = find(event_arr(:,rptlop) == 24,1,'first');
                        cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                        blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                        begtimdum = time_arr(perbegind,rptlop);
                        endtimdum = time_arr(perendind,rptlop);
                        if endtimdum > begtimdum
                            valrptcnt = valrptcnt + 1;
                            vpcind(valrptcnt)=rptlop;
                            per(valrptcnt).begsmpind = begtimdum;
                            per(valrptcnt).endsmpind = endtimdum;
                            per(valrptcnt).begpos = 1;
                            per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                            per(valrptcnt).blk = event_arr(blknumind,rptlop);
                            per(valrptcnt).allval = event_arr(:,rptlop);
                            per(valrptcnt).alltim = time_arr(:,rptlop);
                            cnd = [cnd event_arr(cndnumind,rptlop)];
                            new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                        end
                    end
                end
            end
            eyedat = [];
            for trlop=1:size(per,2)
                trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
                horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
                vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
                picstart=per(trlop).alltim(find(per(trlop).allval==23,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture start time relative to eye scan start
                picend=per(trlop).alltim(find(per(trlop).allval==24,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture end time relative to eye scan start
                
                if picend > length(horeog)*5
                    picend =  length(horeog)*5;
                end
                if picend < picstart
                    eyedat{trlop}(1,:) = NaN;
                    eyedat{trlop}(2,:) = NaN;
                else
                    eyedat{trlop}(1,:) = (horeog(ceil(picstart/5):floor(picend/5)));
                    eyedat{trlop}(2,:) = (horeog(ceil(picstart/5):floor(picend/5)));
                    
                end
            end
            
            img_dir = [Matlab_dir image_dirs{task}];
            if strcmpi(tasks{task},'SCM')
                for cnd = 1:length(novelconditions)
                    img = imread([img_dir imgset '\' num2str(cnd) '.bmp']);
                    img = rgb2gray(img);
                    ent = entropy(img);%pixel intesnity entropy
                    xedges = imfilter(img,sobelx);
                    yedges = imfilter(img,sobely);
                    edg = mean2(xedges+yedges); %edgineess
                    
                    entropyvalues{task}(1,rcount(1,task)) = ent;
                    edgevalues{task}(1,rcount(1,task))=edg;
                    lookingtime{task}(1,rcount(1,task))=5*size(eyedat{novelconditions(cnd)},2);
                    rcount(1,task) = rcount(1,task)+1;
                end
                for cnd = 1:length(repeatconditions);
                    img = imread([img_dir imgset '\' num2str(cnd) '.bmp']);
                    img = rgb2gray(img);
                    ent = entropy(img);%pixel intesnity entropy
                    xedges = imfilter(img,sobelx);
                    yedges = imfilter(img,sobely);
                    edg = mean2(xedges+yedges); %edgineess
                    
                    entropyvalues{task}(2,rcount(2,task)) =ent;
                    edgevalues{task}(2,rcount(2,task))=edg;
                    lookingtime{task}(2,rcount(2,task))=5*size(eyedat{repeatconditions(cnd)},2);
                    rcount(2,task) = rcount(2,task)+1;
                end
            elseif strcmpi(tasks{task},'LIST')
                for imnum = 1:length(images)
                    ind = find(images == imnum);
                    if ~isempty(ind);
                        if imgset < 10
                            img = imread([img_dir 'List0' num2str(imgset) '\' num2str(imnum) '.bmp']);
                        else
                            img = imread([img_dir 'List' num2str(imgset) '\' num2str(imnum) '.bmp']);
                            
                        end
                        
                        if size(img,3) ~=1
                            img = rgb2gray(img);
                        end
                        ent = entropy(img);%pixel intesnity entropy
                        xedges = imfilter(img,sobelx);
                        yedges = imfilter(img,sobely);
                        edg = mean2(xedges+yedges); %edgineess
                        
                        entropyvalues{task}(1,rcount(1,task)) = ent;
                        edgevalues{task}(1,rcount(1,task))=edg;
                        lookingtime{task}(1,rcount(1,task))=5*size(eyedat{ind(1)},2);
                        rcount(1,task) = rcount(1,task)+1;
                        if length(ind) == 2;
                            entropyvalues{task}(2,rcount(1,task)) =ent;
                            edgevalues{task}(2,rcount(2,task))=edg;
                            lookingtime{task}(2,rcount(2,task))=5*size(eyedat{ind(2)},2);
                            rcount(2,task) = rcount(2,task)+1;
                        end
                    end
                end
            elseif strcmpi(tasks{task},'VPLT')
                for imnum = 1:size(imageprez,2);
                    img = imread([img_dir 'SET' num2str(imgset) '\'  num2str(imnum) '.bmp']);
                    img = rgb2gray(img);
                    ent = entropy(img);%pixel intesnity entropy
                    xedges = imfilter(img,sobelx);
                    yedges = imfilter(img,sobely);
                    edg = mean2(xedges+yedges); %edgineess
                    
                    if ~isnan(imageprez(1,imnum));
                        entropyvalues{task}(1,rcount(1,task)) = ent;
                        edgevalues{task}(1,rcount(1,task))=edg;
                        lookingtime{task}(1,rcount(1,task))=5*size(eyedat{imageprez(1,imnum)},2);
                        rcount(1,task) = rcount(1,task)+1;
                    end
                    if ~isnan(imageprez(2,imnum));
                        entropyvalues{task}(2,rcount(1,task)) = ent;
                        edgevalues{task}(2,rcount(2,task))=edg;
                        lookingtime{task}(2,rcount(2,task))=5*size(eyedat{imageprez(2,imnum)},2);
                        rcount(2,task) = rcount(2,task)+1;
                    end
                end
            end
        end
    end
end
%%
for task = 1:length(tasks)
    figure
    
    lt = lookingtime{task}(1,:);
    et = entropyvalues{task}(1,:);
    ed = edgevalues{task}(1,:);
    
    et(lt < 100) = [];
%     if task == 3
%         lt(end) = [];
%     end
    ed(lt < 100) = [];
    lt(lt < 100) = [];
    
    subplot(2,2,1)
    plot(lt,et,'.')
    xlabel('Image Presentation Times (ms)')
    ylabel('Image Entropy (bits)')
    title('First Presentation')
    
    subplot(2,2,2)
    plot(lt,ed,'.')
    xlabel('Image Presentation Times (ms)')
    ylabel('% Edginess')
    title('First Presentation')
    
    lt = lookingtime{task}(2,:);
    et = entropyvalues{task}(2,:);
    ed = edgevalues{task}(2,:);
    
    et(lt < 100) = [];
%     if task == 3
%         lt(end) = [];
%     end
    ed(lt < 100) = [];
    lt(lt < 100) = [];
    
    subplot(2,2,3)
    plot(lt,et,'.')
    xlabel('Image Presentation Times (ms)')
    ylabel('Image Entropy (bits)')
    title('Second Presentation')
    
    subplot(2,2,4)
    plot(lt,ed,'.')
    xlabel('Image Presentation Times (ms)')
    ylabel('% Edginess')
    title('Second Presentation')
    
    subtitle(tasks{task})
end