function getSFTeyedat(cortexfile,ITMFile,CNDFile,imageX,imageY)
% extracts raw Shift task (SFT) eye data from a cortex files, calibrates
% the eye data, and then uses Cluster Fix to extract fixations and saccades

%ImageX: number of horizontal pixels in image
%ImageY: number of vertical pixels in image
%Only necessary for removing eye tracking data outside of the image


%---read item file data---%
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

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

%---Extract and organize raw cortex data---%
[time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(cortexfile);

numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) < 2500)) ~=0
        if size(find(event_arr(:,rptlop) == 200)) ~=0
            perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
            perendind = find(event_arr(:,rptlop) == 24);
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <= 5000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop)-100;
            endtimdum = time_arr(perendind,rptlop);
            if event_arr(cndnumind,rptlop) <= 1100
                valrptcnt = valrptcnt + 1;
                clrchgind(valrptcnt)=rptlop;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).begpos = 1;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
            end
        end
    end
end

samprate=5;

clear cnd
numrpt = size(per,2);
for rptlop = 1:numrpt
    cnd(rptlop)=per(rptlop).cnd;
end

evnnmb=2:2:size(eog_arr,1);
oddnmb=1:2:size(eog_arr,1);

clear x y
cndlst=unique(cnd);
control = NaN(length(cndlst),2);
for k = 1:length(cndlst)
    C = textscan(itmfil(7+(k-1)*2,:), '%f');
    control(k,:) = C{1}(4:5)';
end

%---For Calibration with Eye tracking data with cp2tform---%
% Create structures x and y of the corresponding average eye data for each trial
% instance (l) of each condition (k)
for k=1:length(cndlst)
    cndind=find(cnd==cndlst(k));
    allind=clrchgind(cndind);
    for l=1:length(allind)
        x{k}(l)=mean(eog_arr(intersect(floor(((per(cndind(l)).begsmpind-1000)/samprate)*2):(floor((per(cndind(l)).endsmpind-1000)/samprate))*2,oddnmb),allind(l)));
        y{k}(l)=mean(eog_arr(intersect(floor(((per(cndind(l)).begsmpind-1000)/samprate)*2):(floor((per(cndind(l)).endsmpind-1000)/samprate))*2,evnnmb),allind(l)));
    end
end

clear meanx meany


for k=1:numel(x)
    xss = x{k};
    low = mean(xss)-std(xss);
    high = mean(xss)+std(xss);
    xss(xss < low) = [];
    xss(xss > high) = [];  
    meanx(k)=median(xss);
end

for k=1:numel(y)
    yss = y{k};
    low = mean(yss)-std(yss);
    high = mean(yss)+std(yss);
    yss(yss < low) = [];
    yss(yss > high) = [];
    meany(k)=median(yss);
end

input = zeros(length(cndlst), 2);

for i = 1:length(cndlst)
    input(i, 1) = meanx(i);
    input(i, 2) = meany(i);
end
tform = cp2tform(control, input,'polynomial',4);
tform.forward_fcn = tform.inverse_fcn;

figure
hold on
for i = 1:length(control);
    plot(control(i,1),control(i,2),'r+')
    [x,y] = tformfwd(tform,meanx(i),meany(i));
    plot(x,y,'*b')
end

numrpt = size(event_arr,2);
valrptcnt = 0;
clear per vpcind
new_eog_arr=[];
for rptlop = 1:numrpt
    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) >= 1010)) ~=0
        if size(find(event_arr(:,rptlop) == 200)) ~=0
            perbegind = find(event_arr(:,rptlop) == 23,1,'first');
            perendind = find(event_arr(:,rptlop) == 24,1,'first');
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2500);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop);
            endtimdum = time_arr(perendind,rptlop);
            if event_arr(cndnumind,rptlop) > 1100
                valrptcnt = valrptcnt + 1;
                vpcind(valrptcnt)=rptlop;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).begpos = 1;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
            end
        end
    end
end


%---get eye data for only when fixation cross or picture is displayed---%
eyedat = [];
for trlop=1:size(per,2)
    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
    picstart=per(trlop).alltim(find(per(trlop).allval==23,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture start time relative to eye scan start
    picend=per(trlop).alltim(find(per(trlop).allval==24,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture end time relative to eye scan start
    
    if picend/5<=length(horeog)
        
        eyedat{trlop}(1,:) = (horeog(round(picstart/5):floor(picend/5)));
        eyedat{trlop}(2,:) = (vrteog(round(picstart/5):floor(picend/5)));
        
    else
        eyedat{trlop}(1,:)=nan;
        eyedat{trlop}(2,:)=nan;
    end
    
end

cnd=[];
numrpt = size(per,2);
for rptlop = 1:numrpt
    cnd(rptlop)=per(rptlop).cnd;
end

%---Recalibrate and automatically scale eye data---%
for eye = 1:length(eyedat)
    x = eyedat{eye}(1,:);
    y = eyedat{eye}(2,:);
    [x y] = tformfwd(tform,x,y);
    eyedat{eye} = [x;y];
end


%---Post processing that I used to organize my data probably not necessary---%
images = NaN(length(eyedat),1);
initialcross = NaN(length(eyedat),2);
for i = 1:length(eyedat);
    cndline = cnd(i)-1000+1;
    C = textscan(cndfil(cndline,:),'%d');
    itm = C{1}(5)+6;
    cross = C{1}(4)+6;
    C = textscan(itmfil(cross,:),'%f');
    initialcross(i,:) = C{1}(4:5);
    imgind = strfind(itmfil(itm,:),'\');
    img = itmfil(itm,imgind(end)+1:imgind(end)+3);
    images(i) = str2double(img);
    C = textscan(itmfil(itm,:),'%f');
    shift = 25.2*double(C{1}(3));
    x = 25.2*eyedat{i}(1,:)+imageX/2-shift; 
    y = -25.2*eyedat{i}(2,:)+imageY/2; %flip to image coordinates and offset correct
    eyedat{i} = [x;y];
    initialcross(i,:) = [25.2*initialcross(i,1)+imageX/2-shift,-25.2*initialcross(i,2)+imageY/2];
end

repeat_ind = []; 
for img = 1:max(images);
    ind = find(images == img);
    if length(ind) == 2
        repeat_ind = [repeat_ind ind(2)];
    end
end
%remove the 2nd presentation of each image
eyedat(repeat_ind) = [];
initialcross(repeat_ind,:) = [];
images(repeat_ind) = [];

%---Removes eye data outside of image---%
for i = 1:size(eyedat,2);
    x = eyedat{i}(1,:);
    y = eyedat{i}(2,:);
    badx = find(x < -50 | x > imageX+50); %~2 dva leave margin of error
    x(badx) = []; y(badx) = [];
    bady = find(y < -50 | y > imageY+50); %~2 dva margin of error
    x(bady) = []; y(bady) = [];
    eyedat{i} = [x;y];
end

[fixationstats] = ClusterFixation_Final(eyedat);
save([cortexfile(1:end-2) '_' cortexfile(end) '-fixation.mat'],...
    'fixationstats','initialcross','images')
end