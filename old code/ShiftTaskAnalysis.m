clear,clc
itemFileName = 'sft22.itm';
dataFileName = 'WR130412.2';

% load files and analyze

itmfil=[];
[fid,message]=fopen(itemFileName, 'r');
if fid<0
    disp(message);
else
    while 1
        tline = fgetl(fid);
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
end
fclose(fid);

cndfil=[];
[fid,message]=fopen('SFT_CND.cnd', 'r');
if fid<0
    disp(message);
else
    while 1
        tline = fgetl(fid);
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
end
fclose(fid);

[time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(dataFileName);

numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) < 2000)) ~=0
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
for k=1:length(x)
    meanx(k)=mean(x{k});
end
for k=1:length(y)
    meany(k)=mean(y{k});
end

input = zeros(length(cndlst), 2);

for i = 1:length(cndlst)
    input(i, 1) = meanx(i);
    input(i, 2) = meany(i);
end
tform = cp2tform(control, input,'polynomial',4);
tform.forward_fcn = tform.inverse_fcn;

numrpt = size(event_arr,2);
valrptcnt = 0;
clear per vpcind
new_eog_arr=[];
for rptlop = 1:numrpt
    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) >= 1010)) ~=0
        if size(find(event_arr(:,rptlop) == 200)) ~=0
            perbegind = find(event_arr(:,rptlop) == 23,1,'first');
            perendind = find(event_arr(:,rptlop) == 24,1,'first');
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
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

for eye = 1:length(eyedat)
    x = eyedat{eye}(1,:);
    y = eyedat{eye}(2,:);
    [x y] = tformfwd(tform,x,y);
    eyedat{eye} = [x;y];
end

imageX = 756;
imageY = 378;
images = NaN(length(eyedat),1);
initialcross = NaN(length(eyedat),2);
for i = 1:length(eyedat);
    %     figure
    cndline = cnd(i)-1000+1;
    C = textscan(cndfil(cndline,:),'%d');
    itm = C{1}(5)+6;
    cross = C{1}(4)+6;
    C = textscan(itmfil(cross,:),'%f');
    initialcross(i,:) = C{1}(4:5);
    imgind = strfind(itmfil(itm,:),'\');
    img = itmfil(itm,imgind(end)+1:imgind(end)+2);
    images(i) = str2num(img);
    %     img = itmfil(itm,imgind(1)+1:end);
    %     imshow(imread(img));
    C = textscan(itmfil(itm,:),'%f');
    shift = 25*double(C{1}(3));
    x = 25*eyedat{i}(1,:)+imageX/2-shift;
    y = -25*eyedat{i}(2,:)+imageY/2;
    %     hold on
    %     plot(x,y)
    %     plot(25*initialcross(i,1)+imageX/2-shift,-25*initialcross(i,2)+imageY/2,'+r')
    %     hold off
    %     pause
    %     close
    eyedat{eye} = [x;y];
     initialcross(i,:) = [25*initialcross(i,1)+imageX/2-shift,-25*initialcross(i,2)+imageY/2];
end

[~,ia,~] = unique(images,'last');
eyedat(ia) = [];
initialcross(ia,:) = [];

for i = 1:size(eyedat,2);
    x = eyedat{i}(1,:);
    y = eyedat{i}(2,:);
    badx = find(x < -50 | x > imageX+50); %~2 dva leave margin of error
    x(badx) = []; y(badx) = [];
    bady = find(y < -50 | y > imageY+50); %~2 dva margin of error
    x(bady) = []; y(bady) = [];
    eyedat{i} = [x;y];
end
clearvars -except eyedat initialcross

[fixationstats] = ClusterFixation_Final(eyedat);
save('WR130412_1-fixation','fixationstats','initialcross')
%%
cd 'C:\Users\skoenig\Documents\MATLAB\BCRW Salience Model\Shift Task\sft22';
for i = 1:90
    if i < 10
        getSalienceMap(['0' num2str(i) '.bmp']);
    else
        getSalienceMap([num2str(i) '.bmp']);
    end
end
%%
imageX = 756;
imageY = 378;
avgmap  = zeros(imageY,imageX);
for i = 1:90;
    if i < 10
        load(['0' num2str(i) '-saliencemap.mat'],'fullmap')
    else
        load([num2str(i) '-saliencemap.mat'],'fullmap')
    end
    avgmap = avgmap+fullmap;
end
imagesc(avgmap)
axis off