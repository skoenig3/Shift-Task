function avgreturnmap = Bias_of_Image_Return
% If the monkey leaves the image where does in the image do they return to?
% Essentially is the central bias mediated by a bias to return to the
% center of the image if the monkey looks outside of the image.

sft_image_dir = 'C:\Users\skoenig\Documents\MATLAB\BCRW Salience Model\Shift Task\';
image_sets = {'sft22','sft23','sft24','sft25','sft26','sft27'};


imageX = 756;
imageY = 378;
avgreturnmap = zeros(imageY,imageX);
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
        returnmap = getSFTeyedat(cortexfile,ITMFile,CNDFile,imageX,imageY);
        avgreturnmap = avgreturnmap + returnmap;
    end
end
end

function [returnmap] = getSFTeyedat(cortexfile,ITMFile,CNDFile,imageX,imageY);
% extracts raw Shift task (SFT) eye data from a cortex files, calibrates
% the eye data, and then uses Cluster Fix to extract fixations and saccades

returnmap = zeros(imageY,imageX);
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

[time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(cortexfile);

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
    img = itmfil(itm,imgind(end)+1:imgind(end)+2);
    images(i) = str2double(img);
    C = textscan(itmfil(itm,:),'%f');
    shift = 25.2*double(C{1}(3));
    x = 25.2*eyedat{i}(1,:)+imageX/2-shift;
    y = -25.2*eyedat{i}(2,:)+imageY/2;
    eyedat{i} = [x;y];
    initialcross(i,:) = [25.2*initialcross(i,1)+imageX/2-shift,-25.2*initialcross(i,2)+imageY/2];
end

[~,ia,~] = unique(images,'last');
eyedat(ia) = [];
initialcross(ia,:) = [];

for i = 1:size(eyedat,2);
    x = round(eyedat{i}(1,:));
    y = round(eyedat{i}(2,:));
    outsideind = find(x < 1 | x > imageX | y < 0 | y >= imageY);
    gaps = find(diff(outsideind) > 1);
    if isempty(gaps);
        if length(outsideind) >= 10;
            reind = outsideind(end)+1:outsideind(end)+40;
            reind(reind > length(x)) = [];
            if ~isempty(reind);
                imgind = sub2ind([imageY,imageX],...
                    imageY-y(reind),x(reind));
                returnmap(imgind) = returnmap(imgind)+1;
            end
        end
    else
        for ii = 1:length(gaps);
            if ii == 1;
                if gaps(1) >= 10;
                    if length(gaps) == 1;
                        if outsideind(gaps(1)+1)-outsideind(gaps(1)) > 40;
                            reind = outsideind(1:gaps(1));
                            imgind = sub2ind([imageY,imageX],...
                                imageY-y(reind(end)+1:reind(end)+40),...
                                x(reind(end)+1:reind(end)+40));
                            returnmap(imgind) = returnmap(imgind)+1;
                        end
                    else
                        if outsideind(gaps(ii)+1)-outsideind(gaps(ii)) > 40
                            reind = outsideind(1:gaps(1));
                            imgind = sub2ind([imageY,imageX],...
                                imageY-y(reind(end)+1:reind(end)+40),...
                                x(reind(end)+1:reind(end)+40));
                            returnmap(imgind) = returnmap(imgind)+1;
                        end
                    end
                end
            elseif ii == length(gaps);
                if length(outsideind)-gaps(end) >= 10;
                    reind = outsideind(gaps(end)+1:end);
                    reind = reind(end)+1:reind+40;
                    reind(reind > length(x))= [];
                    if ~isempty(reind)
                    imgind = sub2ind([imageY,imageX],...
                        imageY-y(reind),x(reind));
                    end
                    returnmap(imgind) = returnmap(imgind)+1;
                end
            else
                if gaps(ii)-gaps(ii-1) >= 10;
                    if outsideind(gaps(ii)+1)-outsideind(gaps(ii)) > 40
                        reind = outsideind(gaps(ii-1)+1:gaps(ii));
                        imgind = sub2ind([imageY,imageX],...
                            imageY-y(reind(end)+1:reind(end)+40),...
                            x(reind(end)+1:reind(end)+40));
                        returnmap(imgind) = returnmap(imgind)+1;
                    end
                end
            end
        end
    end
end
end