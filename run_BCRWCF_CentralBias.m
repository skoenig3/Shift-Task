function run_BCRWCF_CentralBias(CombinedViewingBehaviorData,saliencemapfile,tagname,imageX,imageY,plotoptions)
% Function made by Seth Koenig in 3/27/13. Mod of run_BCRWCF. Unlinke the
% original function this runs the BCRW 10 times for 9 different starting
% positions to determine if central bias depends on starting position. The
% 9 positions are the center followed by N, NE, E, SE, S, SW, W, and SW of
% the central starting point with a distance of 1/4 of the image size.
%
% Function runs the individual BCRWs 90 times based off of the viewing
% behavior statitics from an individual monkey. Each BCRW also requires an
% environment to run in which here it is a saliency map. Dt = 5 ms. 


% INPUTS:
%   1) CombinedViewingBehaviorData: A matrix of viewing behavior statitics for
%   a monkey. Taken from allview{tag} from CombinedViewingBehavior.mat
%   2) Saliencemapfile: The name of teh matlab data file of the saliency map
%   3) Tagname: Monkey/Subuject initals as a string e.g. 'MP'
%   4) imageX: horizontal dimensions of the image
%   5) imageY: vertical dimensions of the image
%   6) plotopions...
%           a) plotoptions.runs = plot 'all' or 'none' of the individual BCRW walks 
%           b) plotoptions.type = for individual runs show 'all' image or
%           'sal'[iency] maps.
%           c) plotoptions.probdens = show 'all' or 'none' of the average
%           results for the 100 BCRW walks. Shows images of the position PDF,
%           fixation PDF, and the saliency map

% OUTPUT:
%   1) Save files named [tagname '-' saliencemapfile '-CB.mat'] with variables...
%       a) fixationtimes: a 3D matrix of run number (row) by trial time
%       (column) sampled at at, and the last 2 dimensions where nonzero
%       when a fixation occured at a time at position x,y. The values of x,y 
%       occupie the elements of the last dimension. 
%       b) fixation : a 2D matrix marked with the number of times the BCRW
%       made a fixation at particular position.
%       c) alltrials: a 2D matrix marked with the number of times the BCRW
%       occupided a particular position. 
%       d) fixations order: a 2D matrix marked with the sum of the fixation order 
%       at fixation locations in a particular position. Honestly not
%       particularly useful. 


warning off all
if nargin < 3
    error(['Not enough inputs: function requires CombinedViewingBehaviorData,'...
        'tagname saliencemapfile, imageX, imageY, and plot options'])
end
if nargin < 5
    imageX = 800;
    imageY = 600;
end
if nargin < 6
    plotoptions.runs = 'none';
    plotoptions.probdens = 'all';
    plotoptions.type = 'sal';
end

data = CombinedViewingBehaviorData;
distance = data.distanceprofile;
probdst = histc(distance,1:max(max(distance)),1);
distCDF = cumsum(probdst,1);
sCDF = sum(probdst,1);
distCDF = bsxfun(@rdivide,distCDF,sCDF);
sacend = data.mediansac;
sacdistance = distCDF(:,1:sacend);
fixdistance = distCDF(:,sacend+1:end);
persistence.sac = nanmean(data.persistence(:,1:sacend));
persistence.fix = nanmean(data.persistence(:,sacend+1:end));
sacduration = data.sacduration(1:end); sacduration(isnan(sacduration)) = [];
probsacduration = hist(sacduration,1:max(sacduration));
sacdurationCDF = cumsum(probsacduration)/sum(probsacduration);
fixduration = data.fixduration(1:end); fixduration(isnan(fixduration)) = [];
probfixduration = hist(fixduration,1:max(fixduration));
fixdurationCDF = cumsum(probfixduration)/sum(probfixduration);
sacangle = data.sacangle; sacangle = sacangle(1:end); sacangle(isnan(sacangle))= [];
sacangle = sacangle*180/pi+180;
nang = (0:360);
probsacangle = hist(sacangle,nang);
sacangleCDF = cumsum(probsacangle)/sum(probsacangle);
load(saliencemapfile,'fullmap');
saliencemap = fullmap;
dash = strfind(saliencemapfile,'-');
dash = dash(1);
clear CombinedViewingBehaviorData data fullmap distCDF

filt = fspecial('gauss',128,12);
% saliencemap(saliencemap < 0.1) = 0; %ignore lowest 10% of salience
saliencemap = imfilter(saliencemap,filt);
[fx, fy] = gradient(saliencemap);

imageX = size(saliencemap,2);
imageY = size(saliencemap,1);
%---Parameters determined using a Parameter Sweep---%
IOR_tau = 1/17; %recovery time of IOR
IOR_area = 48; % area of visual space affected by IOR
border_buffer = 24;
border_sacdist = 48;

nn = 900; %run 100 simulations per fixation crosshair location
trialtime = 5; %5 second run time
[rr, cc] = meshgrid(1:imageX,1:imageY);
dt = 0.005; %5 ms temporal resolution
fixations = zeros(imageY,imageX);
fixationtimes = zeros(nn,trialtime/dt,2);
fixationorder = zeros(600,800);
alltrials = zeros(600,800);
for n = 1:nn;
    if strcmpi(plotoptions.runs,'all')
        figure
        if strcmp(plotoptions.type,'image');
            imagesc(imread([num2str(saliencemapfile(1:dash-1)) '.bmp']))
        else
            imagesc(saliencemap)
        end
        hold on
    end
    fixcount = 1;
%     if n <= 100
%         x = imageX/2;
%         y = imageY/2;
%     elseif n <= 200
%         x = imageX/2;
%         y = imageY/2+imageY/4;
%     elseif n <= 300
%         x = imageX/2+imageX/4;
%         y = imageY/2+imageY/4;
%     elseif n <= 400
%         x = imageX/2+imageX/4;
%         y = imageY/2;
%     elseif n <= 500
%         x = imageX/2+imageX/4;
%         y = imageY/2-imageY/4;
%     elseif n <= 600
%         x = imageX/2;
%         y = imageY/2-imageY/4;
%     elseif n <= 700
%         x = imageX/2-imageX/4;
%         y = imageY/2-imageY/4;
%     elseif n <= 800
%         x = imageX/2-imageX/4;
%         y = imageY/2;
%     elseif n <= 900
%         x = imageX/2-imageX/4;
%         y = imageY/2+imageY/4;
%     end
    if n <= 100
        x = 10;
        y = 189;
    elseif n <= 200
        x = 10;
        y = 368;
    elseif n <= 300
        x = 378;
        y = 10; 
    elseif n <= 400
        x = 378;
        y = 189;
    elseif n <= 500
        x = 378;
        y = 368;
    elseif n <= 600
        x =  746;
        y = 10;
    elseif n <= 700
        x = 746;
        y = 189;
    elseif n <= 800
        x = 746;
        y = 368;
    elseif n <= 900
        x = 10;
        y = 10;
    end
    
    xxyy = [[x;y] zeros(2,9)];
    tmr = 0;
    t = 0;
    angold = [];
    
    
    if IOR_tau > 0;
        previous_fixations = NaN(2,1/IOR_tau);
        previous_fixations(:,end) =[x;y];
    end
    
    
    sacdur = find(rand <= sacdurationCDF);
    sacdur = sacdur(1);
    sacdur(sacdur < 2) = 2;
    timewarp = round(linspace(1,sacend,sacdur));
    %saccade amplitude and over time is controlled
    sacdist = sacdistance(:,timewarp);
    persac = persistence.sac(timewarp);
    %select a random fixation duration
    fixdur = find(rand <= fixdurationCDF);
    fixdur = fixdur(1);
    fixdur(fixdur < 5) = 5;
    timewarp = round(linspace(1,size(fixdistance,2),fixdur));
    fixdist = fixdistance(:,timewarp);
    perfix = persistence.fix(timewarp);
    
    fxx = fx;
    fyy = fy;
    C = sqrt((rr-x).^2+(cc-y).^2)<=IOR_area;
    Cind = find(C);
    fxx(Cind) = 0;
    fyy(Cind) = 0;
    while t < trialtime
        if round(tmr*1/dt)+1 == (sacdur + fixdur+1) %end fixation period so reset model
            xy = ceil(mean(xxyy(:,1:5),2));
            if strcmpi(plotoptions.runs,'all')
                plot(xy(1),xy(2),'*k','markersize',6)
            end
            fixations(xy(2),xy(1)) =  fixations(xy(2),xy(1)) + 1;
            fixationorder(xy(2),xy(1)) = fixationorder(xy(2),xy(1)) + fixcount;
            fixationtimes(n,round(t/dt),:) = xy;
            if IOR_tau > 0;
                xyp = previous_fixations(:,1);
                C = sqrt((rr-xyp(1)).^2+(cc-xyp(2)).^2)<=IOR_area;
                fxx(C) = fx(C);
                fyy(C) = fy(C);
                previous_fixations = [previous_fixations(:,2:end) xy];
            end
            C = sqrt((rr-xy(1)).^2+(cc-xy(2)).^2)<=IOR_area;
            Cind = find(C);
            fxx(Cind) = 0;
            fyy(Cind) = 0;
            tmr = 0;
            xxyy = [[x;y] zeros(2,9)];
            angold = [];
            fixcount = fixcount + 1;
            
            sacdur = find(rand <= sacdurationCDF);
            sacdur = sacdur(1);
            sacdur(sacdur < 2) = 2;
            timewarp = round(linspace(1,sacend,sacdur));
            sacdist = sacdistance(:,timewarp);
            persac = persistence.sac(timewarp);
            fixdur = find(rand <= fixdurationCDF);
            fixdur = fixdur(1);
            fixdur(fixdur < 5) = 5;
            timewarp = round(linspace(1,size(fixdistance,2),fixdur));
            fixdist = fixdistance(:,timewarp);
            perfix = persistence.fix(timewarp);
        end
        if round(tmr*1/dt)+1 <= sacdur;
            dhr = find(sacdist(:,round(tmr*1/dt)+1) >= rand);
            dh = dhr(1);
            b = persac(round(tmr*1/dt)+1)/2;
        else
            dhr = find(fixdist(:,round(tmr*1/dt)-sacdur+1) >= rand);
            dh = dhr(1);
            b = perfix(round(tmr*1/dt)-sacdur+1);
        end
        if tmr == 0; %just starting this simuliation
            angr = find(sacangleCDF >= rand);
            ang = nang(angr(1));
        else
            if x > imageX-border_buffer || x < border_buffer || ...
                    y > imageY-border_buffer || y < border_buffer
                if round(tmr*1/dt)+1 <= sacdur
                    [dh, ang] = border1(dh,x,y,imageX,imageY,dt,tmr,sacdur,...
                        border_sacdist,border_buffer);
                end
            else
                if abs(fyy(y,x)) == 0 && abs(fxx(y,x)) == 0
                    if round(tmr*1/dt)+1 <= sacdur
                        angh = angold;
                    else
                        angr = randi(length(nang));
                        angh = nang(angr);
                    end
                elseif  abs(fyy(y,x)) > 0 && abs(fxx(y,x)) > 0
                    angh = atand(fyy(y,x)/fxx(y,x));
                    if fxx(y,x) < 0
                        angh = angh + 180;
                    end
                elseif abs(fyy(y,x)) > 0
                    if fyy(y,x) < 0;
                        angh = 90;
                    else
                        angh = 270;
                    end
                elseif abs(fxx(y,x)) > 0
                    if fxx(y,x) < 0;
                        angh = 180;
                    else
                        angh = 0;
                    end
                end
                if angold > 360
                    angold = angold - 360;
                elseif angold < -360
                    angold = angold + 360;
                end
                if abs(angh-angold) > 180
                    if angh < angold
                        angh = angh+360;
                    elseif angold < angh
                        angold = angold+360;
                    end
                end
                ang = angold*(1-b) + b*angh;
            end
        end
        xn = round(x + dh*cos(ang*pi/180));
        yn = round(y + dh*sin(ang*pi/180));
        angold = ang;
        if (xn > imageX || xn < 1 || yn < 1 || yn > imageY)
            [xn, yn, angold] = border2(x,xn,y,yn,dh,imageX,imageY,dt,tmr,sacdur,...
                border_sacdist,border_buffer);
        end
        if strcmpi(plotoptions.runs,'all')
            plot([x xn],[y yn],'m')
            plot(xn,yn,'.m','markersize',3)
            %             pause(0.01)
        end
        alltrials(yn,xn) = alltrials(yn,xn) + 1;
        x = xn;
        y = yn;
        tmr = tmr + dt;
        t = t+dt;
        xxyy =  [[x;y]  xxyy(:,1:9)];
    end
    if strcmpi(plotoptions.runs,'all')
        close
    end
end

filt = fspecial('gauss',32,6);
if strcmpi(plotoptions.probdens,'all')
    figure
    f = 5*alltrials;
    f = imfilter(f,filt,'replicate');
    f(f<1) = 1;
    imagesc(log(f))
    imagesc(alltrials)
    title('PDF: All Positions')
    figure
    f = 100*fixations;
    f = imfilter(f,filt,'replicate');
    f(f<1) = 1;
    imagesc(log(f))
    title('PDF: Fixations')
    figure,imagesc(saliencemap)
    title('Saliencemap')
end
    function   [dh, ang] = border1(dh,x,y,imageX,imageY,dt,tmr,sacdur,...
            border_sacdist,border_buffer)
        if x >= imageX-border_buffer && y >= imageY-border_buffer
            ang = 225;
            dh(dh < border_sacdist*sqrt(2)) = border_sacdist*sqrt(2);
        elseif x >= imageX-border_buffer && y <= border_buffer
            dh(dh < border_sacdist*sqrt(2)) = border_sacdist*sqrt(2);
            ang = 135;
        elseif x <= border_buffer && y >= imageY-border_buffer
            ang = 315;
            dh(dh < border_sacdist*sqrt(2)) = border_sacdist*sqrt(2);
        elseif x <= border_buffer && y <= border_buffer
            ang = 45;
            dh(dh < border_sacdist*sqrt(2)) = border_sacdist*sqrt(2);
        elseif x <= border_buffer
            ang = 0;
            dh(dh < border_sacdist) = border_sacdist;
        elseif y <= border_buffer
            ang = 90;
            dh(dh <  border_sacdist) = border_sacdist;
        elseif x >= imageX-border_buffer
            ang = 180;
            dh(dh < border_sacdist) = border_sacdist;
        elseif y >= imageY-border_buffer
            ang = 270;
            dh(dh < border_sacdist) = border_sacdist;
        end
        %         end
    end
    function [xn, yn, angold] = border2(x,xn,y,yn,dh,imageX,imageY,dt,tmr,sacdur,...
            border_sacdist,border_buffer)
        xx = [x xn]; yy = [y yn];
        p = polyfit([x xn],[y yn],1);
        if xn > imageX
            y = round(p(1)*imageX + p(2));
            x = imageX;
        elseif xn < 1
            y = round(1*p(1)+p(2));
            x = 1;
        elseif yn < 1
            if any(isinf(p));
                x = x; %#ok
                y = 1;
            else
                x = round((0-p(2))/p(1));
                y = 1;
            end
        elseif yn > imageY
            if  any(isinf(p));
                x = x; %#ok
                y = imageY;
            else
                x = round((imageY-p(2))/p(1));
                y = imageY;
            end
        else
            x = xn; y = yn;
        end
        x(x > imageX) = imageX;
        y(y > imageY) = imageY;
        [dhn, angn] = border1(dh,x,y,imageX,imageY,dt,tmr,sacdur,...
            border_sacdist,border_buffer);
        xn = round(x + dhn*cos(angn*pi/180));
        yn = round(y + dhn*sin(angn*pi/180));
        xn(xn > imageX) = imageX; xn(xn < 1) = 1;
        yn(yn > imageY) = imageY; yn(yn < 1) = 1;
        angold = angn;
    end

save([tagname '-' saliencemapfile(1:dash-1) '-CB.mat'],'alltrials',...
    'fixations','fixationorder','fixationtimes')
end