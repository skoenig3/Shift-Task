%---[5] Extract Viewing Behavior---%
sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\';
image_sets = {'sft5','sft6','sft8','sft9','sft22','sft24','sft27'};




%define diagonals for 1-back,2-back,3-back
id0 =  eye(40);
id0 =  [id0(2:end,:); zeros(1,40)];
id1 = eye(40);
id1 = [id1(3:end,:); zeros(2,40)];
id2 = eye(40);
id2 = [id2(4:end,:); zeros(3,40)];
id3 = eye(40);
id3 = [id3(5:end,:); zeros(4,40)];

all_nov_recurence_map = zeros(40,40);

imgdur = [];

distance_threshold = 48;
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
        load(matfiles.mat{eyefile});% loads fixaqtionstats
        
        for cndlop = 1:length(fixationstats); %only uses novel viewing since images changes could alter natural behavior
            fixations = fixationstats{cndlop}.fixations;
            if ~isempty(fixations)
                fixationtimes = fixationstats{cndlop}.fixationtimes;
                saccadetimes =  fixationstats{cndlop}.saccadetimes;
                xy = fixationstats{cndlop}.XY;
                if fixations(1,1) > initialcross(cndlop,1)-100 && fixations(1,1) < initialcross(cndlop,1)+100 &&...
                        fixations(2,1) < initialcross(cndlop,2)+100 && fixations(2,1) > initialcross(cndlop,2)-100
                    fixations(:,1) = [];
                    fixationtimes(:,1) = [];
                end
                
                
                nov_recurence_map = zeros(40,40);
                N=size(fixations,2);
                N(N > 40) = 40;
                [x,y]=meshgrid(1:N);
                i=find(ones(N)); %forms pairs except for self-pairing
                i=[x(i), y(i)];
                dist =sqrt((fixations(1,i(:,1))-fixations(1,i(:,2))).^2 +...
                    (fixations(2,i(:,1))-fixations(2,i(:,2))).^2);
                dind = find(dist <= distance_threshold);
                for d = 1:length(dind);
                    nov_recurence_map(i(dind(d),1),i(dind(d),2)) = nov_recurence_map(i(dind(d),1),i(dind(d),2))+1;
                end
                
                all_nov_recurence_map = all_nov_recurence_map+nov_recurence_map;
                
            end
        end
    end
end
%%
rm_nov = all_nov_recurence_map;
for r = 1:length(rm_nov);
     fix_count = rm_nov(r,r); 
     for i = 1:r-1
        rm_nov(r,i) = rm_nov(r,i)/fix_count;
        rm_nov(i,r) = rm_nov(i,r)/fix_count;
     end
     rm_nov(r,r) = 1;
end

rm_nov(eye(40) == 1) = NaN;
rm_nov(id0 == 1) = NaN;
rm_nov(id0' == 1) = NaN;
rm_nov(id1== 1) = NaN;
rm_nov(id1' == 1) = NaN;
%%
rm_nov = rm_nov(1:20,1:20);