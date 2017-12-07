%check image complexity

% sobel filters detect edges!
sobelx = [1     2   1;
    0     0   0;
    -1    -2  -1;];

sobely = [1     2   1;
    0     0   0;
    -1    -2  -1;];

entropyvalues = NaN(1,720);
entropyvaluessal = NaN(1,720);
edgevalues =NaN(1,720);
img_count = 1;

sft_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Shift Task\';
image_sets = {'sft5','sft6','sft8','sft9','sft22','sft24','sft27'};
number_images = [90 120 120 120 90 90 90];

for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([sft_image_dir SETNUM])
    
    for img = 1:number_images(SET)
        if number_images(SET)== 90
            if img < 10
                imgr = imread(['0' num2str(img) '.bmp']);
                load(['0' num2str(img) '-saliencemap.mat'],'fullmap')
            else
                imgr = imread([num2str(img) '.bmp']);
                load([num2str(img) '-saliencemap.mat'],'fullmap')
            end
        else
            if img < 10
                imgr = imread(['00' num2str(img) '.bmp']);
                load(['00' num2str(img) '-saliencemap.mat'],'fullmap')
            elseif img < 100
                imgr = imread(['0' num2str(img) '.bmp']);
                load(['0' num2str(img) '-saliencemap.mat'],'fullmap')
            else
                imgr = imread([num2str(img) '.bmp']);
                load([num2str(img) '-saliencemap.mat'],'fullmap')
            end
        end
        if size(img,3) > 1
            img = rgb2gray(img);
        end
        entropyvalues(img_count) = entropy(imgr);%pixel intesnity entropy
        entropyvaluessal(img_count) = entropy(fullmap);
        xedges = imfilter(imgr,sobelx);
        yedges = imfilter(imgr,sobely);
        edgevalues(img_count) = mean2(xedges+yedges); %edgineess
        
        img_count = img_count+1;
    end
end