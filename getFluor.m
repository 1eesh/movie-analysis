function [YFPmean, YFPmedian, totalCellYFPFluor, RFPmean, RFPmedian, totalCellRFPFluor] = getFluor(tracked, position, numFrames)

% Tracked is the output of Yaron's code.

totalCellYFPFluor  = cell(1,numFrames);
totalCellRFPFluor  = cell(1,numFrames);
YFPmedian = zeros(1, numFrames);
RFPmedian = zeros(1, numFrames);
YFPmean = zeros(1, numFrames);
RFPmean = zeros(1, numFrames);

for i=1:numFrames
    
    numCells = length(tracked{i}.cells);
    
    fileNameBegin = regexp(tracked{i}.filename, 'w','split');
    fileNameEnd = regexp(fileNameBegin{1}(2), 's', 'split');
    YFPfileName = [fileNameBegin{1}{1}, 'w3YFP_s', fileNameEnd{1}{2}];
    RFPfileName = [fileNameBegin{1}{1}, 'w2RFP_s', fileNameEnd{1}{2}];
    
    imY = imread([tracked{i}.dirname, YFPfileName]);
    %     imYsort = sort(imY(:));
    %     imYbg = mean(imYsort(1:20000));
    [imY,fitresultY]=MamRemoveBackground2(imY, []);
    
    
    imR = imread([tracked{i}.dirname, RFPfileName]);
    %     imRsort = sort(imR(:));
    %     imRbg = mean(imRsort(1:50000));
    [imR,fitresultR]=MamRemoveBackground2(imR, []);
    
    for j=1:numCells
        rect = [tracked{i}.cells{j}.BBox];
        
        croppedimY = imcrop(imY, rect);
        imCellY = croppedimY(tracked{i}.cells{j}.mask==1);
        totalCellYFPFluor{i}(j) = sum(imCellY(:)); % - length(imCellY)*imYbg ;
        
        croppedimR = imcrop(imR, rect);
        imCellR = croppedimR(tracked{i}.cells{j}.mask==1);
        totalCellRFPFluor{i}(j) = sum(imCellR(:)); % - length(imCellR)*imRbg ;
    end
    
    YFPmean(i) = mean(totalCellYFPFluor{i}(:));
    RFPmean(i) = mean(totalCellRFPFluor{i}(:));
    YFPmedian(i) = median(totalCellYFPFluor{i}(:));
    RFPmedian(i) = median(totalCellRFPFluor{i}(:));
    
end

eval([position, '_YFPmean = YFPmean']);
eval([position, '_RFPmean = RFPmean']);
eval([position, '_YFPmedian = YFPmedian']);
eval([position, '_RFPmedian = RFPmedian']);


save(position, [position,'_YFPmean'], [position, '_RFPmean'], [position, '_YFPmedian'], [position, '_RFPmedian']);


end

