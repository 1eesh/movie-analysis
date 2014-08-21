function [Tracked,AllImg]=DoSeg(directory,channels)
%
% [Tracked,AllImg] = DoTrackSeg(dirname,channel)
%       dirname = name of directory containing all files.
%           e.g.'\Users\ayaron\Movies\2012-01-27\stage1'
%       channel = cell array of strings containing the segmentation color
%           name. e.g. {'YFP'}
%

%global variables
Tracked=[];
AllImg=[];
%

if nargin >2
    error('Wrong number of arguments')
end

if nargin==0
    %open a directory
    directory=uigetdir(pwd,'Select a directory to load files');
end
if ischar(directory) %should start from files
    if ~ischar(directory) || ~isdir(directory)
        error('First argument is not a valid directory name or a Tracked structure')
    end
    if ~exist('channels','var')
        %try to find the possible colors
        olddir=pwd;
        cd(directory);
        f=[dir('*TIF'),dir('*tif')];
        [~,~,~,chlist]=regexp({f.name},'_w[^_.]*');
        chlist=[chlist{:}];
        chlist=cellfun(@(x) x(3:end), unique(chlist),'unif',0);
        chsel=listdlg('ListString', chlist, 'SelectionMode','single', 'ListSize',[160,80]);
        channels=chlist(chsel);
        cd(olddir);
    elseif ~iscellstr(channels)
        error('Second argument is not a valid cell array of strings')
    end
    %do the segmentation
    ParSegment(directory,channels);
else
    error('Wrong arguments')
end


    function ParSegment(directory,channels)
        
        olddir=pwd;
        cd(directory);
        directory=[pwd filesep];
        
        imagefiles={};
        for cchannel=channels
            filelist=[dir(['*' char(cchannel) '*.tif']) dir(['*' char(cchannel) '*.TIF'])];
            filelist=unique({filelist.name});
            [~,~,~,~,fieldpos]=regexp(filelist,'_t(\d+)\.');
            
            fieldpos=[fieldpos{:}];fieldpos=[fieldpos{:}];
            [~,I]=sort(cellfun(@(x) str2num(x),fieldpos));
            imagefiles(end+1,:)=filelist(I);
        end
        
        Tracked=cell(1,length(imagefiles));
        for imagenb=1:length(imagefiles)
            Tracked{imagenb}.filename=imagefiles(:,imagenb);
        end
        
        framestart=1;
        framestop=length(imagefiles);
        
        parallel=matlabpool('size');
        AllImg=cell(1,length(Tracked));
        if parallel==0
            h=waitbar(0,'Segmenting Images...');
        else
            h=[];
            disp('Segmenting Images...')
        end
        parfor imagenb=framestart:framestop
            if parallel==0
                waitbar((imagenb-framestart)/(framestop-framestart),h,sprintf('Segmenting Image %d',imagenb));
            else
                disp(sprintf('Segmenting Image %d',imagenb));
            end
            [ccell,cimg,cfit]=SegmentFile(imagefiles(:,imagenb));
            AllImg{imagenb}=cimg;
            Tracked{imagenb}.cells=num2cell(ccell);
            Tracked{imagenb}.filename=imagefiles(:,imagenb);
            Tracked{imagenb}.dirname=directory;
            Tracked{imagenb}.fit=cfit;
            
        end
        if parallel==0
            delete(h)
        end
        cd(olddir)
    end

end

function [cells,im_bs,cfit]=SegmentFile(filename)
[im_bs,cfit]=loadImage(filename);
cells=SegmentImage(im_bs);
end

function [im_bs,cfit]=loadImage(filename)
if ~exist('cfit','var')
    cfit=[];
end
im={};
for cfilename=filename'
    im{end+1}=imread(char(cfilename));
end
im=sum(cat(3, im{:}),3);
[im_bs,cfit]=RemoveBackground(im,cfit);
end
function [imbs,fitresult]=RemoveBackground(im,cfit)
[imy,imx]=size(im);
if isempty(cfit)
    %fit background
    sx=1:5:imx;
    sy=1:5:imy;
    ss=im(1:5:imy,1:5:imx);
    [xInput, yInput, zOutput] = prepareSurfaceData( sx, sy, ss );
    bkgI=zOutput<prctile(zOutput,10);
    xInput=xInput(bkgI);
    yInput=yInput(bkgI);
    zOutput=zOutput(bkgI);
    % Set up fittype and options.
    ft = fittype( 'p00+p10*x+p01*y+p20*x^2+p11*x*y+p02*y^2', 'indep', {'x', 'y'}, 'depend', 'z' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf];
    opts.StartPoint = [300 0 0 0 0 0];
    opts.Upper = [Inf Inf Inf Inf Inf Inf];
    % Fit model to data.
    [fitresult, gof] = fit( [xInput, yInput], zOutput, ft, opts );
else
    fitresult=cfit;
end
[x,y]=meshgrid(1:imx,1:imy);
%malmmalian
L=fitresult(x,y)-300;
L=max(L,max(max(L))/2)+300;
%185 is about the dark value. maybe 200?
imbs=im./L -185./L+185./min(L(:)) -1;


end

function cells=SegmentImage(im_bs)
im=medfilt2(im_bs,[5,5],'symmetric');

[imy,imx]=size(im);

noise=2*iqr(im(:)-im_bs(:));
segim=zeros(size(im));
ll=localsegment(im);% | find_blob(im,33);
[ll2,nn]=bwlabel(ll,4);

stats = regionprops(ll2, im, 'MeanIntensity','area');
mmean=[stats.MeanIntensity];
area=[stats.Area];
%sstd=arrayfun(@(x) std(x.PixelValues), stats)';
%mmin=arrayfun(@(x) min(x.PixelValues), stats)';

%calculate statistic about the neighborhood of the cells
stats2 = regionprops(ll2, imdilate(im.*(1-ll),strel('diamond',2)), 'MeanIntensity');
stats3 = regionprops(ll2-imerode(ll2,strel('diamond',1)),'area');

Ssignal=(mmean).*sqrt(area)/noise;
Sedge=mmean./[stats2.MeanIntensity];
Scircular=[stats3.Area]./sqrt(area);

signify=Scircular<6;
signify=signify &((Ssignal>6 & Sedge>3) | (Ssignal>100 & Sedge>2));

for cseg=find(signify)
    %segim=segim+imclose((ll2==cseg),strel('disk',3));
    segim=segim+(ll2==cseg);
end

%cells=mask2cells(imclose(segim,strel('disk',5)));
%cells=mask2cells(imopen(imclose(segim,strel('disk',10)),strel('disk',10)));
cells=mask2cells(segim);
end

function imseg=localsegment(im)

%average_size, about the cell size, reducing noise
average_size=30;%15;

%a more uniform image, keep 1e-5 resolution
Th=imtophat(im,strel('disk',average_size));
Th=round(Th*1e5)/1e5;

%where to put the threshold?
%find where there are the most objects (noise segments)
%and take twice that value
localthreshold=2*(fminsearch(@(t) -max(max(bwlabel((Th+1)>t))),1+mode(Th(Th>0))))-2;
%for every cell take it until it is .5 of its maximum
imseg=Th>localthreshold & (Th./imfilter(Th,fspecial('gauss',50,50)))>0.5;

%now use watershed to split neighbouring cells
localmax=imdilate(Th,strel('diamond',10));
localmax=max(min(localmax(localmax>0)),localmax);
imseg=imseg>0;
imseg=waterSegmentImage(Th./localmax,imseg,0.4)>0;
end
function Iseg = waterSegmentImage(I,  Ibin , minima_depth_thresh )

I = double(I);

%Get an image with basins corresponding to white areas in the bin image
Iseg = -I + max(max(I));

%Suppress all minima below a specified value 'minima_depth_thresh'
Iseg  = imhminlocal(Iseg, minima_depth_thresh);

%Each basin is charachterised by a unique label (i.e. number)
Iseg = watershed(Iseg);

%changing all the labels that were 0 in the original image back to 0.
Iseg(Ibin == 0) = 0;
end
function I2 = imhminlocal(varargin)
% Adaptive h-mnima transform, based on imhmin
% instead of substracting a constant h from I, use h=h(I) a monotone
% non-increasing non-negative function of the grey level in I
% see: ftp://doc.nit.ac.ir/cee/y.baleghi/Advanced%20Image%20Processing/Reference%20Books/3D%20Images%20of%20Materials%20Structures%20Processing%20and%20Analysis.pdf
%
% In practice it substract it h-I/3
% so if h is smaller than max(max(I)) it is set to that value
%
%IMHMIN H-minima transform.
%   I2 = IMHMIN(I,H) suppresses all minima in I whose depth is less than
%   H.  I is an intensity image and H is a nonnegative scalar.
%
%   Regional minima are connected components of pixels with the same
%   intensity value, t, whose external boundary pixels all have a value
%   greater than t.
%
%   By default, IMHMIN uses 8-connected neighborhoods for 2-D images and
%   26-connected neighborhoods for 3-D images.  For higher dimensions,
%   IMHMIN uses CONNDEF(NDIMS(I),'maximal').
%
%   I2 = IMHMIN(I,H,CONN) computes the H-minima transform, where CONN
%   specifies the connectivity.  CONN may have the following scalar
%   values:
%
%       4     two-dimensional four-connected neighborhood
%       8     two-dimensional eight-connected neighborhood
%       6     three-dimensional six-connected neighborhood
%       18    three-dimensional 18-connected neighborhood
%       26    three-dimensional 26-connected neighborhood
%
%   Connectivity may be defined in a more general way for any dimension by
%   using for CONN a 3-by-3-by- ... -by-3 matrix of 0s and 1s.  The 1-valued
%   elements define neighborhood locations relative to the center element of
%   CONN.  CONN must be symmetric about its center element.
%
%   Class support
%   -------------
%   I can be of any nonsparse numeric class and any dimension.  I2 has
%   the same size and class as I.
%
%   Example
%   -------
%       a = 10*ones(10,10);
%       a(2:4,2:4) = 7;  % minima 3 lower than surround
%       a(6:8,6:8) = 2;  % minima 8 lower than surround
%       b = imhmin(a,4); % only the deeper minima survive
%
%   See also CONNDEF, IMEXTENDEDMIN, IMHMAX, IMRECONSTRUCT,
%   IMREGIONALMIN.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.6.4.3 $  $Date: 2005/03/31 16:31:28 $

% Testing notes
% -------------
% I       - N-D, real, full
%         - empty ok
%         - Inf ok
%         - NaNs not allowed
%         - logical flag ignored if present
%
% h       - Numeric scalar; nonnegative; real
%         - Inf ok (doesn't make much sense, though)
%         - NaNs not allowed
%
% conn    - valid connectivity specifier
%
% I2      - same class and size as I

[I,H,conn] = ParseInputs(varargin{:});

%H =  0.1*(max_I+(filter2(filt_gauss,max_I-I)));
%figure;imagesc(H);figure;
%%%%

I = imcomplement(I);

%I2 = imreconstruct(imsubtract(I,h), I, conn);
I2 = imreconstruct(imsubtract(I,H), I, conn);

I2 = imcomplement(I2);

%%%
%%% ParseInputs
%%%
    function [I,H,conn] = ParseInputs(varargin)
        
        if verLessThan('matlab', '7.12.1')
            iptchecknargin(2,3,nargin,mfilename);
        else
            narginchk(2,3)
        end
        
        I = varargin{1};
        H = varargin{2};
        
        if verLessThan('matlab', '7.12.1')
            iptcheckinput(I, {'numeric'}, {'real' 'nonsparse'}, mfilename, 'I', 1);
            iptcheckinput(H, {'numeric'}, {'real'}, mfilename, 'H', 2);
        else
            validateattributes(I, {'numeric'}, {'real' 'nonsparse'}, mfilename, 'I', 1);
            validateattributes(H, {'numeric'}, {'real'}, mfilename, 'H', 2);
        end
        
        H = double(H);
        
        if nargin < 3
            conn = conndef(ndims(I),'maximal');
        else
            conn = varargin{3};
            iptcheckconn(conn, mfilename, 'CONN', 3);
        end
        
    end
end

function cells=mask2cells(newsegmask)
[imy,imx]=size(newsegmask);
[newsegmaskl,num]=bwlabel((newsegmask)>0);
newsegclose=imclose(newsegmask,strel('diamond',1));
cells=struct('mask',{},'pos',{},'size',{},'progenitor',{},'descendants',{});
for cseg=1:num
    csegclose=imdilate(newsegmaskl==cseg,strel('diamond',1)).*newsegclose;
    bbox=regionprops(csegclose,'BoundingBox');
    cellslice={floor(bbox.BoundingBox(2))+1:floor(bbox.BoundingBox(2)+bbox.BoundingBox(4)),...
        floor(bbox.BoundingBox(1))+1:floor(bbox.BoundingBox(1)+bbox.BoundingBox(3))};
    cells(end+1).mask=double(csegclose(cellslice{1},cellslice{2}));
     cells(end).pos=floor(bbox.BoundingBox(2:-1:1))+1;
     cells(end).BBox = bbox.BoundingBox;
    cells(end).size=size(cells(end).mask);
    if any(cells(end).pos==1) || any(cells(end).pos-[1,1]+cells(end).size==[imy,imx])
        cells(end)=[];
        continue
    end
end

end
