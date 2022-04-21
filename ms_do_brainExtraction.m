function [ Pout, finalBrainSize ] = ms_do_brainExtraction( P, BrSize, flag )
%ms_do_brainExtraction Summary of this function goes here
%   PCNN is an iterative procedure; it activates neurons/pixels based on
%   certain image properties (starting with the highest value and then influencing neighbours). 
%   the algorithm performs the PCNN stepwise and collects actived pixels; these pixels serve as mask;
%   idea: find the biggest cluster of activated pixels and define that as brain mask
%   you will see that the brain mask follows a decreasing slope till it
%   "breaks" and the brain mask combines with 'outer clusters' 
%   check the iterations with 'ms_gui_checkBrainMasks'
%   ** flag **
%   you can set some parameters using the 'flag' variable (see "set some
%   'standard' options" below). Of most interest is probably..
%   - the p value which influences the smoothing of the mask; p=5 seems to be a good value -> but sometimes 3 is better (more details);
%   - r influences the PCNN itself; best to keep it 3 (or 2)
%   - preprocStyle is usually 'full' -> 'minimal' can be better if your image is already a nice one (e.g. bias-corrected, volume coil,..)
% 
%   the best parameter set is maybe more like an art.. but try to make the 'plateau' as wide as possible 

% add path if needed
if isempty(which('findit')); addpath(genpath([fileparts(which('ms_do_brainExtraction')) filesep 'PCNN3D_matlab'])); end

%check if the mex version works
mexworks=1;
try
    calcPCNN_mex(1,1,1,1,1,1);
catch
    fprintf('mex version of calcPCNN does not work. Using m-file instead\n');
    mexworks=0;
    addpath(genpath([fileparts(which('ms_do_brainExtraction')) filesep 'origCalcPCNN']));
end
    

if nargin <1
    P=spm_select;
end

if nargin <2
    if ~isempty(strfind(P,'ZI_M')) % that should work; maybe better to calculate the whole volume size
        BrSize=[350,650]*1000; % brain size range for mouse (mm3). If you use 10x data, this should be 10x too (for each dimension!)
    else
        BrSize=[1200,3000]*1000; % brain size range for RAT (mm3)
    end
end
%% load and get voxel dimensions
V=spm_vol([P ',1']); %V=V(1);
Mtx=spm_read_vols(V); 
voxdim=spm_imatrix(V.mat); voxdim=abs(voxdim(7:9));
[d,name,ext]=fileparts(P);
totVol=numel(Mtx)*prod(voxdim)/1000;

%% set some 'standard' options
% r is for PCNN; p for morphological smoothing; %in the paper they mention r=3 and p=4; but it seems r=p (more or less confirmed)
if prod(voxdim)<10 % low voxel dimensions indicate anatomical image
    r=3; p=5; % for anatomical r/p = 3/5 seems okay; if not.. lower p (p=3); seems to be better for rats sometimes
else
    r=2; p=3; % for EPIs it's more 2/2-3
end   
verb = 0; % verbose mode
ToSkip = 0; % iterations to skip the actual brain mask calculation; but PCNN runs! NOT recommended!
showResult = 1; % final check_reg and iteration figure
preprocStyle = 'full'; % the image preprocessing can help. But not always.. 'full'/'minimal'
if nargin == 3 % flag is used to set options
    fprintf('User settings:\n')
    fn = fieldnames(flag);
    for i=1:length(fn)
        eval([fn{i} '= flag.' fn{i} ';'])
        if(isnumeric(flag.(fn{i}))); fprintf([fn{i} ' =  %d\n'], flag.(fn{i})); else; fprintf([fn{i} ' =  %s\n'], flag.(fn{i})); end
    end
end

%% check if BrSize is set correctly
if BrSize(2)>=length(find(Mtx))*prod(voxdim)
    BrSize(2)=floor(length(find(Mtx))*prod(voxdim));
    fprintf('upper Limit of BrainSize is bigger than the provided image!\n Setting BrainSize to %.2f\n', BrSize(2)/1000)
end

%% preprocess the image
% ms_preprocImage tries to preprocess the image in a way so that the PCNN
% algorithm (or brain extraction itself) works better
[S, Vproc] = ms_preprocImage(Mtx, V, voxdim, preprocStyle);
%% "load/define" stuff we need for PCNN and biological smoothing
% M will be used for PCNN calculation
M = ms_3Dgaussian(r, voxdim./min(voxdim), 0.5); % more or less the same as scaledgauss
% NH (or se) is used for biological smoothing
NH = ms_3Dsphere(p, voxdim./min(voxdim));
se = strel('arbitrary', NH); % adapt it to the actual voxel dimensions

% setup the PCNN variables
Y=zeros(size(S)); F = Y; L = F;  T = Y + 1; %U = F;
% setup the brain mask
A = logical(Y); BW=false([60, size(S,1), size(S,2), size(S,3)]); % probably we won't need more than 60 iterations
% initialize some parameters
t = 1; % to keep track of brainmask evolution
brainSize = 0; % currently calculated brain size -> break condition 
tic
%% perform the "actual" brain extraction
while brainSize<BrSize(2)/1000 % BrSize(2) is the set upper limit we expect for the brain size
    if(mexworks)
        [F,L,T,Y] = calcPCNN_mex(F,L,T,Y,M,S); % calc the next PCNN step
    else
        [F,L,T,Y] = calcPCNN(F,L,T,Y,M,S); % calc the next PCNN step
    end
    A = A | logical(Y); % keep track of "activated" voxels
    if t>ToSkip % to save time we could skip the first iterations
        % perform biological smoothing and find the cluster most reasonable as being the brain mask 
        [BW(t+1,:,:,:)] = calcBrainMask(A,se); % this takes the most time ~1-2s for high res images
        % calculate the current brain size
        vox(t+1)=length(find(BW(t+1,:,:,:)))*prod(voxdim)/1000; brainSize=vox(t+1);
        fprintf('#%.0f: BrainSize: %.2f mm^3 (%.2f of total Volume)\n',t, vox(t+1), vox(t+1)/totVol)
        % in case we want to see details -> plot it
        if verb; plotSomething(squeeze(BW(t+1,:,:,:)), A, Y, T,F,L,vox,t,Mtx); end
        % it can happen that the algorithm is "trapped".. but very rare
        if t>30 && vox(t)>0
            if vox(t)<=vox(t-10)*1.01
                vox(end)=BrSize(2)/1000*1.1; warndlg(sprintf('I was trapped and stopped the iteration!\n%s',name)); break; 
            end
        end
    end
    t=t+1;
end
toc
%% find the "optimal" iteration step; that is more a guess.. in case you're not happy use 'ms_gui_checkBrainMasks'
optIteration = findOptIter(vox, BrSize, BW, S, showResult, name);
mask=squeeze(BW(optIteration,:,:,:)); %mask(~mask)=NaN;
finalBrainSize = vox(optIteration);

%% write a brain and mask
Vbrain = V;
Vbrain.fname = [d filesep name '_brain' ext]; Pout=Vbrain.fname;
spm_write_vol(Vbrain, Mtx.*mask);
Vmask = V;
Vmask.fname = [d filesep name '_brainmask' ext]; 
spm_write_vol(Vmask, mask);
if showResult
%     spm_check_registration(V,Vproc,Vbrain, Vmask)
    spm_check_registration(char(V.fname,Vproc.fname,Vbrain.fname, Vmask.fname)); % this works with spm8 and spm12
end

%% save the masks
BW(t+2:end,:,:,:)=[];
BC=V.fname; I_border = BW; G_I = vox; optG = optIteration; % renaming to be in line with the gui 
save([d filesep 'BrainMasks_' name '.mat'], 'I_border', 'G_I', 'optG','BC');

end

function plotSomething(BW, A, Y, T,F,L,vox,t,IM)
figure(22);
sl=floor(size(A,3)/2)+4;
subplot(2,4,1); imagesc(squeeze(BW(:,:,sl))); title(['Mask, n=' num2str(t)]);
subplot(2,4,2); imagesc(squeeze(A(:,:,sl))); title(['A']);
subplot(2,4,3); imagesc(squeeze(Y(:,:,sl))); title(['Y']);
subplot(2,4,4); imagesc(squeeze(T(:,:,sl))); title(['T']);
subplot(2,4,5); imagesc(double(squeeze(IM(:,:,sl)))); title(['Image']);
subplot(2,4,6); plot(vox); title(['voxel']);
subplot(2,4,7); imagesc(squeeze(F(:,:,sl))); title(['F']);
subplot(2,4,8); imagesc(squeeze(L(:,:,sl))); title(['L']);
pause(0.1)
end

function [BW] = calcBrainMask(A, se)
% this function "detects"/calculates the brain mask using the activated
% voxels A and se as biological smoothing parameter
% set the mask to false
BW = false(size(A));
% erode the activated voxel map
A = imerode(A,se);
% find all cluster in A
CC = bwconncomp(A,26); % 6 or 26? it doesn't make a big difference but it seems in the paper they use 26
numPixels = cellfun(@numel,CC.PixelIdxList);
% take the biggest cluster and see this as "brain"
[biggest,idx] = max(numPixels);
if ~isempty(idx)
    BW(CC.PixelIdxList{idx}) = true; 
    BW=imdilate(BW,se); % dilate the brain cluster
    % remove holes within the cluster
    for ix=1:size(BW,3); BW(:,:,ix)=imfill(BW(:,:,ix),'holes'); end % that is quicker, not really "correct"; but it works okay
end
% MS_Oct20: it can happen, that the "biggest" is not the brain but stuff around..
if(biggest>10000)
    midPointImage = zeros(size(A)); midPointImage(round(size(A,1)/2)-5:round(size(A,1)/2)+5, round(size(A,2)/2)-5:round(size(A,2)/2)+5, round(size(A,3)/2)-5:round(size(A,3)/2)+5)=1;
    if(~sum(sum(sum(midPointImage.*BW))))
        fprintf('"biggest" does NOT cover the middle of the image! I take the second biggest instead\n')
        [~,I]=sort(numPixels,'descend');
        idx=I(2);
        BW(CC.PixelIdxList{idx}) = true; 
        BW=imdilate(BW,se); for ix=1:size(BW,3); BW(:,:,ix)=imfill(BW(:,:,ix),'holes'); end
    end
end

end

function optIteration = findOptIter(vox, BrSize, BW, S, showResult, name)
%% after the PCNN iteration is done we have to find the "best" iteration step (or the best mask)
% a well performed iteration shows an increasing brain size whose slope is
% getting flatter and flatter till a "breaking point" in which the brain
% size "explodes" -> the best brain mask definition is usually 1-2 steps
% before this break point.

% can we use the gradient field to get rid of these little blobs? they
% often disappear some iterations before: works sometimes..
mask=squeeze(BW(end,:,:,:));
[Gmag] = (imgradient3(S.*mask)); Gmag=Gmag./max(Gmag(:)); se = strel('sphere',1);
for ix=0:10
    mask=squeeze(BW(end-ix,:,:,:)); Bou =  mask - imerode(mask,se); % get the boundary of the mask
    tmp=Gmag(logical(Bou)); %tmp=S(logical(Bou));
    OL(ix+1)=mean(tmp(tmp>0));%/voxnorm(optIteration-ix); % the highest overlap value should do the trick
end
[~,idx]=sort(OL,'descend'); optIteration=length(vox)-idx(1)-1;
fprintf('Gradient operation says %.0f or maybe %.0f\n', optIteration,length(vox)-idx(2)-1)
% paper approach
[paper]=findit(vox*1000,BrSize); % that is the paper approach
fprintf('The orig. paper says %.0f [%.0f %.0f]\n', paper(2), paper(1),paper(3)); 

% another approach: find the plateau
voxsmall=vox(vox>BrSize(1)/1000);
stdvox=stdfilt(diff(voxsmall)./voxsmall(2:end)); pl=find(stdvox<min(stdvox)*6); pl=pl+(length(vox)-length(voxsmall)); % I don't like this "6".. but don't know how to make it better
if ~isempty(pl)
    optIteration=pl(1)+round((pl(end)-pl(1))/2);
else
    fprintf('Own approach failed! Using the paper approach.\n');
    optIteration=paper(2); pl =[paper(1) paper(3)];
end
fprintf('I guess the best Iteration is number %.0f; Brain Volume: %.2f\n', optIteration, vox(optIteration))
if showResult
    figure(10); plot(vox); hold on; plot(optIteration, vox(optIteration), 'ro'); line([pl(1) pl(1)], [0 BrSize(2)/1000], 'Color', 'g'); line([pl(end) pl(end)], [0 BrSize(2)/1000], 'Color', 'g'); title(name, 'Interpreter', 'none'); hold off;
end
% save('vox.mat', 'vox')
end
