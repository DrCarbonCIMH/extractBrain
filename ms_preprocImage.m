function [ S, Vproc ] = ms_preprocImage( Mtx, V, voxdim, style )
%MS_PREPROCIMAGE Summary of this function goes here
%   The image processing is not always helpful!
%   In case you have already a 'clean' image (bc, volumecoil, etc.) it
%   sometimes better just to do the bare minimum
if nargin<4
    style='full';
end

fprintf('ms_preprocImage: performing image preprocessing.\n');
if(strcmp(style,'full'))
    S=Mtx; S=S./max(S(:));
    %% what about the denoising?
    % [ noise ] = ms_DenoiseMR(S); S=S-0.2*noise;
    %% image sharpening; works well.. a ~ 0.5
    if prod(voxdim)<10; a=0.5; S = S + a*S - a*imgaussfilt3(S,1); end % not sure if also good for EPIs
    %% can we identify noise? don't know if this is helpful!
    % [N,X]=hist(S(:),1024);
    % noise=X(find(cumsum(N)<0.6*numel(S),1,'last')); S(S<=noise)=0;
    %% kmeans clustering helpful? somewhat.. but not really!
    % [IDX, C] = kmeans(S(:), 2);
    % IDX=reshape(IDX,size(S)); [~,idx]=sort(C); S(IDX==idx(1))=0;
    %% laplacian filter
    
    %% make a pseudo bias correction; keep the gaussian filter low 1-3;
    if prod(voxdim)<10; S=S./(imgaussfilt3(S,2)+0.001); end % this gaussian approach works well for biased data, but not for EPIs!
    S = reshape(imadjust(S(:)),size(S)); % lowers the needed iterations
    %% gradient image helpful? yeah.. it is.. sth around -0.5*Gmag seems okay; esp. for EPIs
    % [Gmag] = (imgradient3(S)); Gmag=Gmag./max(Gmag(:)); % maybe the imgradient should be performed on Mtx?
    % [Gmag] = (imgradient3(Mtx./max(Mtx(:)))); Gmag=Gmag./max(Gmag(:));
    adj=reshape(imadjust(Mtx(:)./max(Mtx(:))),size(Mtx)); [Gmag] = (imgradient3(adj)); Gmag=Gmag./max(Gmag(:));
    if prod(voxdim)<10; S=S-0.4*Gmag; else; S=S-0.6*Gmag; end
    
    %% bring it back to max = 1
    S=S./max(max(max(S)));
end
if(strcmp(style,'minimal'))
    % bare minimum
    S=Mtx./max(max(max(Mtx)));
    S = reshape(imadjust(S(:)),size(S));
end
%% write the processed image
Vproc = V; Vproc.dt(1)=16;
[d,name,ext]=fileparts(V.fname);
Vproc.fname = [d filesep name '_procInp' ext]; 
spm_write_vol(Vproc, S);

%%
% S = reshape(imadjust(S(:), [0.06 0.9],[]),size(S));

% S=imcomplement(imgradient3(imgradient3(S))); S=S./max(max(max(S))); 

% lap = ones(3,3,3).*-1; lap(2,2,2)=8; %// Change - Centre is now positive
% resp = imfilter(S, lap, 'conv'); %// Change
% %// Change - Normalize the response image
% minR = min(resp(:));
% maxR = max(resp(:));
% resp = (resp - minR) / (maxR - minR);
% S = S + resp;
%// Change - Normalize the sharpened result
% minA = min(S(:));
% maxA = max(S(:));
% S = (S - minA) / (maxA - minA);
% S=S./max(max(max(S)));


% S = reshape(histeq(S(:)),size(S)); % doesn't work
% S(S<0.1)=0;
% get the image details
% FFT=ifftshift(ifftn(S));
% [i1,i2,i3]=ind2sub(size(S),find(abs(FFT)==max(max(max(abs(FFT))))));
% FFT(i1-2:i1+2,i2-2:i2+2,i3-2:i3+2)=zeros(5,5,5);
% E=abs(fftn(fftshift(FFT)));
% get the edges
% for ix=1:size(S,3)
%     E(:,:,ix)=edge(S(:,:,ix),'canny');
% end
% for ix=1:size(S,2)
%     E(:,ix,:)=squeeze(E(:,ix,:))+edge(squeeze(S(:,ix,:)),'canny');
% end

% S=medfilt3(S,[2 2 2]);
%% additional stuff
% get the edges
% S=Mtx.*mask; S=S./max(S(:));
% cThr=[0.2 0.9];
% S=mask.*T; cThr=[];
% for ix=1:size(S,3)
%     E(:,:,ix)=edge(S(:,:,ix),'canny',cThr);
% end
% for ix=1:size(S,2)
%     E(:,ix,:)=squeeze(E(:,ix,:))+edge(squeeze(S(:,ix,:)),'canny',cThr);
% end
% for ix=1:size(S,1)
%     E(ix,:,:)=squeeze(E(ix,:,:))+edge(squeeze(S(ix,:,:)),'canny',cThr);
% end
% % E=~E; E(mask==0)=0; se = strel('sphere',1); E = imerode(E,se);
% Vedge = V;
% Vedge.fname = [d filesep name '_brainEdge' ext]; 
% spm_write_vol(Vedge, double(E));
% 
% % get boundary of image
% se = strel('sphere',1);
% Bou =  mask - imerode(mask,se);
% Vedge = V;
% Vedge.fname = [d filesep name '_brainBoundary' ext];
% spm_write_vol(Vedge, double(Bou));
% % get gradient of image
% [Gmag, Gaz, Gelev] = imgradient3(imgradient3(S)); Gmag=Gmag./max(Gmag(:)); %Gmag=reshape(imadjust(Gmag(:),[0.2 0.8],[]),size(Gmag));
% Vedge = V; Vedge.fname = [d filesep name '_brainGradient' ext];
% Gmag=imcomplement(Gmag); Gmag(logical(Bou))=0;
% spm_write_vol(Vedge, double(Gmag));
% 
% test=E;
% test=Gmag; test(test==0)=2; test(test<2)=0;
% CC = bwconncomp(test,6);
% numPixels = cellfun(@numel,CC.PixelIdxList);
% [biggest,idx] = max(numPixels);
% idx
% numPixels
% BW = zeros(size(E));
% BW(CC.PixelIdxList{idx}) = 1;
% Vedge = V; Vedge.fname = [d filesep name '_brainGradient2' ext];
% spm_write_vol(Vedge, double(BW));
end

