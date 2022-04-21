function [ F ] = ms_getBiasField( P )
%MS_GETBIASFIELD Summary of this function goes here
%   Detailed explanation goes here
if nargin <1
    P = spm_select;
end

load(P)

nbas           = [size(T) 1];
nbas           = nbas(1:3);
B1             = spm_dctmtx(V(1).dim(1),nbas(1));
B2             = spm_dctmtx(V(1).dim(2),nbas(2));
B3             = spm_dctmtx(V(1).dim(3),nbas(3));


for p=1:V.dim(3)
%     M   = spm_matrix([0 0 p]);
%     img = spm_slice_vol(V, M, V.dim(1:2), 1);
    t   = reshape(T,  nbas(1)*nbas(2), nbas(3));
    t   = reshape(t*B3(p,:)', nbas(1), nbas(2));
%     img = img.*exp(B1*t*B2');
%     img = img./exp(B1*t*B2');
    F(:,:,p)=exp(B1*t*B2');
%     if nargout==0,
%         VO  = spm_write_plane(VO,img,p);
%     else
%         VO.dat(:,:,p) = single(img);
%     end;
end

%% write sth
% Vout=V;
% [d,name,ext]=fileparts(V.fname);
% Vout.fname= [d filesep 'BiasField_' name ext];
% Vout.dt = [spm_type('float32') spm_platform('bigend')];
% spm_write_vol(Vout,F);

end

