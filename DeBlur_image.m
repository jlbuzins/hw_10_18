function [X,XKP, kc,kr,kp] = DeBlur_image(Ac,Ar,B, varargin)
% [B,Ac,Ar,X,Btrue,E] = Blur_image(varargin)
%        B = Btrue + E, Btrue=Ac X Ar^T
% Input variables
% Ac - column blur
% Ar - row blur
% B  - blurred image
% varargin are as listed in this order
% varargin{1} tolc -tolerance for the column SVs- else set to .01 of max
% varargin{2} tolr  tolerance for the row SVs - else set to .01 of max
% varargin{3} tolkp tolerance on the KP SVs - else set to .01 of max
% Output variables:
% X     - the restored image using tolc and tolr
% XKP   - resorted image using tolKP
% kc    - truncation parameter of Sc
% kr    - truncation parameter of Sc
% kp    - truncation parameter of KPSVS
%%
% Find the SVD
[Uc,Sc,Vc]=svd(Ac);
[Ur,Sr,Vr]=svd(Ar);
Sc=diag(Sc); %extract into diagonal matrix
Sr=diag(Sr); %extract into diagonal matrix
KPSVs=Sc*Sr'; %KP singular values page 50 HNO book
[~, sortKP]=sort(KPSVs(:),'descend');% index into sorted values of the KPs
% Set defaults or find values
if nargin<4 
    tolc=max(Sc)*.01;
else
    tolc=varargin{1};
end
if nargin<5
    tolr=max(Sr)*.01;
    else
    tolr=varargin{2};
end
if nargin<6 
    tolkp=max(KPSVs)*.01;
    else
    tolkp=varargin{3};
end
kc=find(Sc<tolc,1);
kr=find(Sr<tolr,1);
kp=find(sort(KPSVs(:),'descend')<tolkp,1);
if isempty(kc), kr=length(Sc);end
if isempty(kr), kr=length(Sr);end
if isempty(kp), kr=length(KPSVs(:));end
[indSc]=find(Sc>=tolc);
[indSr]=find(Sr>=tolr);
[indKP]=find(KPSVs>=tolkp);
CoeffB=Uc'*B*Ur;
Sctrunc=zeros(size(Sc));
Srtrunc=zeros(size(Sr));
Sctrunc(Sc>tolc)=Sc(Sc>tolc);
Srtrunc(Sr>tolr)=Sr(Sr>tolr);
Kptrunc=Sctrunc*Srtrunc';
Kptrunc(Kptrunc>0)=1./Kptrunc(Kptrunc>0);
X=Vc*(CoeffB.*Kptrunc)*Vr';
KPSVpinv=zeros(size(KPSVs));
KPSVpinv(indKP)=1./KPSVs(indKP);
XKP=Vc*(CoeffB.*KPSVpinv)*Vr';
end