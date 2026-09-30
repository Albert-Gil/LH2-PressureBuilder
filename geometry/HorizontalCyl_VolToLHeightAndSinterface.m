function [H, interfaceArea, interfacePerimeter] = HorizontalCyl_VolToLHeightAndSinterface(V, R, L)
%HORIZONTALCYL_VOLTOLHEIGHTANDSINTERFACE Horizontal circular-cylinder geometry.
% SI units. Scalar R,L; V may be an array. Out-of-range volumes are errors.
validateattributes(V,{'numeric'},{'real','finite','nonnegative'});
validateattributes(R,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(L,{'numeric'},{'scalar','real','finite','positive'});
capacity=pi*R^2*L;
if any(V(:)>capacity*(1+1e-12)), error('vap:Volume','Liquid volume exceeds tank capacity.'); end
fraction=min(V/capacity,1);
H=zeros(size(V)); interfaceArea=H; interfacePerimeter=H;
for k=1:numel(V)
    f=fraction(k);
    if f==0, continue; end
    if f==1, H(k)=2*R; continue; end
    % Solve the minor segment angle on a fixed, closed, scale-free bracket.
    % A small-angle series avoids cancellation near an empty/full tank.
    target=pi*min(f,1-f);
    theta=fzero(@(a) segment(a)-target,[0,pi/2],optimset('TolX',1e-14));
    minorHeight=2*R*sin(theta/2)^2;
    if f<=.5, H(k)=minorHeight; else, H(k)=2*R-minorHeight; end
    chord=2*R*sin(theta);
    interfaceArea(k)=chord*L;
    interfacePerimeter(k)=2*(chord+L);
end
end
function a=segment(theta)
if abs(theta)<1e-3
    a=(2/3)*theta^3-(2/15)*theta^5+(4/315)*theta^7;
else
    a=theta-sin(theta)*cos(theta);
end
end
