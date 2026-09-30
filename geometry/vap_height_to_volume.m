function V = vap_height_to_volume(H, D, Length)
%VAP_HEIGHT_TO_VOLUME Stable circular-segment volume, including the ends.
validateattributes(D,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(Length,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(H,{'numeric'},{'real','finite','>=',0,'<=',D});
R=D/2; minor=min(H,D-H);
theta=2*asin(sqrt(minor/(2*R)));
area=theta-sin(theta).*cos(theta);
small=theta<1e-3;
a=theta(small);
area(small)=(2/3)*a.^3-(2/15)*a.^5+(4/315)*a.^7;
V=R^2*Length*area;
upper=H>R;
V(upper)=pi*R^2*Length-V(upper);
end
