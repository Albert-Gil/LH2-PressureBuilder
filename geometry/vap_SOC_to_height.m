function H = vap_SOC_to_height(SOC, D, Length, vol_m3)
%VAP_SOC_TO_HEIGHT Volume fill fraction -> height using the local helper.
validateattributes(SOC,{'numeric'},{'real','finite','>=',0,'<=',1});
validateattributes(D,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(Length,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(vol_m3,{'numeric'},{'scalar','real','finite','positive'});
capacity=pi*D^2*Length/4;
if abs(capacity-vol_m3)>1e-10*capacity
    error('vap:Geometry','vol_m3 must match the circular-cylinder geometry.');
end
[H,~,~]=HorizontalCyl_VolToLHeightAndSinterface(SOC*vol_m3,D/2,Length);
end
