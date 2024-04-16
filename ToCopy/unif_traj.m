%% supporting function
function traj = unif_traj(RO,NProj)
x=-RO/2:RO/2-1;

phin1=0;
for ii = 1:NProj
    hn = -1 + 2*(ii-1)/(NProj-1);
    theta = acos(hn);
    if (ii == NProj || ii == 1)
        phi = 0;
    else
        phi = mod(phin1+3.6/sqrt(NProj*(1-hn^2)),2*pi);
    end
    phin1=phi;
    
    xA = cos(phi) * sin(theta);
    yA = sin(phi)*sin(theta);
    zA = cos(theta);
  
    traj(1,:,ii) = x*xA;
    traj(2,:,ii) = x*yA;
    traj(3,:,ii) = x*zA;
end
end