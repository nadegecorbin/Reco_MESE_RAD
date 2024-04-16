% write basis file 

function write_basis_shared(TE,nT2,filename,nComp,nB1,comp, loga)

if nargin<5
nB1=80;
end

if nargin<6
    comp=0;
end

if nargin<7
    loga=1
end

% range T2 and B1
if loga 
T2=exp(log(10)*linspace(0.7,2.3,nT2)); %3.3
else
T2=linspace(5,200,nT2)
end
%T2=5:10:2300;
B1=linspace(0.7,1,nB1);


ESP=TE(2)-TE(1);
T1=1000%800:100:2000;
nT1=length(T1);
TR=1000;

% get signal curves for every T2 and B1
parfor t2=1:length(T2)
   for b1=1:nB1
       for t1=1:nT1
        XBase(:,t2,b1,t1)=MonoExpEPG(TE,T2(t2),B1(b1),ESP,T1(t1),TR,comp); %last argument is 1 for complex data and 0 for magnitude data
        end
    end
end

% svd
[U S V]=bart('svd',reshape(XBase,length(TE),length(T2)*nB1*length(T1))); %%input=U*diag(S)*V

writecfl(filename,permute(U(:,1:nComp),[3 4 5 6 7 1 2]));

end 



function Sig=MonoExpEPG(t,T2,B1,ESP,T1,TR,comp)
        
        nEcho=length(t);
        theta=[B1*pi/2 B1*repmat([pi],1,nEcho)]; 

        Sig=zeros(size(t));
        
        [F0,Fn,Zn,F] = EPG_MESE(theta,ESP,T1,T2,2,TR);
        if nargin<7
        Sig=abs(F0(end-nEcho+1:end));
        else
            if comp==1
             Sig=(F0(end-nEcho+1:end));
            else
             Sig=abs(F0(end-nEcho+1:end));
            end
        end
  end



  