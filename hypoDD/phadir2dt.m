function mct = phadir2dt(phadir,jE)
% function mct = phadir2dt(phadir,jE)
%
% 2020-05-05
% This function makes a list of "catalog" differential times
% from a list of picks. The differential times (DT) are T1a-T1b, where 
% "1" is a station number and a/b are events. 
%
% INPUTS
%
%     phadir = [EV,STA,PHA,T]
%        jE  = list of event pairs to use [EV1,EV2]
%
% OUTPUTS
%
%     mct = [EVa,EVb,STA,PHA,DT,(CC)]
%            note that the CC column is just zero in this case.


mct = zeros(1e7,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- First deal with P-waves
A = phadir(phadir(:,3)==1,:);

usta = unique(A(:,2));
jj   = 0;

for jss = usta'
    
    Asta = A(A(:,2)==jss,:);
    Ne   = size(Asta,1);
    
    if Ne < 2
        continue
    end
   
    [sw1,j1] = ismember(jE(:,1),Asta(:,1));
    [sw2,j2] = ismember(jE(:,2),Asta(:,1));
    jb = find(sw1.*sw2);
    N = size(jb,1);
    
    if N 
        dt = Asta(j1(jb),4)-Asta(j2(jb),4);
        mct(jj+[1:N]',:) = [Asta(j1(jb),1),Asta(j2(jb),1),repmat([jss,1],N,1),dt,ones(N,1)];
        jj = jj + N;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Repeat for S-waves
A = phadir(phadir(:,3)==2,:);

usta = unique(A(:,2));

for jss = usta'
    
    Asta = A(A(:,2)==jss,:);
    Ne   = size(Asta,1);
    
    if Ne < 2
        continue
    end
   
    [sw1,j1] = ismember(jE(:,1),Asta(:,1));
    [sw2,j2] = ismember(jE(:,2),Asta(:,1));
    jb = find(sw1.*sw2);
    N = size(jb,1);
    
    if N 
        dt = Asta(j1(jb),4)-Asta(j2(jb),4);
        mct(jj+[1:N]',:) = [Asta(j1(jb),1),Asta(j2(jb),1),repmat([jss,2],N,1),dt,ones(N,1)];
        jj = jj + N;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Sort and cull

mct(jj+1:end,:) = [];
mct = sortrows(mct);


