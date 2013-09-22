function [A,B] = JointDistAlt(AllCa,AllCb)

% this is why I hate combinatorics
% new, shinier and hopefully more reasonable solution 
% (August 22, 2012)


AllC = AllCa+AllCb;
N = size(AllC,2); %number of series in collection
L = size(AllC,1); %number of samples per series


AC = sum(AllC,2); %total activity distribution
p = sum(AC)/(N*L);
NAct = hist(AC,0:N);
ACa = sum(AllCa,2);
pa = sum(ACa)/(N*L); %estimate of likelyhood of event a
ACb = sum(AllCb,2);
pb = sum(ACb)/(N*L); %estimate of likelyhood of event b
    pab = pb/p; % if there is an action, the likelyhood of it being b


Model = zeros(N+1);
altModel = zeros(N+1);

for s = 0:N %activity sum
    set = zeros(s+1,1);
   for r = 0:s
      set(r+1) = nchoosek(s,r)*(pa^(s-r))*(pb^r)/p^s;
   end
   set = NAct(s+1)*set/sum(set);
   altModel(1:s+1,s+1) = set;
   for r = 0:s
      Model(r+1,s+1-r) = set(r+1);
   end
end


 B = Model;
A=zeros(N+1);
for j=0:N
    A(j+1,:)=hist(ACb(ACa==j),0:N); %make sure this isn't backwards
end

 
 
% mapping from s+1 values evenly to N+1 spots
mapper=cell(N+1,1);

mapper{1} = round(0.5*(N))+1;

for s = 1:N
    mapper{s+1} = round(N*(0:s)/s)+1;
end
 
% switch A and B to this form

Aalt = zeros(N+1);
Balt = Aalt;

for s = 0:N
    v = hist(ACb(AC==s),0:s);
    if ~isempty(v)
        Aalt(mapper{s+1},s+1)= v;
    end
   Balt(mapper{s+1},s+1)=altModel(1:s+1,s+1);
end

A = Aalt;
B = Balt;
 
