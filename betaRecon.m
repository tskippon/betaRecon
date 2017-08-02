function [newgrains newebsd]=betaRecon(ebsd,grains)
tic

for i=[unique(ebsd(ebsd.phase>0).phase)]'
    if strcmp(ebsd(ebsd.phase==i).CS.lattice,'hexagonal')
        alpha=ebsd.phase==i;
    elseif strcmp(ebsd(ebsd.phase==i).CS.lattice,'cubic')
        beta=ebsd.phase==i;
    end
end


%Set up parallel directions for alpha/beta relationship
a1=Miller(0,0,0,2,ebsd(alpha).CS);
b1=Miller(1,1,0,ebsd(beta).CS);
a2=Miller(1,1,-2,0,ebsd(alpha).CS);
b2=Miller(1,1,1,ebsd(beta).CS);


%set misorientation threshold (in degrees) for determining if two alpha grain are from the same beta
%grain
threshold=2;

%get lsit of pairs of neighboring grains
[counts,pairs]=grains.neighbors;

%preallocate memory for reconstructed grains
reconstructedGrains=grains;

%set up mapping for transformation from alpha to beta
ba=symmetrise(orientation('map',b1,a1,b2,a2,ebsd(beta).CS,ebsd(alpha).CS));


%get phase labels for alpha and beta data in ebsd
betaPhaseNum=unique(ebsd(beta).phase);
alphaPhaseNum=unique(ebsd(alpha).phase);


%look through list of pairs and for each pair, check if they can both be
%transformed into beta grains with misorientation between them less than
%threshold
fprintf('Reconstructing beta grains from alpha grain pairs.');

g1Votes=zeros(length(pairs),length(ba));
g2Votes=zeros(length(pairs),length(ba));

for i=1:length(pairs)
    if(mod(i,5000)==0)
        fprintf('\nFinished %i of %i grain pairs...\n',i,length(pairs));
        toc
    end
    %If both grains are alpha, then transform both to beta
    if reconstructedGrains(pairs(i,1)).phase==alphaPhaseNum && reconstructedGrains(pairs(i,2)).phase==alphaPhaseNum
        
        %g1 and g2 are grain 1 and grain 2
        %here we transform both of them using different alpha to
        %beta variants and take only the unique results
        g1=unique(repmat(grains(pairs(i,1)).meanOrientation,length(ba),1).*ba);
        g2=unique(repmat(grains(pairs(i,2)).meanOrientation,length(ba),1).*ba);
        
        
        %compare all unique g1 transformations to all unique g2 transformations and
        %find the minimum misorientation angle we can achieve
        g1 = g1( ceil( (1:length(g2)*length(g1))/length(g2) ) );
        g2 = repmat(g2,length(g2),1);
        angles=angle(g1,g2)/degree;
        minAngle=min(angles);
        
        %If that misorientation angle is < threshold, then save the
        %versions of g1 and g2 that gave that minimum angle to
        %reconstructedGrains
        if minAngle<threshold
            voteTemp=ceil(find(angles==minAngle)/length(ba));
            weight=sum(grains.boundary(grains.boundary.hasGrain(pairs(i,1),pairs(i,2))).segLength);
            g1Votes(i,voteTemp)=g1Votes(i,voteTemp)+weight;
            voteTemp=mod(find(angles==minAngle),length(ba));
            voteTemp(voteTemp==0)=length(ba);
            g2Votes(i,voteTemp)=g2Votes(i,voteTemp)+weight;
            %reconstructedGrains(pairs(i,1:2)).phase=betaPhaseNum;
            %reconstructedGrains(pairs(i,1:2)).meanOrientation=[g1(find(angles==min(angles),1)); g2(find(angles==min(angles),1))];
        end
        
    %If first grain is beta and second grain is alpha, then transform
    %second grain to beta using all variants and compare orientations
    %to the first grain
    elseif reconstructedGrains(pairs(i,1)).phase==betaPhaseNum && reconstructedGrains(pairs(i,2)).phase==alphaPhaseNum
        g1=reconstructedGrains(pairs(i,1)).meanOrientation;
        g2=unique(repmat(grains(pairs(i,2)).meanOrientation,length(ba),1).*ba);
        angles=angle(g1,g2)/degree;
        minAngle=min(angles);
        if minAngle<threshold
%             reconstructedGrains(pairs(i,2)).phase=betaPhaseNum;
%             reconstructedGrains(pairs(i,2)).meanOrientation=g2(find(angles==min(angles),1));
            voteTemp=mod(find(angles==minAngle),length(ba));
            voteTemp(voteTemp==0)=length(ba);
            weight=sum(grains.boundary(grains.boundary.hasGrain(pairs(i,1),pairs(i,2))).segLength);
            g2Votes(i,voteTemp)=g2Votes(i,voteTemp)+weight;
        end
    %If first grain is alpha and second is beta, then transform first grain
    %to beta using all 6 variants and compare orientations to second grain
    elseif reconstructedGrains(pairs(i,1)).phase==alphaPhaseNum && reconstructedGrains(pairs(i,2)).phase==betaPhaseNum
        
        g1=unique(repmat(grains(pairs(i,1)).meanOrientation,length(ba),1).*ba);
        g2=reconstructedGrains(pairs(i,2)).meanOrientation;
        
        angles=angle(g1,g2)/degree;
        minAngle=min(angles);
        if minAngle<threshold
%             reconstructedGrains(pairs(i,1)).phase=betaPhaseNum;
%             reconstructedGrains(pairs(i,1)).meanOrientation=g1(find(angles==min(angles),1));
            voteTemp=ceil(find(angles==minAngle)/length(ba));
            weight=sum(grains.boundary(grains.boundary.hasGrain(pairs(i,1),pairs(i,2))).segLength);
            g1Votes(i,voteTemp)=g1Votes(i,voteTemp)+weight;

        end
     end
end


%now we have reconstructedGrains, in which the alpha grains have been
%transformed into their parent beta grains where possible

reconstructedGrains.prop.vtype=zeros(length(reconstructedGrains),1);
for i=unique(pairs)'
    votes=sum([g1Votes(pairs(:,1)==i,:); g2Votes(pairs(:,2)==i,:)]);
    winner=find(votes==max(votes) & votes>0);
    if ~isempty(winner)
        if length(winner>1)
            winner=winner(1);
        end
        reconstructedGrains(i).vtype=winner;
        reconstructedGrains(i).phase=betaPhaseNum;
        reconstructedGrains(i).meanOrientation=grains(i).meanOrientation*ba(winner);
    end
end


%temp is the size of the indexed part of ebsd and has data from
%reconstructedGrains.  It will let us quickly transfer reconstructed grain
%info into ebsd data
  temp=reconstructedGrains(ebsd('indexed').grainId);
  
%create newebsd, whiCSch is equal to the indexed parts of ebsd
  newebsd=ebsd('indexed');
  
%set phase for each point in newebsd to be equal to the phase of the
%reconstructedGrain it is located in.  This means that points in
%transformed grains will now have their phase set to beta
  newebsd.phase=temp.phase;
  
%set orientation data of points to be equal to the mean orientation of the
%reconstructedGrain it is located in.  This means that we lose any subgrain
%orientation variation and all points in a grain will have the same
%orientation, but we probably don't care.
  newebsd(temp.phase==betaPhaseNum).orientations=temp(temp.phase==betaPhaseNum).meanOrientation;
  
%finally, we make MTEX run grain identification on newebsd.  This will let
%it find the grain boundaries between the transformed beta grains
  [newgrains,newebsd.grainId,newebsd.mis2mean]=calcGrains(newebsd('indexed'),'angle',5*degree);

fprintf('\nGrain Reconstruction Complete\n',i,length(pairs));
toc
end