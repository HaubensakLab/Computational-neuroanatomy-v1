% Copyright (C) 2016 VRVis.
% All rights reserved.
% Contact: VRVis Forschungs-GmbH (office@vrvis.at)
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the VRVis Forschungs-GmbH.
% 4. Neither the name of the VRVis Forschungs-GmbH nor the
%    names of its contributors may be used to endorse or promote products
%    derived from this software without specific prior written permission.
%

%The first part generates atlasRegions. atlasRegions is a vector with a length of 3702, where every element contains the brain region id

urlwrite('http://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,products[abbreviation$eq%27HumanMA%27],rma::include,specimens[parent_id$eqnull](alignment3d),rma::options[only$eq%27donors.id,donors.name,specimens.id%27]','..//storage//tempDonors.json'); 
fid = fopen('..//storage//tempDonors.json');
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
donorJson = parse_json(str);


urlwrite('http://api.brain-map.org/api/v2/data/query.xml?num_rows=10&criteria=model::Probe,rma::criteria,[probe_type$eq%27DNA%27],products[abbreviation$eq%27HumanMA%27],rma::options[only$eq%27probes.id%27]','..//storage//temp1ProbeInDB.xml');
root = xmlread('..//storage//temp1ProbeInDB.xml'); 

elems = root.getElementsByTagName('id');

urlwrite(strcat('http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$eq',char(elems.item(0).getTextContent),']'),'..//storage//tempSamples.json');

fid = fopen('..//storage//tempSamples.json');
raw = fread(fid,inf);
str = char(raw');
fclose(fid);

indizes = strfind(str,'mri');
coordinatesPos = 0;

coordinates = zeros(size(indizes,2),3);
indexX = zeros(size(indizes,2),1);
indexY = zeros(size(indizes,2),1);
indexZ = zeros(size(indizes,2),1);
atlasRegions = zeros(size(indizes,2),1);

for j = indizes
    coordinatesPos = coordinatesPos +1;
    actpos = j;
    while not(str(actpos) == ']')
        actpos=actpos+1;
    end
    coordinates(coordinatesPos,:) =str2num(str((j+6):(actpos-1)));
end

indizes = strfind(str,'donor');
donorPos = 0;
donorIDs = zeros(size(indizes,2),1);

for j = indizes
    donorPos = donorPos +1;
    actpos = j;
    while not(str(actpos) == ',')
        actpos=actpos+1;
    end
    donorID = str2num(str((j+13):(actpos-1)));
    donorIDs(donorPos)=donorID;
    foundDonor = false;
    
   
    for k = 1:size(donorJson{1}.msg,2)
        if donorJson{1}.msg{k}.id == donorID
          
           x = donorJson{1}.msg{k}.specimens{1}.alignment3d;
           T = [ x.tvr_00 x.tvr_01 x.tvr_02 x.tvr_09;
                 x.tvr_03 x.tvr_04 x.tvr_05 x.tvr_10;
                 x.tvr_06 x.tvr_07 x.tvr_08 x.tvr_11;
                 0        0        0        1];
           newCoords = T*[coordinates(donorPos,:),1]';
           coordinates(donorPos,:) = newCoords(1:3)'; 
          
           foundDonor = true;
        end
    end
    
    
    if foundDonor == false
        disp('Donor not found!');
    end

end


indizes = strfind(str,'"structure"');
structuresPos = 0;

coordinatesPos = 0;

for j = indizes
    structuresPos = structuresPos +1;
    coordinatesPos = coordinatesPos +1;
    actpos = j;
    while not(str(actpos) == ',')
        actpos=actpos+1;
    end
    indexX(coordinatesPos)=coordinates(coordinatesPos,1);
    indexY(coordinatesPos)=coordinates(coordinatesPos,2);
    indexZ(coordinatesPos)=coordinates(coordinatesPos,3);
    atlasRegions(coordinatesPos)=str2num(str((j+18):(actpos-1)));
 
end
 
save('../storage//atlasRegionsHuman.mat','atlasRegions');
save('../storage//atlasRegionsHuman.mat','indexX','-append');
save('../storage//atlasRegionsHuman.mat','indexY','-append');
save('../storage//atlasRegionsHuman.mat','indexZ','-append');

disp('Atlas generation complete!')

%The second part downloads the gene-expression for every gene in the AHBA ony biopsy-site level (corresponds atlasRegions
urlwrite('http://api.brain-map.org/api/v2/data/query.json?num_rows=20000&criteria=model::Probe,rma::criteria,[probe_type$eq%27DNA%27],products[abbreviation$eq%27HumanMA%27]','..//storage//allGeneIdsInDB1.json');
urlwrite('http://api.brain-map.org/api/v2/data/query.json?start_row=20000&num_rows=20000&criteria=model::Probe,rma::criteria,[probe_type$eq%27DNA%27],products[abbreviation$eq%27HumanMA%27]','..//storage//allGeneIdsInDB2.json');
urlwrite('http://api.brain-map.org/api/v2/data/query.json?start_row=40000&num_rows=20000&criteria=model::Probe,rma::criteria,[probe_type$eq%27DNA%27],products[abbreviation$eq%27HumanMA%27]','..//storage//allGeneIdsInDB3.json');

fid = fopen('..//storage//allGeneIdsInDB1.json');
raw = fread(fid,inf);
str = char(raw');
fclose(fid);

fid = fopen('..//storage//allGeneIdsInDB2.json');
raw = fread(fid,inf);
str = [str,char(raw')];
fclose(fid);

fid = fopen('..//storage//allGeneIdsInDB3.json');
raw = fread(fid,inf);
str = [str,char(raw')];
fclose(fid);

indizes = strfind(str,'gene_id');
disp(sprintf('%d Entrez-IDs downloaded!',size(indizes,2)))
genePos = 0;
geneIDsOfRows = zeros(size(indizes,2),1);

for j = indizes
    genePos = genePos +1;
    actpos = j;
    while not(str(actpos) == ',')
        actpos=actpos+1;
    end
    geneIDsOfRows(genePos) =str2num(str((j+9):(actpos-1)));
end

geneIDsOfRows = unique(geneIDsOfRows);
entrezIDsOfRows = zeros(size(geneIDsOfRows,1),1);

actGeneIDIndex = 0;

expressionMatrix =sparse(length(atlasRegions),length(geneIDsOfRows));

for actGeneID = geneIDsOfRows'
    actGeneIDIndex = actGeneIDIndex + 1;
    
    disp(sprintf('Load %d/%d',actGeneIDIndex,length(geneIDsOfRows)))
    probes = 0;
    errorfound=true;
            while errorfound
                try
               
                    urlwrite(strcat('http://api.brain-map.org/api/v2/data/query.json?criteria=model::Probe,rma::criteria,[probe_type$eq%27DNA%27],products[abbreviation$eq%27HumanMA%27],gene[gene_id$eq%27',sprintf('%d',actGeneID),'%27],rma::options[only$eq%27probes.id%27]'),'..//storage//tempProbe.json','Timeout',10);
                       fid = fopen('..//storage//tempProbe.json');
                        raw = fread(fid,inf);
                        str = char(raw');
                        fclose(fid);

                        str = str(83:end);
                        indizes = strfind(str,'id');

                        probePos = 0;
                        probes = zeros(size(indizes,2),1);

                        for j = indizes
                            probePos = probePos +1;
                            actpos = j;
                            while not(str(actpos) == '}')
                                actpos=actpos+1;
                            end
                            probes(probePos) =str2num(str((j+4):(actpos-1)));
                        end
    
                      errorfound=false;
                catch err
                      errorfound=true;
                      disp('try again');
                end
            end
   
 
    for probe = probes'
        errorfound=true;
        while errorfound
                try
                     urlwrite(strcat('http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$eq%27',sprintf('%d',probe),'%27]'),'..//storage//tempSample.json','Timeout',10);
                        fid = fopen('..//storage//tempSample.json');
                        raw = fread(fid,inf);
                        str = char(raw');
                        fclose(fid);

                        indizes = strfind(str,'entrez-id');
                        probePos = 0;
                        for j = indizes
                            probePos = probePos +1;
                            actpos = j;
                            while not(str(actpos) == ',')
                                actpos=actpos+1;
                            end
                            entrezIDsOfRows(actGeneIDIndex)=str2num(str((j+11):(actpos-1)));
                        end
                         indizes = strfind(str,'z-score');
                        probePos = 0;
                        for j = indizes
                            probePos = probePos +1;
                            actpos = j;
                            while not(str(actpos) == ']')
                                actpos=actpos+1;
                            end
                            a=str2num(strrep(str((j+10):(actpos-1)),'"',''))';
                        end
    
                                    errorfound=false;
                      catch err
                                  errorfound=true;
                                  disp('try again');
                            end
                        end
                        indizes = strfind(str,'z-score');
                        probePos = 0;
                        for j = indizes
                            probePos = probePos +1;
                            actpos = j;
                            while not(str(actpos) == ']')
                                actpos=actpos+1;
                            end
                            expressionMatrix(:,actGeneIDIndex) = expressionMatrix(:,actGeneIDIndex) + str2num(strrep(str((j+10):(actpos-1)),'"',''))';
                        end
            
    end
     expressionMatrix(:,actGeneIDIndex) = expressionMatrix(:,actGeneIDIndex)/length(probes);
end

regionIDs = atlasRegions;
save(sprintf('../storage/all_genes_expression_Human.mat'),'expressionMatrix');
save(sprintf('../storage/all_genes_expression_Human.mat'),'entrezIDsOfRows','-append');
save(sprintf('../storage/all_genes_expression_Human.mat'),'regionIDs','-append');


disp('Computation complete!')


