% In Gecko, proteins are just like peptides in ETFL. Each reaction in
% Gecko, represents a pair of enzyme complex and a reaction. So, by
% iterating over the reactions, one can enumerate enzymes. Each enyme is
% different from the others based on its composition, i.e. constituting
% peptides. After finding enzymes, we need to find reactions coupled to
% this enzyme with different kcats. Per each different kcat-enzyme, a new
% element should be defined. It's enough to do it for only one of the
% peptides per each complex.

enzymes = find_enzymes(ecModel);
enzyme_kcat_pairs = couple_enzyme_kcat(ecModel, enzymes);
std_gene_set = convert_to_standard(enzyme_kcat_pairs);
enzyme_kcat_pairs(:,1) = (std_gene_set');
bkwd_kcats = find_different_bkwd(ecModel, enzymes);
std_gene_set = convert_to_standard(bkwd_kcats);
bkwd_kcats(:,1) = (std_gene_set');

function enzymes = find_enzymes(model)
    % model is a Gecko model after manual curation
    enzymes = {0}; % a list of unique complexes
    for i = 1:length(model.rxns)
        %Find set of proteins present in rxn:
        S        = model.S;
        subs_pos = find(S(:,i) < 0);
        prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
        int_pos  = intersect(subs_pos,prot_pos);
        prot_set = cell(size(int_pos)); % the peptides for this complex
        gene_set = cell(size(int_pos)); % the peptides for this complex
        for j = 1:length(int_pos)
            met_name    = model.mets{int_pos(j)};
            prot_set{j} = met_name(6:end);
            gene_set{j} = model.enzGenes(ismember(model.enzymes,prot_set{j}));
        end
        hasMatch = any(cellfun(@isequal, enzymes(:,1), repmat({prot_set}, size(enzymes(:,1)))));
        if ~ hasMatch && length(prot_set) ~= 0
            enzymes{length(enzymes(:,1))+1,1} = prot_set;
            enzymes{length(enzymes(:,1)),2} = gene_set;
        end
    end
    enzymes = enzymes(2:end,:);
end

function enzyme_kcat_pairs = couple_enzyme_kcat(model, enzymes)
    enzyme_kcat_pairs = {0}; % a list of each tuple of unique rxn,enz,kcat
    for i=1:length(enzymes)
       complex = enzymes{i,1};
       gene_ids = enzymes{i,2};
       this_pep_1 = complex{1}; % take the 1st peptide
       this_pep_2 = complex{end}; % take the last to make sure all these peptides are related to this reaction
       [kcats1,rxnIdx1] = getKcat(model,this_pep_1);
       [kcats2,rxnIdx2] = getKcat(model,this_pep_2);
       [rxnIdx, idx, ~] = intersect(rxnIdx1, rxnIdx2);
       kcats = kcats1(idx);
       if length(unique(kcats)) == 1 %kcats are similar
           enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1))+1,1} = gene_ids;
           enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1)),3} = kcats;
       else
           [~,ia,iu] = unique(kcats); 
           original = mode(iu); % the most usual kcat
           original = ia(original);
           exceptional = setdiff(ia,original); % the less usual kcats
           enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1))+1,1} = gene_ids;
           enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1)),3} = kcats(original); % we don't set rxn id for the most usual kcat
           %
           for k=1:length(exceptional)
               rxn_id = model.rxns(rxnIdx(exceptional(k)));
               rxn_id = rxn_id{:};
               enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1))+1,1} = gene_ids;
               enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1)),2} = rxn_id(1:6);
               enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1)),3} = kcats(exceptional(k));
           end
       end
    end
    enzyme_kcat_pairs = enzyme_kcat_pairs(2:end,:);
    for i=1:length(enzyme_kcat_pairs)
       if length(enzyme_kcat_pairs{i,3}) > 1
          helper = enzyme_kcat_pairs{i,3};
          enzyme_kcat_pairs{i,3} = helper(1);
       end
    end
end

function enzyme_kcat_pairs = find_different_bkwd(model, enzymes)
    enzyme_kcat_pairs = {0}; % a list of each tuple of unique rxn,enz,kcat
    for i=1:length(enzymes)
       complex = enzymes{i,1};
       gene_ids = enzymes{i,2};
       this_pep_1 = complex{1}; % take the 1st peptide
       this_pep_2 = complex{end}; % take the last to make sure all these peptides are related to this reaction
       [kcats1,rxnIdx1, rxnName1] = getKcat(model,this_pep_1);
       [kcats2,rxnIdx2, rxnName2] = getKcat(model,this_pep_2);
       [rxnIdx, idx, ~] = intersect(rxnIdx1, rxnIdx2);
       kcats = kcats1(idx);
       rxnName = rxnName1(idx);
       if length(unique(kcats)) == 1 %kcats are similar
           continue
       else
           for k=1:length(rxnName)
               if ~ contains(rxnName(k),'(reversible)') % it's not kcat_bwd
                   continue
               end
               rxn_id = model.rxns(rxnIdx(k));
               rxn_id = rxn_id{:};
               enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1))+1,1} = gene_ids;
               enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1)),2} = rxn_id(1:6);
               enzyme_kcat_pairs{length(enzyme_kcat_pairs(:,1)),3} = kcats(k);
           end
       end
    end
    enzyme_kcat_pairs = enzyme_kcat_pairs(2:end,:);
    for i=1:length(enzyme_kcat_pairs)
       if length(enzyme_kcat_pairs{i,3}) > 1
          helper = enzyme_kcat_pairs{i,3};
          enzyme_kcat_pairs{i,3} = helper(1);
       end
    end
end

function std_gene_set = convert_to_standard(enzyme_kcat_pairs)
    std_gene_set = {};
    gene_set = enzyme_kcat_pairs(:,1);
    for i=1:length(gene_set)
       h_1 = gene_set{i};
       h_2 = '';
       for j=1:length(h_1)
          if isempty(h_2)
              h_2 = strcat(h_2,h_1{j});
          else
              h_2 = strcat(h_2,'//' ,h_1{j});
          end
       end
       std_gene_set(i) = h_2;
    end
end