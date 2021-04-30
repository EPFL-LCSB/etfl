% add inner mitochondrial space
model.mets(end+1)        = {'h_i'};
model.metNames(end+1)    = {'H+ [mitochondrion inner membrane]'};
model.metFormulas(end+1) = {'H'};
model.metCharges(end+1)  = 1;
model.S(end+1,:) = zeros(1,length(model.rxns));
model.csense(end+1) = 'E';
model.b(end+1) = 0;
 
ind_h_i = length(model.mets);
ind_h_c = find(ismember(model.metNames,'H+ [cytoplasm]')); 
 
% ATP synthase
stoich_coef = model.S(ind_h_c,find(ismember(model.rxns,'r_0226')));
model.S(ind_h_c,find(ismember(model.rxns,'r_0226'))) = 0;
model.S(ind_h_i,find(ismember(model.rxns,'r_0226'))) = stoich_coef;
 
% ferrocytochrome-c:oxygen oxidoreductase
stoich_coef = model.S(ind_h_c,find(ismember(model.rxns,'r_0438')));
model.S(ind_h_c,find(ismember(model.rxns,'r_0438'))) = 0;
model.S(ind_h_i,find(ismember(model.rxns,'r_0438'))) = stoich_coef;
 
% ubiquinol:ferricytochrome c reductase
stoich_coef = model.S(ind_h_c,find(ismember(model.rxns,'r_0439')));
model.S(ind_h_c,find(ismember(model.rxns,'r_0439'))) = 0;
model.S(ind_h_i,find(ismember(model.rxns,'r_0439'))) = stoich_coef;