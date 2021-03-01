import cobra
from cobra.flux_analysis import single_gene_deletion
from gene_knockout import single_gene_knockout
import pandas as pd


# model = cobra.io.read_sbml_model('../reproducibility/yeast-GEM-latest_version/ModelFiles/xml/yeastGEM.xml')
model = cobra.io.load_matlab_model('../thermo_curation/yeast8_thermo_curated.mat')


# setting the medium
def complete_medium(model):
    constrainedUptake = ['r_1604','r_1639','r_1873','r_1879','r_1880',
            'r_1881','r_1671','r_1883','r_1757','r_1891','r_1889','r_1810',
            'r_1993','r_1893','r_1897','r_1947','r_1899','r_1900','r_1902',
            'r_1967',
            #'r_1918', Y7 doesn't have a linoleic acid exchange rxn, which
            #was substituted for octadecenoate and octadecynoate for Y5 and Y6,
            #so leave out.
            'r_1903','r_1548','r_1904','r_2028',
            'r_2038','r_1906','r_2067','r_1911','r_1912','r_1913','r_2090',
            'r_1914','r_2106']   # potassium exchange
    glucoseExchange = ['r_1714']
    unconstrainedUptake = ['r_1672','r_1654', # ammonium exchange
                        'r_1992', # oxygen exchange
                        'r_2005', # phosphate exchange
                        'r_2060', # sulphate exchange
                        'r_1861', # iron exchange, for test of expanded biomass def
                        'r_1832', # hydrogen exchange
                        'r_2100', # water exchange
                        'r_4593', # chloride exchange
                        'r_4595', # Mn(2+) exchange
                        'r_4596', # Zn(2+) exchange
                        'r_4597', # Mg(2+) exchange
                        'r_2049', # sodium exchange
                        'r_4594', # Cu(2+) exchange
                        'r_4600', # Ca(2+) exchange
                        'r_2020' ]
    for rxn_id in constrainedUptake:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = -0.5
        
    for rxn_id in glucoseExchange:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = -20
        
    for rxn_id in unconstrainedUptake:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = -1000

if __name__ == '__main__':
    complete_medium(model)
    # single_gene_deletion function for all genes is not reproducible!
    # model.solver.problem.Params.FeasibilityTol = 1e-7
    deletion_results = dict()
    for gene in model.genes:
        deletion_results[gene.id] = single_gene_knockout(model,[gene.id])
        # model.solver.problem.reset()  #Discard any solution information.  The next optimize() call starts from scratch.
    
    helper = [v for _,v in deletion_results.items()]
    deletion_results = pd.concat(helper)
    
    results = deletion_results['growth']
    results.index = [list(x)[0] for x in results.index]
    
    results.to_csv('outputs/gene_essentiality_FBA.csv')