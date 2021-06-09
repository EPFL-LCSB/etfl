"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Utilities to create a small model from a 1 reaction model in FBA

"""

import cobra
from cobra.test import create_test_model

from pytfa.io import        load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from ..core.memodel import MEModel
from ..core.thermomemodel import ThermoMEModel
from ..core.allocation import add_dna_mass_requirement, \
    add_protein_mass_requirement, add_rna_mass_requirement
from ..optim.config import standard_solver_config


from ..data.ecoli import   get_model, get_thermo_data, get_coupling_dict, \
                        get_mrna_dict, get_rib, get_rnap, get_monomers_dict, \
                        get_nt_sequences, get_ratios, get_neidhardt_data, \
                        get_mrna_metrics, get_enz_metrics, \
                        remove_from_biomass_equation, get_ecoli_gen_stats, \
                        get_essentials, get_lloyd_coupling_dict

from optlang.exceptions import SolverError

import logging


CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'

DEFAULT_SOLVER = GLPK
aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides = get_monomers_dict()
essentials = get_essentials()

def create_fba_model(solver = DEFAULT_SOLVER):

    # g6p_c = cobra.Metabolite(id = 'g6p_c', formula = 'C6H13O9P')
    # f6p_c = cobra.Metabolite(id = 'f6p_c', formula = 'C6H13O9P')
    #
    # g6p_ex = cobra.Reaction(id='EX_g6p')
    # g6p_ex.add_metabolites({g6p_c: -1})
    #
    # f6p_ex = cobra.Reaction(id='EX_g6p')
    # f6p_ex.add_metabolites({f6p_c: -1})
    #
    # pgi = cobra.Reaction(id='PGI')
    # pgi.add_metabolites({g6p_c: -1, f6p_c: 1})
    # pgi.gene_reaction_rule = 'b4025'
    #
    # pgi.lower_bound = -1000
    # g6p_ex.lower_bound = -1000
    # f6p_ex.lower_bound = -1000


    # the_model = cobra.Model()
    the_model = create_test_model('textbook')
    # the_model.add_reactions([g6p_ex, pgi, f6p_ex])

    pi = the_model.metabolites.pi_c
    h = the_model.metabolites.h_c
    h2o = the_model.metabolites.h2o_c

    ppi = cobra.Metabolite(id='ppi_c', formula='P2O7', compartment = 'c')
    ppi_hydrolysis = cobra.Reaction(id='pi_to_ppi')
    ppi_hydrolysis.add_metabolites({ppi:-1, h2o:-1,pi:2, h:2 })

    gdp = cobra.Metabolite(id='gdp_c', formula='C10H15N5O11P2', compartment = 'c' )
    gdp_ex = cobra.Reaction(id='EX_g6p')
    gdp_ex.add_metabolites({gdp: -1})

    the_model.add_reactions([ppi_hydrolysis,gdp_ex])

    add_e_metabolites(the_model)

    for met in ['dttp', 'dctp', 'datp', 'dgtp']:
        the_ex = the_model.reactions.get_by_id('EX_'+met+'_c')
        the_ex.lower_bound = -1
    # add_growth(the_model)

    the_model.solver = solver
    return the_model

def add_e_metabolites(model):
    """
    Adds the metabolites necessary for the expression partof the problem:
    Amino acids, (d)N(M/T)Ps, PPi
    :param model:
    :return:
    """

    new_mets =  list(aa_dict.values()) + \
                list(rna_nucleotides.values()) + \
                list(rna_nucleotides_mp.values()) + \
                list(dna_nucleotides.values())


    for metid in new_mets:

        try:
            the_met = model.metabolites.get_by_id(metid)
        except KeyError:
            the_met = cobra.Metabolite(id = metid, compartment = 'c')
            model.add_metabolites([the_met])

            ex_rxn = cobra.Reaction('EX_{}'.format(metid))

            model.add_reactions([ex_rxn])

            ex_rxn.add_metabolites({the_met:-1})
            ex_rxn.lower_bound = -1
            ex_rxn.upper_bound = 10


def create_etfl_model(has_thermo, has_neidhardt,
                      n_mu_bins = 64,
                      mu_max = 3,
                      optimize = True,
                      solver=DEFAULT_SOLVER,
                      ):
    #------------------------------------------------------------
    # Initialisation
    #------------------------------------------------------------

    # this hack works because we are using the solver switch to update the var
    # names in the solver but really we should not do this
    # TODO: clean up model.sanitize_varnames

    growth_reaction_id = 'Biomass_Ecoli_core'

    vanilla_model = create_fba_model()

    vanilla_model.objective = growth_reaction_id
    fba_sol = vanilla_model.slim_optimize()
    mu_0 = fba_sol
    mu_range = [0, mu_max]
    n_mu_bins = n_mu_bins

    coupling_dict = get_coupling_dict(vanilla_model, mode = 'kmax', atps_name = 'ATPS4r')
    # coupling_dict = get_lloyd_coupling_dict(vanilla_model)



    # Initialize the model

    name = 'small_model_T{:1}E{:1}N{:1}_{}_enz_{}_bins.json'.format(
        has_thermo,
        True,
        has_neidhardt,
        len(coupling_dict),
        n_mu_bins)

    # for k,v in coupling_dict.items():
    #     for enz in v:
    #         enz.kcat_fwd = enz.kcat_bwd = 1e9


    if has_thermo:

        thermo_data, lexicon, compartment_data = get_thermo_data()

        ecoli = ThermoMEModel(thermo_data, model = vanilla_model,
                              growth_reaction = growth_reaction_id,
                              mu_range = mu_range,
                              n_mu_bins = n_mu_bins,
                              name = name,
                              )
    else:
        ecoli = MEModel(model = vanilla_model,
                        growth_reaction = growth_reaction_id,
                        mu_range = mu_range,
                        n_mu_bins = n_mu_bins,
                        name = name,
                        )

    ecoli.name = name
    ecoli.logger.setLevel(logging.WARNING)
    ecoli.sloppy = True

    # Solver settings
    ecoli.solver = solver
    standard_solver_config(ecoli)

    if has_thermo:
        # Annotate the cobra_model
        annotate_from_lexicon(ecoli, lexicon)
        apply_compartment_data(ecoli, compartment_data)

        # TFA conversion
        ecoli.prepare()
        ecoli.convert()#add_displacement = True)


    nt_sequences = get_nt_sequences()
    mrna_dict = get_mrna_dict(ecoli)
    rnap = get_rnap()
    rib = get_rib()


    all_peptides = set([x for enzymes in coupling_dict.values()
                        for enz in enzymes
                        for x in enz.composition])
    prune_to_genes = lambda the_dict:{k:v for k,v in the_dict.items() \
                    if  k in vanilla_model.genes or
                        k in rib.rrna_composition or
                        k in rib.composition or
                        k in rnap.composition or
                        k in all_peptides}

    nt_sequences = prune_to_genes(nt_sequences)
    mrna_dict = prune_to_genes(mrna_dict)

    # Remove nucleotides and amino acids from biomass reaction as they will be
    # taken into account by the expression

    remove_from_biomass_equation(model = ecoli,
                                 nt_dict = rna_nucleotides,
                                 aa_dict = aa_dict,
                                 essentials_dict=essentials)

    ##########################
    ##    MODEL CREATION    ##
    ##########################

    ecoli.add_nucleotide_sequences(nt_sequences)
    ecoli.add_essentials(  essentials=essentials,
                           aa_dict=aa_dict,
                           rna_nucleotides=rna_nucleotides,
                           rna_nucleotides_mp = rna_nucleotides_mp
                           )
    ecoli.add_mrnas(mrna_dict.values())

    ecoli.add_ribosome(rib,0.2)
    ecoli.add_rnap(rnap,0.75)

    ecoli.build_expression()
    ecoli.add_enzymatic_coupling(coupling_dict)
    
    nt_ratios, aa_ratios = get_ratios()
    chromosome_len, gc_ratio = get_ecoli_gen_stats()
    kdeg_mrna, mrna_length_avg  = get_mrna_metrics()
    kdeg_enz,  peptide_length_avg   = get_enz_metrics()
    ecoli.add_dummies(nt_ratios=nt_ratios,
                          mrna_kdeg=kdeg_mrna,
                          mrna_length=mrna_length_avg,
                          aa_ratios=aa_ratios,
                          enzyme_kdeg=kdeg_enz,
                          peptide_length=peptide_length_avg)

    if has_neidhardt:
        neidhardt_mu, neidhardt_rrel, neidhardt_prel, neidhardt_drel = get_neidhardt_data()

        
        add_protein_mass_requirement(ecoli,neidhardt_mu, neidhardt_prel)
        add_rna_mass_requirement(ecoli,neidhardt_mu, neidhardt_rrel)
        add_dna_mass_requirement(ecoli,mu_values=neidhardt_mu,
                                       dna_rel=neidhardt_drel,
                                       gc_ratio=gc_ratio,
                                       chromosome_len=chromosome_len,
                                       dna_dict=dna_nucleotides)

    # Need to put after, because dummy has to be taken into account if used.
    ecoli.populate_expression()
    ecoli.add_trna_mass_balances()


    ecoli.print_info()

    need_relax = False
    
    ### a problem with the solver, which can be temporarily solved
    ecoli.constraints.MB_b3855.set_linear_coefficients({ecoli.variables.b3855_degradation:-1e-5})

    ecoli.repair()

    if optimize:
        try:
            ecoli.optimize()

            print('Objective            : {}'.format(ecoli.solution.objective_value))
            print(' - Glucose uptake    : {}'.format(ecoli.reactions.EX_glc__D_e.flux))
            print(' - Growth            : {}'.format(ecoli.growth_reaction.flux))
            print(' - Ribosomes produced: {}'.format(ecoli.ribosome.X))
            print(' - RNAP produced: {}'.format(ecoli.rnap.X))

        except (AttributeError, SolverError):
            pass

    return ecoli


def create_simple_dynamic_model():
    from cobra import Reaction, Metabolite, Model

    from etfl.core import Enzyme, MEModel
    from etfl.optim.constraints import TotalCapacity, ForwardCatalyticConstraint, EnzymeMassBalance

    glc  = Metabolite('glucose')
    lcts = Metabolite('lactose')
    mint = Metabolite('intermediate')
    bio  = Metabolite('biomass')

    glc_yield = 1
    lcts_yield = 0.95

    n_glc_C  = 6 #  number of carbons in glucose
    n_lcts_C = 12 # number of carbons in lactose
    n_mint_C = 3 # number of carbons in metabolic intermediate (pyr for example)

    ex_glc  = Reaction('EX_glc',  lower_bound = -10)
    ex_glc. add_metabolites({glc:-1})
    ex_lcts = Reaction('EX_lcts', lower_bound = -n_glc_C/n_lcts_C * 10)
    ex_lcts.add_metabolites({lcts:-1})

    glc2mint  = Reaction('glc2mint')
    glc2mint.add_metabolites({glc :-1,
                              mint:glc_yield  * n_glc_C /n_mint_C})
    lcts2mint = Reaction('lcts2mint')
    lcts2mint.add_metabolites({lcts:-1,
                               mint:lcts_yield * n_lcts_C/n_mint_C})
    mint2bio  = Reaction('mint2bio')
    mint2bio.add_metabolites({mint:-1,
                              #bio:1
                              })

    model = Model('smol_model')
    model = MEModel(model)

    model.add_reactions([glc2mint, lcts2mint, mint2bio, ex_glc, ex_lcts])

    kcat = 10
    kdeg = 1
    delta_v_ub = 0.5
    Eg = Enzyme(id='Eg', kcat=kcat, kdeg=kdeg, composition = {})
    Eg._scaling_factor = 1
    Eg.complexation = Reaction('Eg_cplx', upper_bound = delta_v_ub)
    Eg.degradation  = Reaction('Eg_deg' , upper_bound = delta_v_ub)
    Eg.complexation.scaling_factor = 1
    El = Enzyme(id='El', kcat=kcat, kdeg=kdeg, composition = {})
    El._scaling_factor = 1
    El.complexation = Reaction('El_cplx', upper_bound = delta_v_ub)
    El.degradation  = Reaction('El_deg' , upper_bound = delta_v_ub)
    El.complexation.scaling_factor = 1

    model.add_enzymes([Eg,El], prep = False) # We don't want synthesis et al.
    model.add_reactions([Eg.complexation, El.complexation, Eg.degradation, El.degradation])
    Eg.init_variable()
    El.init_variable()

    # coupling_dict = {glc2mint.id:[Eg], lcts2mint:[El]}

    model.add_constraint(kind=ForwardCatalyticConstraint,
                         hook=glc2mint,
                         expr=glc2mint.forward_variable - Eg.kcat_fwd*Eg.variable,
                         ub=0)
    model.add_constraint(kind=ForwardCatalyticConstraint,
                         hook=lcts2mint,
                         expr=lcts2mint.forward_variable - El.kcat_fwd*El.variable,
                         ub=0)
    model.add_constraint(kind=TotalCapacity,
                         hook=model,
                         id_='total_E',
                         expr = Eg.variable + 3*El.variable,
                         lb=0.1,
                         ub=0.1)
    model.add_constraint(kind=EnzymeMassBalance,
                         hook=Eg,
                         expr=Eg.complexation.forward_variable - Eg.degradation.forward_variable,
                         ub=0,lb=0)
    model.add_constraint(kind=EnzymeMassBalance,
                         hook=El,
                         expr=El.complexation.forward_variable - El.degradation.forward_variable,
                         ub=0,lb=0)

    model.objective = mint2bio
    model.growth_reaction = mint2bio.id

    model.repair()
    model.optimize()

    def get_yield():
        return model.solution.objective_value / (lcts2mint.flux * n_lcts_C + glc2mint.flux * n_glc_C)

    # model.mu_bins = [[x/10,(max(0,x/10-0.05), x/10+0.05)] for x in range(0,100)]
    model.mu_bins = [[0,(max(0,x/10-0.05), x/10+0.05)] for x in range(0,100)]

    print('Yield:', get_yield())
    return model

if __name__ == '__main__':
    model = create_etfl_model(False, False)