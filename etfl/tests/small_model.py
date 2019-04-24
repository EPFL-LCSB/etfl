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

solver = GLPK
aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides = get_monomers_dict()
essentials = get_essentials()

def create_fba_model(solver = GLPK):

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
                                 atp_id=essentials['atp'],
                                 adp_id=essentials['adp'],
                                 pi_id=essentials['pi'],
                                 h2o_id=essentials['h2o'],
                                 h_id=essentials['h'],
                                 )

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

    ecoli.add_ribosome(rib,free_ratio=0.2)
    # http://bionumbers.hms.harvard.edu/bionumber.aspx?id=102348&ver=1&trm=rna%20polymerase%20half%20life&org=
    # Name          Fraction of active RNA Polymerase
    # Bionumber ID  102348
    # Value 	    0.17-0.3 unitless
    # Source        Bremer, H., Dennis, P. P. (1996) Modulation of chemical composition and other parameters of the cell by growth rate.
    #               Neidhardt, et al. eds. Escherichia coli and Salmonella typhimurium: Cellular
    #                       and Molecular Biology, 2nd ed. chapter 97 Table 1

    ecoli.add_rnap(rnap, free_ratio=0.75)

    ecoli.build_expression()
    ecoli.add_enzymatic_coupling(coupling_dict)

    if has_neidhardt:

        nt_ratios, aa_ratios = get_ratios()
        chromosome_len, gc_ratio = get_ecoli_gen_stats()
        kdeg_mrna, mrna_length_avg  = get_mrna_metrics()
        kdeg_enz,  peptide_length_avg   = get_enz_metrics()
        neidhardt_mu, neidhardt_rrel, neidhardt_prel, neidhardt_drel = get_neidhardt_data()

        ecoli.add_interpolation_variables()
        ecoli.add_dummies(nt_ratios=nt_ratios,
                          mrna_kdeg=kdeg_mrna,
                          mrna_length=mrna_length_avg,
                          aa_ratios=aa_ratios,
                          enzyme_kdeg=kdeg_enz,
                          peptide_length=peptide_length_avg)
        ecoli.add_protein_mass_requirement(neidhardt_mu, neidhardt_prel)
        ecoli.add_rna_mass_requirement(neidhardt_mu, neidhardt_rrel)
        ecoli.add_dna_mass_requirement(mu_values=neidhardt_mu,
                                       dna_rel=neidhardt_drel,
                                       gc_ratio=gc_ratio,
                                       chromosome_len=chromosome_len,
                                       dna_dict=dna_nucleotides)

    # Need to put after, because dummy has to be taken into account if used.
    ecoli.populate_expression()
    ecoli.add_trna_mass_balances()


    ecoli.print_info()

    need_relax = False

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


if __name__ == '__main__':
    model = create_etfl_model(False, False)