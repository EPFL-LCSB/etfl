"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Utilities to create a small model from a 1 reaction model in FBA

"""

import cobra
from cobra.test import create_test_model

from pytfa.io import        load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from ..core.memodel import MEModel
from ..core.thermomemodel import ThermoMEModel

from ..data.ecoli import   get_model, get_thermo_data, get_coupling_dict, \
                        get_mrna_dict, get_rib, get_rnap, get_monomers_dict, \
                        get_nt_sequences, get_ratios, get_neidhardt_data, \
                        get_mrna_metrics, get_enz_metrics, \
                        remove_from_biomass_equation, get_ecoli_gen_stats

from optlang.exceptions import SolverError

import logging


CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'

solver = GUROBI
aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides = get_monomers_dict()


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
            ex_rxn.lower_bound = -0.1


def create_etfl_model(has_thermo, has_neidhardt,
                      prot_scaling=1e3,
                      mrna_scaling=1e3,
                      n_mu_bins = 64,
                      mu_max = 3,
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

    coupling_dict = get_coupling_dict(vanilla_model)


    # Initialize the model

    name = 'small_model_T{:1}E{:1}N{:1}_{}_enz_{}_bins.json'.format(
        has_thermo,
        True,
        has_neidhardt,
        len(coupling_dict),
        n_mu_bins)


    if has_thermo:

        thermo_data, lexicon, compartment_data = get_thermo_data()

        ecoli = ThermoMEModel(thermo_data,model = vanilla_model,
                              growth_reaction = growth_reaction_id,
                              mu_range = mu_range,
                              n_mu_bins = n_mu_bins,
                              max_enzyme_concentration = 1000,
                              prot_scaling = prot_scaling,
                              mrna_scaling = mrna_scaling,
                              name = name,
                              )
    else:
        ecoli = MEModel(model = vanilla_model,
                        growth_reaction = growth_reaction_id,
                        mu_range = mu_range,
                        n_mu_bins = n_mu_bins,
                        max_enzyme_concentration = 1000,
                        prot_scaling = prot_scaling,
                        mrna_scaling = mrna_scaling,
                        name = name,
                        )

    ecoli.name = name
    ecoli.logger.setLevel(logging.WARNING)


    # Solver settings
    ecoli.solver = solver
    ecoli.solver.configuration.verbosity = 1
    ecoli.solver.configuration.tolerances.feasibility = 1e-9
    if solver == 'optlang-gurobi':
        ecoli.solver.problem.Params.NumericFocus = 3
    ecoli.solver.configuration.presolve = True



    if has_thermo:
        # Annotate the cobra_model
        annotate_from_lexicon(ecoli, lexicon)
        apply_compartment_data(ecoli, compartment_data)

        # TFA conversion
        ecoli.prepare()
        ecoli.convert()#add_displacement = True)


    nt_sequences = get_nt_sequences()
    mrna_dict = get_mrna_dict(ecoli)
    rnap, rnap_genes = get_rnap()
    rib, rrna_genes, rprot_genes = get_rib()

    prune_to_genes = lambda the_dict:{k:v for k,v in the_dict.items() \
                    if k in vanilla_model.genes or \
                     k in rrna_genes or k in rprot_genes.values or k in rnap_genes}

    nt_sequences = prune_to_genes(nt_sequences)
    mrna_dict = prune_to_genes(mrna_dict)

    # Remove nucleotides and amino acids from biomass reaction as they will be
    # taken into account by the expression

    remove_from_biomass_equation(model = ecoli,
                                 nt_dict = rna_nucleotides,
                                 aa_dict = aa_dict,
                                 atp_id='atp_c',
                                 adp_id='adp_c',
                                 pi_id='pi_c',
                                 h2o_id='h2o_c',
                                 h_id='h_c',
                                 )

    ##########################
    ##    MODEL CREATION    ##
    ##########################

    ecoli.add_rnap(rnap)
    ecoli.add_ribosome(rib, free_ratio=0.2)
    ecoli.add_nucleotide_sequences(nt_sequences)
    ecoli.add_mrnas(mrna_dict.values())

    ecoli.build_expression(aa_dict = aa_dict,
                           rna_nucleotides= rna_nucleotides,
                           atp='atp_c',
                           amp='amp_c',
                           gtp='gtp_c',
                           gdp='gdp_c',
                           pi='pi_c',
                           ppi='ppi_c',
                           h2o='h2o_c',
                           h='h_c',
                           rnap_genes = rnap_genes,
                           rrna_genes = rrna_genes,
                           rprot_genes= rprot_genes
                           )
    ecoli.add_enzymatic_coupling(coupling_dict)
    ecoli.populate_expression()
    ecoli.add_degradation(rna_nucleotides_mp = rna_nucleotides_mp,
                          h2o='h2o_c',
                          h='h_c')

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
                          peptide_length=peptide_length_avg,
                          gtp='gtp_c',
                          gdp='gdp_c',
                          h2o='h2o_c',
                          h='h_c',
                          ppi='ppi_c')
        ecoli.add_protein_mass_requirement(neidhardt_mu, neidhardt_prel)
        ecoli.add_rna_mass_requirement(neidhardt_mu, neidhardt_rrel)
        ecoli.add_dna_mass_requirement(mu_values=neidhardt_mu,
                                       dna_rel=neidhardt_drel,
                                       gc_ratio=gc_ratio,
                                       chromosome_len=chromosome_len,
                                       dna_dict=dna_nucleotides)
    # Need to put after, because dummy has to be taken into account if used.
    ecoli.add_trna_mass_balances()


    ecoli.print_info()

    need_relax = False

    ecoli.repair()

    try:
        ecoli.optimize()

        print('Growth               : {}'.format(ecoli.solution.f))
        print(' - Ribosomes produced: {}'.format(ecoli.solution.x_dict.EZ_rib))
        print(' - RNAP produced: {}'.format(ecoli.solution.x_dict.EZ_rnap))

    except (AttributeError, SolverError):
        pass

    return ecoli


if __name__ == '__main__':
    model = create_etfl_model(False, False)