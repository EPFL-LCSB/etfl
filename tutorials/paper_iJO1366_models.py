

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#       E. Coli. Debug Model
#5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Values reported are roughly estimated for debugging the model. These are not
# actual reviewed values.


import logging


from pytfa.io import        load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.optim.relaxation import relax_dgo

from therme.core import Enzyme, Ribosome, RNAPolymerase, ThermoMEModel, MEModel

from therme.io.json import save_json_model
from pytfa.utils.logger import get_timestr

from therme.data.ecoli import   get_model, get_thermo_data, get_coupling_dict, \
                        get_mrna_dict, get_rib, get_rnap, get_monomers_dict, \
                        get_nt_sequences, get_ratios, get_neidhardt_data, \
                        get_mrna_metrics, get_enz_metrics, \
                        remove_from_biomass_equation, get_ecoli_gen_stats

from optlang.exceptions import SolverError


data_dir = '../organism_data/info_ecoli'

# solver = 'optlang-cplex'
solver = 'optlang-gurobi'

# McCloskey2014 values
glc_uptake = 7.54
glc_uptake_std = 0.56
observed_growth = 0.61 - 0.02

growth_reaction_id = 'BIOMASS_Ec_iJO1366_WT_53p95M'


def create_model(has_thermo, has_expression, has_neidhardt,
                 prot_scaling = 1e6, mrna_scaling=1e6):
    #------------------------------------------------------------
    # Initialisation
    #------------------------------------------------------------

    assert has_expression == True

    # this hack works because we are using the solver switch to update the var
    # names in the solver but really we should not do this
    # TODO: clean up model.sanitize_varnames
    vanilla_model = get_model('optlang-glpk')
    vanilla_model.reactions.EX_glc__D_e.lower_bound = -1 * glc_uptake - glc_uptake_std
    vanilla_model.reactions.EX_glc__D_e.upper_bound = -1 * glc_uptake + glc_uptake_std

    vanilla_model.objective = growth_reaction_id
    fba_sol = vanilla_model.slim_optimize()
    mu_0 = fba_sol
    mu_range = [0, 4]
    n_mu_bins = 256

    time_str = get_timestr()

    coupling_dict = get_coupling_dict(vanilla_model)
    aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides = get_monomers_dict()

    # Initialize the model

    name = 'iJO1366_T{:1}E{:1}N{:1}_{}_enz_{}_bins_{}.json'.format(
        has_thermo,
        has_expression,
        has_neidhardt,
        len(coupling_dict),
        n_mu_bins,
        time_str)

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

        ecoli.reactions.GLUDy.thermo['computed'] = False
        # ecoli.reactions.DHAtpp.thermo['computed'] = False
        ecoli.reactions.MLTP2.thermo['computed'] = False

        ecoli.convert()#add_displacement = True)


    mrna_dict = get_mrna_dict(ecoli)
    nt_sequences = get_nt_sequences()
    rib, rrna_genes, rprot_genes = get_rib()
    rnap, rnap_genes = get_rnap()

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

    # Need to be put after, because dummy has to be taken into account if used.
    ecoli.add_trna_mass_balances()


    ecoli.print_info()
    ecoli.growth_reaction.lower_bound = observed_growth

    need_relax = False

    ecoli.repair()

    try:
        ecoli.optimize()
    except (AttributeError, SolverError):
        need_relax = True

    # from ipdb import set_trace; set_trace()

    if has_thermo and need_relax:
        final_model, slack_model, relax_table = relax_dgo(ecoli)
    else:
        final_model = ecoli


    final_model.growth_reaction.lower_bound = 0
    solution = final_model.optimize()
    print('Objective            : {}'.format(final_model.solution.f))
    print(' - Growth            : {}'.format(final_model.solution.x_dict[growth_reaction_id]))
    print(' - Ribosomes produced: {}'.format(final_model.solution.x_dict.EZ_rib))
    print(' - RNAP produced: {}'.format(final_model.solution.x_dict.EZ_rnap))

    filepath = 'models/{}'.format(final_model.name)
    save_json_model(final_model, filepath)

    return final_model


def make_fba_model():
    ecoli = get_model(solver)
    ecoli.reactions.EX_glc__D_e.lower_bound = -1 * glc_uptake - glc_uptake_std
    ecoli.reactions.EX_glc__D_e.upper_bound = -1 * glc_uptake + glc_uptake_std

    # ecoli.objective = growth_reaction_id
    ecoli.optimize()

    from cobra.io.json import save_json_model

    save_json_model(ecoli, 'models/iJO1366_T0E0N0_{}.json'.format(get_timestr()))


def make_thermo_model():
    vanilla_model = get_model(solver)

    name = 'iJO1366_T1E0N0_{}'.format(get_timestr())

    vanilla_model.reactions.EX_glc__D_e.lower_bound = -1 * glc_uptake - glc_uptake_std
    vanilla_model.reactions.EX_glc__D_e.upper_bound = -1 * glc_uptake + glc_uptake_std

    vanilla_model.objective = growth_reaction_id
    vanilla_model.optimize()

    thermo_data, lexicon, compartment_data = get_thermo_data()

    ecoli = ThermoMEModel(thermo_data, model=vanilla_model,
                          name=name, )

    need_relax = False

    try:
        ecoli.optimize()
    except AttributeError:
        need_relax = True

    if need_relax:
        final_model, slack_model, relax_table = relax_dgo(ecoli)
    else:
        final_model = ecoli

    from pytfa.io.json import save_json_model

    save_json_model(final_model, 'models/'+name)

if __name__ == '__main__':

    models = dict()

    # Models defined by Thermo - Expression - Neidhardt
    model_calls = [
                        (False,  True,   False),
                        (False,  True,   True),
                        # (True,   True,   False),
                        # (True,   True,   True),
                        ]

    for mc in model_calls:
        models[mc] = create_model(*mc)

    # Make thermo model
    # make_thermo_model()

    # Save FBA model
    # make_fba_model()