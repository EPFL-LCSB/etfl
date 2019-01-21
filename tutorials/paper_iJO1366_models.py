

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

from etfl.core import Enzyme, Ribosome, RNAPolymerase, ThermoMEModel, MEModel

from etfl.io.json import save_json_model
from pytfa.utils.logger import get_timestr

from etfl.data.ecoli import   get_model, get_thermo_data, get_coupling_dict, \
                        get_mrna_dict, get_rib, get_rnap, get_monomers_dict, \
                        get_nt_sequences, get_ratios, get_neidhardt_data, \
                        get_mrna_metrics, get_enz_metrics, \
                        remove_from_biomass_equation, get_ecoli_gen_stats, \
                        get_lloyd_coupling_dict, \
                        get_essentials

from etfl.optim.config import standard_solver_config

from optlang.exceptions import SolverError

from multiprocessing import Pool


data_dir = '../organism_data/info_ecoli'

# solver = 'optlang-cplex'
solver = 'optlang-gurobi'

# McCloskey2014 values
glc_uptake = 7.54
glc_uptake_std = 0.56
observed_growth_std = 0.02
observed_growth = 0.61

growth_reaction_id = 'BIOMASS_Ec_iJO1366_WT_53p95M'


def create_model(has_thermo, has_expression, has_neidhardt, n_mu_bins = 128):
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
    mu_range = [0, 3.5]
    n_mu_bins = n_mu_bins

    time_str = get_timestr()

    coupling_dict = get_coupling_dict(vanilla_model,
                                      mode = 'kmax',
                                      # mode = 'kcat',
                                      atps_name = 'ATPS4rpp',
                                      infer_missing_enz=True)
    # coupling_dict = get_lloyd_coupling_dict(vanilla_model)
    aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides = get_monomers_dict()
    essentials = get_essentials()

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

    ecoli.solver = solver
    standard_solver_config(ecoli)

    if has_thermo:
        # Annotate the cobra_model
        annotate_from_lexicon(ecoli, lexicon)
        apply_compartment_data(ecoli, compartment_data)

        # TFA conversion
        ecoli.prepare()

        # ecoli.reactions.GLUDy.thermo['computed'] = False
        # ecoli.reactions.DHAtpp.thermo['computed'] = False
        # ecoli.reactions.MLTP2.thermo['computed'] = False
        # ecoli.reactions.G3PD2.thermo['computed'] = False
        ecoli.reactions.MECDPS.thermo['computed'] = False
        ecoli.reactions.NDPK4.thermo['computed'] = False
        ecoli.reactions.TMDPP.thermo['computed'] = False
        ecoli.reactions.ARGAGMt7pp.thermo['computed'] = False

        ecoli.convert()#add_displacement = True)


    mrna_dict = get_mrna_dict(ecoli)
    nt_sequences = get_nt_sequences()
    rnap = get_rnap()
    rib = get_rib()

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
    ecoli.growth_reaction.lower_bound = observed_growth - 3*observed_growth_std

    need_relax = False

    ecoli.repair()

    try:
        ecoli.optimize()
    except (AttributeError, SolverError):
        need_relax = True

    # from ipdb import set_trace; set_trace()

    if has_thermo and need_relax:
        final_model, slack_model, relax_table = relax_dgo(ecoli)
        # final_model, slack_model, relax_table = relax_dgo(ecoli, in_place = True)
    else:
        final_model = ecoli


    final_model.growth_reaction.lower_bound = 0
    solution = final_model.optimize()
    print('Objective            : {}'.format(final_model.solution.objective_value))
    print(' - Glucose uptake    : {}'.format(final_model.reactions.EX_glc__D_e.flux))
    print(' - Growth            : {}'.format(final_model.growth_reaction.flux))
    print(' - Ribosomes produced: {}'.format(final_model.ribosome.X))
    print(' - RNAP produced: {}'.format(final_model.rnap.X))
    try:
        print(' - DNA produced: {}'.format(final_model.solution.raw.DN_DNA))
    except AttributeError:
        pass

    filepath = 'models/{}'.format(final_model.name)
    # save_json_model(final_model, filepath)

    final_model.logger.info('Build complete for model {}'.format(final_model.name))

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
        (False, True, True),
        # (True, True, False),
        # (True, True, True),
        # (False, True, False),
        ]

    for mc in model_calls:
        models[mc] = create_model(*mc)

    #
    # pool = Pool()
    #
    # for mc in model_calls:
    #     def this_callback(result, mc=mc):
    #         models[mc] = result
    #     pool.apply_async(create_model, args=mc, callback=this_callback)
    #
    # pool.close()
    # pool.join()

    print('Completed')


    # Make thermo model
    # make_thermo_model()

    # Save FBA model
    # make_fba_model()