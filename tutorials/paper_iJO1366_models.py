

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#       E. Coli. Debug Model
#5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Values reported are roughly estimated for debugging the model. These are not
# actual reviewed values.

import cobra

import logging


from pytfa.io import        load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.optim.relaxation import relax_dgo

from therme.core import Enzyme, Ribosome, RNAPolymerase, ThermoMEModel, MEModel

from therme.io.json import save_json_model
from pytfa.utils.logger import get_timestr

from etfl_data import   get_model, get_thermo_data, get_coupling_dict, \
                        get_mrna_dict, get_rib, get_rnap, get_monomers_dict, \
                        get_nt_sequences, get_ratios, get_neidhardt_data, \
                        get_mrna_metrics, get_enz_metrics, remove_from_biomass_equation

from pytfa.thermo import ThermoModel

data_dir = '../organism_data/info_ecoli'

# solver = 'optlang-cplex'
solver = 'optlang-gurobi'

# McCloskey2014 values
glc_uptake = 7.54
glc_uptake_std = 0.56
observed_growth = 0.61 - 0.02

growth_reaction_id = 'BIOMASS_Ec_iJO1366_WT_53p95M'


def create_model(has_thermo, has_expression, has_neidhardt):
    #------------------------------------------------------------
    # Initialisation
    #------------------------------------------------------------

    assert has_expression == True

    # this hack works because we are using the solver switch to update the var
    # names in the solver but really we should not do this
    # TODO: clean up model.sanitize_varnames
    vanilla_model = get_model('optlang_glpk')
    vanilla_model.reactions.EX_glc__D_e.lower_bound = -1 * glc_uptake - glc_uptake_std
    vanilla_model.reactions.EX_glc__D_e.upper_bound = -1 * glc_uptake + glc_uptake_std

    vanilla_model.objective = growth_reaction_id
    fba_sol = vanilla_model.slim_optimize()
    mu_0 = fba_sol
    mu_range = [0, 4]
    n_mu_bins = 256

    time_str = get_timestr()

    coupling_dict = get_coupling_dict(vanilla_model)
    aa_dict, rna_nucleotides = get_monomers_dict()

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
                        scaling = 1e3,
                        name = name,
                        )
    else:
        ecoli = MEModel(model = vanilla_model,
                        growth_reaction = growth_reaction_id,
                        mu_range = mu_range,
                        n_mu_bins = n_mu_bins,
                        max_enzyme_concentration = 1000,
                        scaling = 1e3,
                        name = name,
                        )

    ecoli.name = name
    ecoli.logger.setLevel(logging.WARNING)

    # Set the solver
    ecoli.solver = solver
    ecoli.solver.configuration.verbosity = 1
    ecoli.solver.configuration.tolerances.feasibility = 1e-9
    if solver == 'optlang_gurobi':
        ecoli.solver.problem.Params.NumericFocus = 3
    ecoli.solver.configuration.presolve = True

    if has_thermo:
        # Annotate the cobra_model
        annotate_from_lexicon(ecoli, lexicon)
        apply_compartment_data(ecoli, compartment_data)

        # TFA conversion
        ecoli.prepare()
        ecoli.convert()#add_displacement = True)


    mrna_dict = get_mrna_dict(ecoli)
    nt_sequences = get_nt_sequences()
    rib, rrna_genes, rprot_genes = get_rib()
    rnap, rnap_genes = get_rnap()

    # Remove nucleotides and amino acids from biomass reaction as they will be
    # taken into account by the expression

    remove_from_biomass_equation(model = ecoli,
                                 nt_dict = rna_nucleotides,
                                 aa_dict = aa_dict)

    ##########################
    ##    MODEL CREATION    ##
    ##########################

    ecoli.add_rnap(rnap)
    ecoli.add_ribosome(rib, free_ratio=0.2)
    ecoli.add_nucleotide_sequences(nt_sequences)
    ecoli.add_mrnas(mrna_dict.values())

    ecoli.build_expression( aa_dict = aa_dict,
                            nt_dict = rna_nucleotides,
                            atp='atp_c',
                            amp='amp_c',
                            gtp='gtp_c',
                            gdp='gdp_c',
                            ppi='ppi_c',
                            h2o='h2o_c',
                            h='h_c',
                            rnap_genes = rnap_genes,
                            rrna_genes = rrna_genes,
                            rprot_genes= rprot_genes
                            )
    ecoli.add_enzymatic_coupling(coupling_dict)
    ecoli.populate_expression()
    ecoli.add_degradation()

    if has_neidhardt:

        nt_ratios, aa_ratios = get_ratios()
        kdeg_mrna, mrna_length_avg  = get_mrna_metrics()
        kdeg_enz,  peptide_length_avg   = get_enz_metrics()
        neidhardt_mu, neidhardt_rrel, neidhardt_prel = get_neidhardt_data()

        ecoli.add_interpolation_variables()
        ecoli.add_dummies(nt_ratios=nt_ratios,
                          mrna_kdeg=kdeg_mrna,
                          mrna_length=mrna_length_avg,
                          aa_ratios=aa_ratios,
                          enzyme_kdeg=kdeg_enz,
                          peptide_length=peptide_length_avg)
        ecoli.add_protein_mass_requirement(neidhardt_mu, neidhardt_prel)
        ecoli.add_rna_mass_requirement(neidhardt_mu, neidhardt_rrel)


    ecoli.print_info()
    ecoli.growth_reaction.lower_bound = observed_growth

    try:
        ecoli.optimize()
    except AttributeError:
        need_relax = True

    if has_thermo and need_relax:
        final_model, slack_model, relax_table = relax_dgo(ecoli)
    else:
        final_model = ecoli


    final_model.growth_reaction.lower_bound = 0
    solution = final_model.optimize()
    print('Growth               : {}'.format(final_model.solution.f))
    print(' - Ribosomes produced: {}'.format(final_model.solution.x_dict.EZ_rib))
    print(' - RNAP produced: {}'.format(final_model.solution.x_dict.EZ_rnap))

    filepath = 'models/{}'.format(ecoli.name)
    save_json_model(final_model, filepath)


def make_fba_model():
    ecoli = get_model('optlang_glpk')
    ecoli.reactions.EX_glc__D_e.lower_bound = -1 * glc_uptake - glc_uptake_std
    ecoli.reactions.EX_glc__D_e.upper_bound = -1 * glc_uptake + glc_uptake_std

    # ecoli.objective = growth_reaction_id
    ecoli.optimize()

    from cobra.io.json import save_json_model

    save_json_model(ecoli, 'models/iJO1366_T0E0N0_{}'.format(get_timestr()))


def make_thermo_model():
    vanilla_model = get_model('optlang_glpk')

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

    # Models defined by Thermo - Expression - Neidhardt
    model_calls = [
                        (False,  True,   False),
                        (True,   True,   False),
                        (False,  True,   True),
                        (True,   True,   True),
                        ]

    for mc in model_calls:
        create_model(*mc)

    # Make thermo model
    make_thermo_model()

    # Save FBA model
    make_fba_model()