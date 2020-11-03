from pytfa.optim.reformulation import linearize_product
from .genes import ExpressedGene
from ..optim.variables import RNAPUsage, FreeEnzyme, BinaryActivator
from pytfa.optim import LinearizationVariable, LinearizationConstraint
from ..optim.constraints import GeneConstraint, EnzymeRatio, \
    EnzymeConstraint, TotalCapacity, InterpolationConstraint
from ..core.allocation import interpolate_growth_data

from pytfa.optim.reformulation import petersen_linearization
from pytfa.optim.variables import ModelVariable
from pytfa.optim.utils import symbol_sum

from numpy import mean, median, diff


class RNAPBindingEquilibrium(GeneConstraint):
    prefix = 'PBEq_'

class SigmaBindingEquilibrium(EnzymeConstraint):
    prefix = 'SBEq_'

class HoloFreeRNAPRatio(ModelVariable):
    prefix = 'HFR_'


def add_rnap_binding_equilibrium_constraints(model, the_rnap, Kb):
    """

    :param model:
    :param the_rnap:
    :param Kb:
    :return:
    """

    #1. Find if there is a ratio constraint on the RNAP - this is incompatible
    remove_rnap_ratio_constraints(model)

    #2. Create the linearized bilinear product delta_u*RNAP_F

    RNAP_F = model.get_variables_of_type(FreeEnzyme).get_by_id(the_rnap.id).variable
    activation_vars = model.get_variables_of_type(BinaryActivator)

    zs = {delta_u.name: linearize_product(model, delta_u, RNAP_F, queue=True)
          for delta_u in activation_vars}

    # TODO: cleanup
    dna_cons = model.interpolation_constraint.DNA_interpolation.constraint
    dna_hat = {k.name:v for k,v in
               dna_cons.get_linear_coefficients(
                   activation_vars.list_attr('variable')).items()}

    free_rnap_times_gene_conc = sum(zs[delta_u]*dna_hat[delta_u] for delta_u in zs)

    sorted_dna_conc = sorted(dna_hat.values())
    dna_hat_resolution = median([x-y for x,y in zip(sorted_dna_conc[1:], sorted_dna_conc[:-1])])


    all_rnap_usage = model.get_variables_of_type(RNAPUsage)

    #3. For each gene:

    for the_gene in model.genes:
        if not isinstance(the_gene, ExpressedGene):
            continue

        copy_number = the_gene.copy_number

        #3.1. get rnap usage

        RNAP_l = all_rnap_usage.get_by_id(the_gene.id)

        #3.2. create constraint

        expr = Kb * RNAP_l - copy_number * free_rnap_times_gene_conc

        # scaling by ~1/1e-7
        sigma = the_rnap.scaling_factor
        expr /= sigma

        #3.3 Add it to the model

        model.add_constraint(kind=RNAPBindingEquilibrium,
                             hook=the_gene,
                             expr=expr,
                             ub= copy_number*dna_hat_resolution/(2*sigma),
                             lb=-copy_number*dna_hat_resolution/(2*sigma),
                             queue=True)
    model._push_queue()
    model.regenerate_constraints()
    model.regenerate_variables()


def remove_rnap_ratio_constraints(model):
    all_ratio_cons = model.get_constraints_of_type(EnzymeRatio)
    for rnap_id in model.rnap:
        try:
            cons = all_ratio_cons.get_by_id(rnap_id)
            # Make equality into inequality
            # 1.0*EF_rnap - 0.75*EZ_rnap == 0 --> >= 0
            # cons.constraint.ub = None
            model.remove_constraint(cons)
        except KeyError:
            pass


def add_sigma_factor(model, rnap_id, sigma_factor, holoenzyme,
                     rnap_activity_array, Kb, max_holo_to_free_ratio=10):
    """
    Transforms the model to take sigma factor activation of the RNAP into
    account

    :param rnap_id: the ID of the RNAP (should be included in the model)
    :param sigma_factor: {<rnap_id> : (sigma_factor_object, holoenzyme_object)}
    :param holoenzyme: {<rnap_id> : (sigma_factor_object, holoenzyme_object)}
    :return:
    """

    #0. If necessary, remove the RNAP ratio constraint
    remove_rnap_ratio_constraints(model)

    #1. Get the proper rnap
    rnap = model.rnap[rnap_id]
    free_rnap = model.free_enzyme.get_by_id(rnap_id)

    #2. Add holoenzyme and sigma factor to the model
    model.add_enzymes([sigma_factor, holoenzyme])
    model._push_queue()

    # index the new cons and variables
    model.regenerate_constraints()
    model.regenerate_variables()

    #3. Add rnap + sigma -> rnap reaction
    # rnap_activation = EnzymaticReaction(id='rnap_activation',
    #                                     name='RNAP_activation')
    # self.add_reaction(rnap_activation)
    # v_holo = rnap_activation.forward_variable - rnap_activation.reverse_variable
    # Remove metabolites from the existing complexation reaction
    holoenzyme.complexation.add_metabolites({k:-v for k,v in
                                             holoenzyme.complexation.metabolites.items()},
                                            rescale=False)
    v_holo = holoenzyme.complexation.net

    # Consumes one sigma
    sigma_mb = model.enzyme_mass_balance.get_by_id(sigma_factor.id)
    sigma_mb.change_expr(sigma_mb.expr - v_holo)

    # Consumes one RNAP
    rnap_mb = model.enzyme_mass_balance.get_by_id(rnap_id)
    rnap_mb.change_expr(rnap_mb.expr - v_holo)

    # Makes one holoenzyme
    # holo_mb = self.enzyme_mass_balance.get_by_id(holoenzyme.id)
    # holo_mb.change_expr(rnap_mb.expr + v_holo)

    model.recompute_transcription()
    model.recompute_translation()
    model.recompute_allocation()

    #4. Build the Activity ratio variable and discretize the input

    mu_values, q = rnap_activity_array
    q_vals, q_sum = \
       interpolate_growth_data(model, q, mu_values)
    # Quantify the resolution on the q_var approx
    q_epsilon = mean(diff(q_vals))

    q_var = model.add_variable(HoloFreeRNAPRatio,
                               id_ = rnap.id,
                               hook=model,
                               lb=0,
                               ub=max_holo_to_free_ratio,
                               queue=False)

    q_expr = q_var - q_sum
    model.add_constraint(kind=InterpolationConstraint,
                         hook=model,
                         id_='holo_vs_free_ratio',
                         expr=q_expr,
                         lb=-q_epsilon/2, #account for approximation resolution
                         ub= q_epsilon/2,
                         queue = True,
                         )

    #5. Build the sigma factor equilibrium constraint
    # Detail the calculations in function documentation
    # q = holoRNAP / RNAP_free
    # Equilibrium: [RNAP_free][sigma] = K_b * [holoRNAP]
    #                         [sigma] - K_b * q = 0
    equilibrium_expr =  sigma_factor.concentration - Kb * q_var
    model.add_constraint(kind = SigmaBindingEquilibrium,
                         hook = sigma_factor,
                         expr = equilibrium_expr/Kb, #for scaling
                         ub = 0, #q_var already accounts for the resolution error
                         lb = 0,
                         queue = True)
    model._push_queue()

    #6. Build the linearized variable
    # Since q = holoRNAP / RNAP_free
    # We need to linearize q * RNAP_free and write the constraint
    # holoRNAP =  q * RNAP_free

    # Build z =   λ_0 * q_0 * [E]
    #           + λ_1 * q_1 * [E]
    #           + ...
    #           + λ_N * q_N * [E]

    # They are sorted because it's a dictlist, and are created sequentially.
    activation_vars = model.get_variables_of_type(BinaryActivator)

    elements = list()
    for i, lambda_i in enumerate(activation_vars):
        # Linearization step for lambda_i * [E]
        model_z_i = linearize_product(model=model, b=lambda_i, x=free_rnap.variable, queue=True)
        elements.append(q_vals[i] * model_z_i)

    holo_free_expr = holoenzyme.concentration \
                     - free_rnap.scaling_factor * symbol_sum(elements)

    model.add_constraint(EnzymeRatio,
                         hook=holoenzyme,
                         expr=holo_free_expr/free_rnap.scaling_factor, #for scaling
                         ub=0,
                         lb=0,
                         queue=True
                         )

    #7. Replace RNAP variable by holo-RNAP in the total capacity constraint,
    #TODO: edit this up.
    # We want:
    # holoRNAP = sum RNAP_l     (1)


    all_total_capa = model.get_constraints_of_type(TotalCapacity)
    the_rnap_capa_cons = all_total_capa.get_by_id(rnap_id)

    subs_dict = {rnap.variable:holoenzyme.variable}
    the_holo_capa_expr = the_rnap_capa_cons.expr.subs(subs_dict)  #(1)

    # Remove previous rnap capa cons
    # Make new holoRNAP capa cons (1)

    model.remove_constraint(the_rnap_capa_cons)
    model.add_constraint(TotalCapacity,
                         id_ = holoenzyme.id,
                         hook=model,
                         expr=the_holo_capa_expr,
                         ub = 0,
                         lb = 0,
                         queue = True)

    # 8. Build total capa cons
    # RNAP_tot = holoRNAP + RNAP_free
    # Make new RNAP capa cons (2)

    sigma = min(x.scaling_factor for x in (holoenzyme, rnap, free_rnap))
    new_total_rnap_capa = holoenzyme.concentration \
                          + free_rnap.unscaled \
                          - rnap.concentration
    model.add_constraint(TotalCapacity,
                         id_ =rnap.id,
                         hook=model,
                         expr=new_total_rnap_capa/sigma,
                         ub = 0,
                         lb = 0,
                         queue = True)





    # self.rnap[holoenzyme.id] = holoenzyme
    # self.rnap.pop(rnap_id)

    model._push_queue()
    model.regenerate_constraints()
    model.regenerate_variables()
