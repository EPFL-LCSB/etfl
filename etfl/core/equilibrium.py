from .genes import ExpressedGene
from ..optim.variables import RNAPUsage, FreeEnzyme, BinaryActivator, \
    LinearizationVariable
from ..optim.constraints import GeneConstraint, EnzymeRatio,\
    LinearizationConstraint
from pytfa.optim.reformulation import petersen_linearization


class RNAPBindingEquilibrium(GeneConstraint):
    prefix = 'PBEq_'


def linearize_product(model, b, x, queue=False):
    """

    :param model:
    :param b: the binary variable
    :param x: the continuous variable
    :param queue: whether to queue the variables and constraints made
    :return:
    """

    # Linearization step for ga_i * [E]
    z_name = '__MUL__'.join([b.name, x.name])
    # Add the variables
    model_z_u = model.add_variable(kind=LinearizationVariable,
                                  hook=model,
                                  id_=z_name,
                                  lb=0,
                                  ub=x.ub,
                                  queue=False)

    big_m = x.ub*10

    z_u, new_constraints = petersen_linearization(b=b, x=x, M=big_m,
                                                  z=model_z_u)

    # Add the constraints:
    for cons in new_constraints:
        model.add_constraint(kind=LinearizationConstraint,
                             hook=model,
                             id_=cons.name,
                             expr=cons.expression,
                             # expr=new_expression,
                             ub=cons.ub,
                             lb=cons.lb,
                             queue=True)

    model._push_queue()
    return model_z_u

def add_rnap_binding_equilibrium_constraints(model, the_rnap, Kb):
    """

    :param model:
    :param the_rnap:
    :param Kb:
    :return:
    """

    # #1. Find if there is a ratio constraint on the RNAP - this is incompatible
    # all_ratio_cons = model.get_constraints_of_type(EnzymeRatio)
    #
    # for rnap_id in model.rnap:
    #     try:
    #         model.remove_constraint(all_ratio_cons.get_by_id(rnap_id))
    #     except KeyError:
    #         pass

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
    dna_hat_resolution = max([x-y for x,y in zip(sorted_dna_conc[1:], sorted_dna_conc[:-1])])


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

        # scaling
        # expr *= model.dna.scaling_factor

        #3.3 Add it to the model

        model.add_constraint(kind=RNAPBindingEquilibrium,
                             hook=the_gene,
                             expr=expr,
                             ub= copy_number*dna_hat_resolution/2,
                             lb=-copy_number*dna_hat_resolution/2,
                             queue=True)
    model._push_queue()
    model.regenerate_constraints()
    model.regenerate_variables()
