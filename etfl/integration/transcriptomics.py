from ..optim.variables import ModelVariable
from ..optim.constraints import GeneConstraint

class RelativeTranscriptomicsLB(GeneConstraint):
    """
    Represents a lower bound on mRNA ratio in relative transcriptomics
    """

    prefix = 'RTL_'


class RelativeTranscriptomicsUB(GeneConstraint):
    """
    Represents an upper bound on mRNA ratio in relative transcriptomics
    """

    prefix = 'RTU_'

class ReferenceLevel(ModelVariable):
    """
    Represents the reference level for relative transcriptomics
    """

    prefix = 'RL_'


def integrate_relative_transcriptomics(model, lower_bounds, upper_bounds, base=2):
    """

    Integrates log-ratio expression data to mRNA levels in ETFL

    :param model: an ETFL model
    :type model: etfl.core.memodel.MEModel
    :param lower_bounds:
    :type lower_bounds: dict or pandas.Series
    :param upper_bounds:
    :type upper_bounds: dict or pandas.Series
    :return:
    """

    rl = model.add_variable(kind=ReferenceLevel,
                            hook = model,
                            id_ = 'relative_transcriptomics',
                            lb = 0,
                            ub = 100,
                            queue=False)

    for the_mrna in model.mrnas:
        # Let gr be a silent reference gene
        # logb(gi/gr) = K
        # <=> gi = gr*b^K
        # if L <= K <= U:
        # gr*b^L <= gi <= gr*b^U
        #     { gi - gr*b^U <= 0 (expr_u)
        # <=> { and
        #     { gr*b^L - gi <= 0 (expr_l)

        try:
            log_ub_ratio = upper_bounds[the_mrna.id]
            log_lb_ratio = lower_bounds[the_mrna.id]
        except (KeyError, IndexError) as e:
            model.logger.warn('No expression data found for mRNA {}'.format(the_mrna.id))
            continue

        expr_u =      the_mrna.concentration - rl*(base**log_ub_ratio)
        expr_l = -1 * the_mrna.concentration + rl*(base**log_lb_ratio)

        model.add_constraint(kind=RelativeTranscriptomicsUB,
                             hook=the_mrna.gene,
                             expr=expr_u,
                             ub=0,
                             queue=True)

        model.add_constraint(kind=RelativeTranscriptomicsLB,
                             hook=the_mrna.gene,
                             expr=expr_l,
                             ub=0,
                             queue=True)
    # Push the model queue
    # And regenerate constraints and variable properties
    model.repair()
    return True