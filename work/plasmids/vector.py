from etfl.data.ecoli import kdeg_enz, kdeg_mrna, get_average_kcat, get_rnap, \
    get_nt_sequences
from etfl.core.genes import ExpressedGene
from etfl.core.rna import mRNA
from etfl.core.enzyme import Enzyme
from etfl.core.vector import TransModel, Plasmid

from utils import read_seq

from cobra.core import Metabolite, Reaction

def get_debug_plasmid(model, has_rnap=False):

    h2o_c = model.metabolites.h2o_c
    h2o_e = model.metabolites.h2o_e

    # CO2 Exchange
    h2o_tpp = Reaction(id='h2otpp', name='H2O Exchange')
    h2o_tpp.add_metabolites({
        h2o_c:  1,
        h2o_e: -1,
    })

    reaction_list = [h2o_tpp]

    #############
    #   Genes   #
    #############

    DBG = ExpressedGene(id='DBG_gene',
                        name='Debug synthase',
                        sequence='AATTTCGGTTGA'.lower())

    #############
    #  Enzymes  #
    #############

    DBG_enz = Enzyme(
        id='DBG',
        kcat=65 * 3600,
        kdeg=kdeg_enz,
        composition={'DBG_gene': 3}
    )

    ############
    # Coupling #
    ############

    coupling_dict = {'h2otpp': [DBG_enz]}

    ###########
    # Plasmid #
    ###########

    if 1:
        rnap_genes = []
        rnap = None

    gene_list = [DBG] + rnap_genes

    plasmid_seq = DBG.sequence

    my_plasmid = Plasmid(id_='debug',
                         sequence=plasmid_seq,
                         genes=gene_list,
                         reactions=reaction_list,
                         coupling_dict=coupling_dict,
                         rnap=rnap)
    my_plasmid.build_default_mrna(kdeg_mrna)

    return my_plasmid

def get_bdo_plasmid(model, has_rnap=False):

    # Plasmid pZS*-13S-ald-adh
    # US Patent 20,140,371,417.
    # Used in Andreozzi, Stefano, et al.
    # "Identification of metabolic engineering targets for the enhancement of
    # 1,4-butanediol production in recombinant E. coli using large-scale
    # kinetic models."
    # Metabolic engineering 35 (2016): 148-159.

    # (too complicated)

    # 2-gene BDO plasmid
    # Reshamwala, Shamlan MS, Shalini S. Deb, and Arvind M. Lali.
    # "A shortened, two-enzyme pathway for 2, 3-butanediol production in
    # Escherichia coli."
    # Journal of industrial microbiology & biotechnology 44.9 (2017): 1273-1277.

    ###############
    # Metabolites #
    ###############
    acetoin     = Metabolite(id='acetoin_c', name = 'Acetoin, cytosol',
                             formula = 'C4H8O2', compartment='c')
    diacetyl    = Metabolite(id='diacetyl_c', name = 'Diacetyl, cytosol',
                             formula = 'C4H6O2', compartment='c')
    bdo         = Metabolite(id='bdo_c', name='Meso-2,3-butanediol, cytosol',
                             formula = 'C4H10O2', compartment='c')
    bdo_e       = Metabolite(id='bdo_e', name='Meso-2,3-butanediol, extracellular',
                             formula = 'C4H10O2', compartment='e')
    acetolactate= model.metabolites.alac__S_c
    nadh        = model.metabolites.nadh_c
    nad         = model.metabolites.nad_c
    pyr         = model.metabolites.pyr_c
    co2         = model.metabolites.co2_c
    h           = model.metabolites.h_c
    #############
    # Reactions #
    #############

    # ALS: 2 pyr + h => alac + co2
    als = Reaction(id = 'ALS', name = 'Acetolactate synthase [Plasmid](Enterobacter)')
    als.add_metabolites({
        pyr:-2,
        h:-1,
        acetolactate:1,
        co2:1
    })

    # SALD: alac => diacetyl + co2 (spontaneous)
    sald = Reaction(id = 'SALD', name = 'Spontaneous acetolactate decarboxylation [Plasmid]')
    sald.add_metabolites({
        acetolactate:-1,
        diacetyl:1,
        co2:1
    })

    # AR: diacetyl + NADH => acetoin + NAD+#############

    ar = Reaction(id = 'AR', name = 'Acetoin oxidoreductase (NADH) [Plasmid]')
    ar.add_metabolites({
        diacetyl:-1,
        nadh:-1,
        h:-1,
        acetoin:1,
        nad:1
    })

    # BDH: acetoin + NADH => bdo + NAD+
    bdh = Reaction(id = 'BDH', name = 'Butanediol dehydrogenase [Plasmid]')
    bdh.add_metabolites({
        acetoin:-1,
        nadh:-1,
        h:-1,
        bdo:1,
        nad:1
    })

    # BDO Transport
    bdo_tp = Reaction(id='BDOtp', name = 'BDO Transport (extracellular)')
    bdo_tp.add_metabolites({
        bdo:-1,
        bdo_e:1
    })

    # BDO Exchange
    bdo_ex = Reaction(id='EX_BDO', name = 'BDO Exchange')
    bdo_ex.add_metabolites({
        bdo_e:-1,
    })

    reaction_list = [als, sald, ar, bdh, bdo_tp, bdo_ex]

    #############
    #   Genes   #
    #############

    ALS = ExpressedGene(id = 'EBA_als',
                        name = 'acetolactate synthase',
                        sequence = read_seq('data/ALS.txt'))
    AR = ExpressedGene (id = 'EBA_ar',
                        name = 'acetoin reductase',
                        sequence = read_seq('data/AR.txt'))
    fwd_als= 'GGAATTCCATATGAACAGTGAGAAACAGTC'
    rev_als= 'CCGAGCTCTTACAAAATCTGGCTGAGAT'
    fwd_ar = 'GGAATTCCATATGCAAAAAGTTGCTCTCGTAAC'
    rev_ar = 'CCGAGCTCTTAGTTGAACACCATCCCAC'

    pET = read_seq('data/pET.txt')

    #############
    #  Enzymes  #
    #############

    # ALS is 60 kDa, and is a homotrimer
    # Kaushal, Anju, Sunil Pabbi, and Prince Sharma.
    # "Characterization of 2, 3-butanediol-forming and valine-sensitive
    # α-acetolactate synthase of Enterobacter cloacae."
    # World Journal of Microbiology and Biotechnology 19.5 (2003): 487-493.

    # ALS kcat is 26.9 +- 2.5[s-1]
    # Choi, Myeong-Hyeon, et al.
    # "Molecular cloning and expression of Enterobacter aerogenes α-acetolactate
    # decarboxylase in pyruvate decarboxylase-deficient Saccharomyces cerevisiae
    # for efficient 2, 3-butanediol production."
    # Process biochemistry 51.2 (2016): 170-176.

    ALS_enz = Enzyme(
        id = 'ALS',
        kcat=26.9*3600,
        kdeg=kdeg_enz,
        composition = {'EBA_als':3}
    )

    # The third enzyme in the 2,3-butanediol pathway, AR, purified from A.aerogenes
    # consists of four equal-size sub-unitsof 25kDa(40).

    # Blomqvist, Kristina, et al.
    # "Characterization of the genes of the 2, 3-butanediol operons from Klebsiella
    # terrigena and Enterobacter aerogenes."
    # Journal of bacteriology 175.5 (1993): 1392-1404.
    # &
    # Stormer, Fredrik C.
    # "[106] 2, 3-Butanediol biosynthetic system in Aerobacter aerogenes."
    # Methods in enzymology. Vol. 41. Academic Press, 1975. 518-533.

    AR_enz = Enzyme(
        id = 'AR',
        kcat=get_average_kcat(),
        kdeg=kdeg_enz,
        composition = {'EBA_ar':4}
    )

    # Probably the catalytic activities are different
    BDH_enz = Enzyme(
        id = 'BDH',
        kcat=get_average_kcat(),
        kdeg=kdeg_enz,
        composition = {'EBA_ar':4}
    )

    ############
    # Coupling #
    ############

    coupling_dict = {'ALS':[ALS_enz], 'AR':[AR_enz], 'BDH':[BDH_enz]}

    # AR and BDH are catalyzed by the same enzyme

    ###########
    # Plasmid #
    ###########

    if has_rnap:
        rnap = get_rnap()
        rnap.id = 'plasmid_rnap'

        nt_seq_ecoli = get_nt_sequences()
        rnap_genes = list()
        for g in rnap.composition:
            this_gene = ExpressedGene(id='PLASMID_'+g,
                                      name='{}, Plasmid',
                                      sequence = nt_seq_ecoli[g])
            rnap_genes.append(this_gene)
    else:
        rnap_genes = []
        rnap = None

    gene_list = [ALS,AR] + rnap_genes

    plasmid_seq = pET \
                  + fwd_als + ALS.sequence + rev_als \
                  + fwd_ar + AR.sequence + rev_ar \
                  + ''.join(str(g.sequence) for g in rnap_genes)

    my_plasmid = Plasmid(id_ = 'pET-AR-ALS',
                         sequence = plasmid_seq,
                         genes = gene_list,
                         reactions = reaction_list,
                         coupling_dict = coupling_dict,
                         rnap = rnap)
    my_plasmid.build_default_mrna(kdeg_mrna)

    return my_plasmid
