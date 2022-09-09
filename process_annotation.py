from models import ChainInfo, UnitCorrespondence, UnitPairsInteractions, UnitInfo, UnitPairsDistances
from database import db_session
import csv


def remove_ife(ife_list, ife_to_remove):
    """
    This removes the ife from the original list once it has been assigned.

    :param ife_list: list
        A list of ribosome RNA component ifes (tRNA/s, mRNA etc)
    :param ife_to_remove: string
        The ife to remove
    :return: list
        A list of ribosome RNA component ifes, one ife shorter
    """
    try:
        ife_list.remove(ife_to_remove)
    except Exception as e:
        print e
        print 'This ife was not found in the list'

    return ife_list


def get_components_ife(ife):
    """
    This query gets all the RNA chains from the structure of interest except the query
    chain (16S in this case) and builds the corresponding ifes.

    :param ife: string
        Query structure ife
    :return: list
        A list of RNA component ifes in the structure of interest
    """
    
    pdb, _, chain = ife.split('|')
    molecule_type = 'Polyribonucleotide (RNA)'
    chain_ids = []
    
    # Get all RNA chains from the ribosome except the query chain (16S in this example)

    with db_session() as session:
        query = session.query(ChainInfo) \
                       .filter(ChainInfo.pdb_id == pdb) \
                       .filter(ChainInfo.entity_macromolecule_type == molecule_type) \
                       .filter(ChainInfo.chain_name != chain)

        for row in query:
            chain_ids.append(row.chain_name)

    # Build the ifes from the pdb and chain information. Assume model num is 1
    ife_list = ['{}|1|{}'.format(pdb, elem) for elem in chain_ids]
    query_info = (pdb, chain)

    return query_info, ife_list


def get_units_correspondence(units_list, test_info):
    """
    This query returns a list of corresponding unit_ids from a particular chain in a structure based on a
    given list of unit_ids

    #TODO
    Consider symmetry operators

    :param units_list: list
        Query list of unit_ids
    :param test_info: tuple
        A tuple of (pdb, chain)
    :return: units_corr: list
        A list of corresponding unit_ids in the test_info
    """
    pdb, chain = test_info[0], test_info[1]
    units_corr = []

    with db_session() as session:
        for unit in units_list:
            query = session.query(UnitCorrespondence) \
                           .filter(UnitCorrespondence.unit_id_1 == unit) \
                           .filter(UnitCorrespondence.pdb_id_2 == pdb) \
                           .filter(UnitCorrespondence.chain_name_2 == chain)

            for row in query:
                units_corr.append(row.unit_id_2)

    return units_corr


def units_correspondence(reference_units, test_ssu, test_lsu):
    """
    This creates a dict of correspondence in the query structure based on the ref dict. The keys of this dict
    would be the relevant chain-chain pair (ssu-lsu, ssu-mRNA, ssu-atRNA etc) while the values would be the
    list of corresponding units in the query structure. For example, the key 'ssu_lsu' would have a list of
    nucleotides in ssu that should interact with the lsu in the query structure as its value.

    :param reference_units: dict
        A dictionary that contains the reference nucleotides that make specific contacts between ifes
    :param test_ssu: tuple
        A tuple of (pdb, chain) of ssu
    :param test_lsu: tuple
        A tuple of (pdb, chain) of lsu
    :return: dict
        A dict of correspondence in the query structure based on the ref dict
        
    """
    corr_dict = {}
    for key, value in reference_units.items():
        label_1, _, label_3 = key.split('_')
        new_key = '{}_{}'.format(label_1, label_3)
        if label_1 == 'ssu':
            units_list = get_units_correspondence(value, test_ssu)
            corr_dict[new_key] = units_list
        else:
            units_list = get_units_correspondence(value, test_lsu)
            corr_dict[new_key] = units_list

    return corr_dict


def infer_atrna_ife(ifes, ribosome_components):
    """
    Check for the presence of mRNA and assign the A-tRNA ife. This needs to be done because nucleotides (A1492 &
    A1493 Ec) in the decoding loop can interact with mRNA as well.

    :param ifes: list
        A list of possible ifes that might interact with the SSU A-tRNA binding site nts
    :param ribosome_components: dict
        A dict containing ribosome RNA component molecules as keys and the assigned ifes as values
    :return: String
        The assigned A-tRNA ife
    """
    if len(ifes) == 1:
        atrna_ife = list(ifes)[0]
    # assume there could only be two 
    elif len(ifes) > 1:
        ifes.remove(ribosome_components['mRNA'])
        atrna_ife = list(ifes)[0]
    else:
        atrna_ife = None

    return atrna_ife


def infer_etrna_ife(possible_etrna_ife, ribosome_components):
    """
    Check if the E-site trna is the same as the P-site trna. We are checking for the presence of E-tRNA by "looking" at
    chains that could possibly interact with the E-tRNA binding site in the LSU. Since P-tRNA can form P/E state, we
    need to test for this possibility.

    :param possible_etrna_ife: string
        Possible E-tRNA ife
    :param ribosome_components: dict
        A dict containing ribosome RNA component molecules as keys and the assigned ifes as values
    :return: String
        The assigned E-tRNA ife
    """
    if possible_etrna_ife == ribosome_components['peptidyl-trna']:
        etrna_ife = None
    else:
        etrna_ife = possible_etrna_ife

    return etrna_ife


def get_interacting_ife(corr_units, components_ife):
    """
    This query returns the ife that makes pairwise interaction/s with the corresponding nts of interest.

    :param corr_units: list
        A list of unit_ids
    :param components_ife: list
        A list of ribosome RNA components ife
    :return: string
        The interacting ife
    
    interacting_ife = None
    for ife in components_ife:
        # '5J7L|1|AA|%'
        chain2 = '{}|%'.format(ife)
        for unit in corr_units:
            with db_session() as session:
                query = session.query(UnitPairsInteractions) \
                               .filter(UnitPairsInteractions.unit_id_1 == unit) \
                               .filter(UnitPairsInteractions.unit_id_2.like(chain2)) \
                               .count()

                if query >= 1:
                    interacting_ife = ife

            return interacting_ife
    """
    with db_session() as session:
        for ife in components_ife:
            # '5J7L|1|AA|%'
            chain2 = '{}|%'.format(ife)
            for unit in corr_units:
                query = session.query(UnitPairsInteractions) \
                               .filter(UnitPairsInteractions.unit_id_1 == unit) \
                               .filter(UnitPairsInteractions.unit_id_2.like(chain2)) \
                               .count()

                if query:
                    return ife


def get_interacting_ifes(corr_units, components_ife):
    """
    This query returns a list of ifes that makes pairwise interaction/s with the corresponding nts of interest.

    :param corr_units: list
        A list of unit_ids
    :param components_ife: list
        A list of ribosome RNA components ife
    :return: set
        A set of interacting ifes
    """
    with db_session() as session:
        nr_ifes = set()
        for ife in components_ife:
            # '5J7L|1|AA|%'
            unit2 = '{}|%'.format(ife)
            for unit in corr_units:
                query = session.query(UnitPairsInteractions) \
                               .filter(UnitPairsInteractions.unit_id_1 == unit) \
                               .filter(UnitPairsInteractions.unit_id_2.like(unit2))

                for row in query:
                    chain2 = '|'.join(row.unit_id_2.split('|')[:3])
                    nr_ifes.add(chain2)

        return nr_ifes


def infer_interacting_ife(corr_units, component_ifes, ribosome_components=None, chain_type=None):
    """
    This assigns the ife that interact with the corresponding nts of interest.

    :param corr_units: list
        A list of unit_ids
    :param component_ifes: list
        A list of ribosome RNA components ife
    :param ribosome_components: dict
        A dict containing ribosome RNA component molecules as keys and the assigned ifes as values
    :param chain_type: string
        Molecule type. The value is None by default. Do extra processing if it's an A- or a P-tRNA
    :return: string
        The assigned ife
    """
    ife = None

    if chain_type is None:
        ife = get_interacting_ife(corr_units, component_ifes)

    elif chain_type == 'atrna':
        possible_ifes = get_interacting_ifes(corr_units, component_ifes)
        ife = infer_atrna_ife(possible_ifes, ribosome_components)

    elif chain_type == 'etrna':
        possible_ife = get_interacting_ife(corr_units, component_ifes)
        ife = infer_etrna_ife(possible_ife, ribosome_components)

    # remove ife from the original list once it has been assigned
    component_ifes = remove_ife(component_ifes, ife)

    return ife, component_ifes


def check_neighboring_contacts(ife, units_corr):
    """
    This query checks whether the ife makes any interactions with the nucleotides of interest.

    :param ife: string
        The ife of interest
    :param units_corr: list
        A list of unit_ids
    :return: bool
        True if the ife makes interaction/s with the nts of interest
    """

    unit2 = '{}|%'.format(ife)
    contacts_made = False

    with db_session() as session:
        for unit1 in units_corr:
            query = session.query(UnitPairsInteractions) \
                           .filter(UnitPairsInteractions.unit_id_1 == unit1) \
                           .filter(UnitPairsInteractions.unit_id_2.like(unit2)) \
                           .count()

            if query >= 1:
                contacts_made = True

        return contacts_made


def get_trna_acceptor_nts(ife_info):
    """
    This query returns the last three nucleotides (CCA- end, usually position 74-76) of a tRNA molecule
    :param ife_info: tuple
        A tuple of (pdb, chain)
    :return: list
        A list of units_ids
    """
    with db_session() as session:
        acceptor_units = []
        query = session.query(UnitInfo) \
                       .filter(UnitInfo.pdb_id == ife_info[0]) \
                       .filter(UnitInfo.chain.like(ife_info[1])) \
                       .order_by(UnitInfo.chain_index.desc()) \
                       .limit(3)
        
        for row in query:
            acceptor_units.append(row.unit_id)

        return acceptor_units


def get_protein_contacts_ife(ife):
    """
    This reads the RNA-protein contacts file generated by FR3D and returns the protein factor ife that interacts with
    the CCA-end of tRNA in the LSU (using a distance cutoff of 4.5 Angstroms)
    :param ife: string
        The tRNA ife
    :return: string
        The protein factor ife
    """

    pdb, _, chain = ife.split("|")
    ife_info = (pdb, chain)
    acceptor_units = get_trna_acceptor_nts(ife_info)

    chains = set()
    with db_session() as session:
        query = session.query(UnitPairsDistances) \
                       .join(UnitInfo, UnitPairsDistances.unit_id_2 == UnitInfo.unit_id) \
                       .filter(UnitInfo.unit_type_id == 'aa') \
                       .filter(UnitPairsDistances.distance <= 8.5) \
					   .filter(UnitPairsDistances.unit_id_1.in_(acceptor_units))

        for row in query:
            protein_chain = '|'.join(row.unit_id_2.split('|')[:3])
            chains.add(protein_chain)

    """
    nr_ifes = set()
    with open('/Applications/mamp/htdocs/contact_list_rename/contact_list_{}.csv'.format(pdb), 'U') as f:
        csv_reader = csv.reader(f, delimiter=',')
        for lines in csv_reader:
            for unit in acceptor_units:
                if unit == str(lines[0]):
                    protein_ife = '|'.join(lines[3].split('|')[:3])
                    nr_ifes.add(protein_ife)
    # Assume we get only one protein ife
    #return list(nr_ifes)[0]
    """
    if chains:
        return list(chains)[0]



def get_compound_name(ife):
    """
    This query returns the compound name of the ife
    :param ife: string
        The query ife
    :return: string
        compound name
    """
    pdb, _, chain = ife.split('|')
    compound_name = None
    
    with db_session() as session:
        query = session.query(ChainInfo) \
                       .filter(ChainInfo.pdb_id == pdb) \
                       .filter(ChainInfo.chain_name == chain)

        # Should be just one record. Do we need to loop?
        for row in query:
            compound_name = row.compound

        return compound_name


def infer_tRNA_state(ribosome_components, trna_type, corr_dict):
    """
    This infers the state of the tRNAs based on the contacts they make both on SSU and LSU. Here we are assuming the
    tRNA to be a full-length molecule. This is not true since some of them can just be anticodon stem-loops (ASLs) or
    CCA-fragments in the lsu. As such, we need to consider the length as well.

    :param ribosome_components: dict
        A dict containing ribosome RNA component molecules as keys and the assigned ifes as values
    :param trna_type: string
        This could be aminoacyl, peptidyl or exit tRNA
    :param corr_dict: dict
        A dict of correspondence in the query structure based on the ref dict
    :return: string
        The state of tRNA (A/A, A/P ap/AP etc)
    """

    complete_state = None
    if trna_type == 'aminoacyl-trna':

        if ribosome_components[trna_type] is not None:

            LSU_contact = check_neighboring_contacts(ribosome_components[trna_type], corr_dict['lsu_atrna'])
            SSU_neighboring_contact = check_neighboring_contacts(ribosome_components[trna_type], corr_dict['ssu_ptrna'])
            LSU_neighboring_contact = check_neighboring_contacts(ribosome_components[trna_type], corr_dict['lsu_ptrna'])

            # Chimeric state
            if SSU_neighboring_contact is True:
                SSU_state = 'ap'
            # Classic state
            else:
                SSU_state = 'A'

            # Chimeric state
            if LSU_contact is True and LSU_neighboring_contact is True:
                LSU_state = 'AP'
            # Classic state
            elif LSU_contact is True and LSU_neighboring_contact is False:
                LSU_state = 'A'
            # Hybrid state
            elif LSU_contact is False and LSU_neighboring_contact is True:
                LSU_state = 'P'
            # Interacting with protein factor in LSU - A/T. A/R states
            elif LSU_contact is False and LSU_neighboring_contact is False:
                trna_length = get_chain_length(ribosome_components[trna_type])
                # ASL fragment
                if trna_length < 30:
                    LSU_state = '*'
                else:
                    protein_ife = get_protein_contacts_ife(ribosome_components[trna_type])
                    if protein_ife:
                        protein_factor_name = get_compound_name(protein_ife)
                        LSU_state = protein_factor_name
                    else:
                        LSU_state = '*'
            # Undefined
            else:
                LSU_state = '*'

            complete_state = '{}/{}'.format(SSU_state, LSU_state)

    elif trna_type == 'peptidyl-trna':

        if ribosome_components[trna_type] is not None:

            LSU_contact = check_neighboring_contacts(ribosome_components[trna_type], corr_dict['lsu_ptrna'])
            SSU_neighboring_contact = check_neighboring_contacts(ribosome_components[trna_type], corr_dict['ssu_etrna'])
            LSU_neighboring_contact = check_neighboring_contacts(ribosome_components[trna_type], corr_dict['lsu_etrna'])

            # Chimeric state
            if SSU_neighboring_contact is True:
                SSU_state = 'pe'
            # Classic state
            else:
                SSU_state = 'P'

            # Chimeric state
            if LSU_contact is True and LSU_neighboring_contact is True:
                LSU_state = 'PE'
            # Classic state
            elif LSU_contact is True and LSU_neighboring_contact is False:
                LSU_state = 'P'
            # Hybrid state
            elif LSU_contact is False and LSU_neighboring_contact is True:
                LSU_state = 'E'
            # Interacting with protein factor in LSU
            elif LSU_contact is False and LSU_neighboring_contact is False:
                trna_length = get_chain_length(ribosome_components[trna_type])
                # ASL fragment
                if trna_length < 30:
                    LSU_state = '*'
                else:
                    protein_ife = get_protein_contacts_ife(ribosome_components[trna_type])
                    if protein_ife:
                        protein_factor_name = get_compound_name(protein_ife)
                        LSU_state = protein_factor_name
                    else:
                        LSU_state = '*'
            # Undefined
            else:
                LSU_state = '*'

            complete_state = '{}/{}'.format(SSU_state, LSU_state)

    elif trna_type == 'exit-trna':
        # The E-tRNA can only be in the E/E state.
        if ribosome_components[trna_type] is not None:
            complete_state = 'E/E'

    return complete_state


def get_chain_length(chain_name):
    chain_pdb, _, chain_id = chain_name.split("|")

    with db_session() as session:
        query = session.query(ChainInfo) \
                       .filter(ChainInfo.pdb_id == chain_pdb) \
                       .filter(ChainInfo.chain_name == chain_id)

        for row in query:
            chain_length = row.chain_length
        
        return chain_length

"""
def get_units_distances(unit_id, components_dict):
    mrna_ife = '{}|%'.format(components_dict['mRNA'])
    ssu_ife = '{}|%'.format(components_dict['SSU_16S'])

    units_list = []
    query = UnitPairsDistances.query.filter(UnitPairsDistances.unit_id_1 == unit_id) \
        .filter(UnitPairsDistances.unit_id_2.notlike(mrna_ife)) \
        .filter(UnitPairsDistances.unit_id_2.notlike(ssu_ife))

    for row in query:
        unit2 = row.unit_id_2
        split_components = unit2.split('|')
        if len(split_components[3]) == 1:
            units_list.append('|'.join(split_components[:3]))

    return units_list
"""