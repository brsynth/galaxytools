from argparse import ArgumentParser
from libsbml import (
    readSBMLFromFile
)
#from requests import get as r_get
from taxonid import get_taxonid


def get_biomass_rxn(sbml_doc):
    '''
    Returns the biomass reaction of the model
    
    Parameters
    ----------
    sbml_doc: libsbml.SBMLDocument
        SBML model
        
    Returns
    -------
    biomass_rxn: libsbml.Reaction
        Biomass reaction
    '''
    reactions = sbml_doc.getModel().getListOfReactions()
    # Search for 'biomass' keyword in reaction name
    for rxn in reactions:
        if 'biomass' in rxn.getName().lower():
            return rxn
    # Search for 'biomass' keyword in products
    # AND not in reactants
    for rxn in reactions:
        in_reactants = False
        for reac in rxn.getListOfReactants():
            if 'biomass' in reac.getSpecies().lower():
                in_reactants = True
                break
        if not in_reactants:
            for prod in rxn.getListOfProducts():
                if 'biomass' in prod.getSpecies().lower():
                    return rxn
    return None


# def get_taxon_id(hostname):
#     '''
#     Returns the taxonomy ID of the host organism
    
#     Parameters
#     ----------
#     hostname: str
#         Extended name of the host organism
        
#     Returns
#     -------
#     taxid: str
#         Taxonomy ID of the host organism
#     '''
#     taxid = get_taxon_id(hostid)
#     hostname = ''
#     # Extended Name
#     server = 'http://bigg.ucsd.edu/api/v2/models/'
#     ext = hostid
#     r = r_get(server+ext, headers={ "Content-Type" : "application/json"})
#     if not r.ok:
#         print(f"Warning: unable to retrieve host name for id {hostid}")
#     else:
#         try:
#             hostname = r.json()["organism"]
#         except KeyError:
#             print(f"Warning: unable to retrieve host name for id {hostid}")
#     if not hostname:
#         taxid = ''
#     else:
#         # TAXON ID
#         server = 'https://rest.ensembl.org'
#         ext = f'/taxonomy/id/{hostname}?'
#         r = r_get(server+ext, headers={ "Content-Type" : "application/json"})
#         if not r.ok:
#             print(f"Warning: unable to retrieve taxonomy ID for host organism {hostname}")
#         else:
#             try:
#                 taxid = r.json()["id"]
#             except KeyError:
#                 print(f"Warning: unable to retrieve taxonomy ID for host organism {hostname}")
#                 taxid = ''
#     return taxid


def args():
    parser = ArgumentParser('Returns cell informations')
    parser.add_argument(
        'infile',
        type=str,
        help='SBML input file (xml)'
    )
    # argument to tag file from BiGG
    parser.add_argument(
        '--bigg',
        action='store_true',
        help='Tag file from BiGG'
    )
    parser.add_argument(
        '--comp',
        type=str,
        help='Path to store cell compartments'
    )
    parser.add_argument(
        '--biomass',
        type=str,
        help='Path to store biomass reaction ID'
    )
    parser.add_argument(
        '--biomass-id',
        type=str,
        help='ID of biomass reaction'
    )
    parser.add_argument(
        '--host-ext-name',
        type=str,
        help='Extended name of the host organism'
    )
    parser.add_argument(
        '--taxid',
        type=str,
        help='Path to store host taxonomy ID'
    )
    params = parser.parse_args()
    return params


def entry_point():

    params = args()

    sbml_doc = readSBMLFromFile(params.infile)

    compartments = sbml_doc.getModel().getListOfCompartments()
    comp_str = ''
    for comp in compartments:
        comp_str += f'{comp.getId()}\t{comp.getName()}\n'
    if params.comp:
        with open(params.comp, 'w') as f:
            f.write('#ID\tNAME\n')
            f.write(comp_str)
    else:
        print('Compartments:')
        for comp in compartments:
            print(f'{comp.getId()}\t{comp.getName()}'.replace('\n', ' | '))

    if params.biomass_id:
        biomass_rxn = sbml_doc.getModel().getReaction(params.biomass_id)
    else:
        biomass_rxn = get_biomass_rxn(sbml_doc)
    if not biomass_rxn:
        print('Warning: unable to retrieve biomass reaction')
        biomass_id = ''
    else:
        biomass_id = biomass_rxn.getId()
    if params.biomass:
        with open(params.biomass, 'w') as f:
            f.write('#ID\n')
            f.write(f'{biomass_id}\n')
    else:
        print(f'Biomass reaction ID: {biomass_id}')

    taxid = get_taxonid(params.host_ext_name)
    # # Model from BiGG
    # if params.bigg:
    #     taxid = get_taxon_id(params.host_ext_name)
    # # Model from user
    # else:
    #     taxid = params.host_ext_name

    if params.taxid:
        with open(params.taxid, 'w') as f:
            f.write('#ID\n')
            f.write(f'{taxid}\n')
    else:
        print(f'Taxonomy ID: {taxid}')


if __name__ == "__main__":
    entry_point()