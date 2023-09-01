from argparse import ArgumentParser
from libsbml import (
    readSBMLFromFile
)
from requests import get as r_get


def entry_point():
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
        '--hostid',
        type=str,
        help='Extended name of the host organism'
    )
    parser.add_argument(
        '--taxid',
        type=str,
        help='Path to store host taxonomy ID'
    )
    params = parser.parse_args()

    sbml_doc = readSBMLFromFile(params.infile)

    if params.comp:
        compartments = sbml_doc.getModel().getListOfCompartments()
        with open(params.comp, 'w') as f:
            f.write('#ID\tNAME\n')
            for comp in compartments:
                f.write(f'{comp.getId()}\t{comp.getName()}\n')

    if params.biomass:
        reactions = sbml_doc.getModel().getListOfReactions()
        with open(params.biomass, 'w') as f:
            f.write('#ID\n')
            for rxn in reactions:
                if 'biomass' in rxn.getId().lower():
                    f.write(f'{rxn.getId()}\n')

    if params.taxid:
        hostname = ''

        # Model from BiGG
        if params.bigg:
            # Extended Name
            server = 'http://bigg.ucsd.edu/api/v2/models/'
            ext = params.hostid
            r = r_get(server+ext, headers={ "Content-Type" : "application/json"})
            if not r.ok:
                print(f"Warning: unable to retrieve host name for id {params.hostid}")
            else:
                try:
                    hostname = r.json()["organism"]
                except KeyError:
                    print(f"Warning: unable to retrieve host name for id {params.hostid}")
            if not hostname:
                taxid = ''
            else:
                # TAXON ID
                server = 'https://rest.ensembl.org'
                ext = f'/taxonomy/id/{hostname}?'
                r = r_get(server+ext, headers={ "Content-Type" : "application/json"})
                if not r.ok:
                    print(f"Warning: unable to retrieve taxonomy ID for host organism {hostname}")
                else:
                    try:
                        taxid = r.json()["id"]
                    except KeyError:
                        print(f"Warning: unable to retrieve taxonomy ID for host organism {hostname}")
                        taxid = ''

        # Model from user
        else:
            taxid = params.hostid

        with open(params.taxid, 'w') as f:
            f.write('#ID\n')
            f.write(f'{taxid}\n')


if __name__ == "__main__":
    entry_point()