from argparse import ArgumentParser
from libsbml import (
    readSBMLFromFile
)


def entry_point():
    parser = ArgumentParser('Returns cell informations')
    parser.add_argument(
        'infile',
        type=str,
        help='SBML input file (xml)'
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

if __name__ == "__main__":
    entry_point()