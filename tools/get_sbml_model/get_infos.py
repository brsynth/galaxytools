from argparse import ArgumentParser
from libsbml import readSBMLFromFile
from taxonid import get_taxonid
from requests import get as r_get


def get_biomass_rxn(sbml_doc):
    """
    Returns the biomass reaction of the model

    Parameters
    ----------
    sbml_doc: libsbml.SBMLDocument
        SBML model

    Returns
    -------
    biomass_rxn: libsbml.Reaction
        Biomass reaction
    """
    reactions = sbml_doc.getModel().getListOfReactions()
    # Search for 'biomass' keyword in reaction name
    for rxn in reactions:
        if "biomass" in rxn.getName().lower():
            return rxn
    # Search for 'biomass' keyword in products
    # AND not in reactants
    for rxn in reactions:
        in_reactants = False
        for reac in rxn.getListOfReactants():
            if "biomass" in reac.getSpecies().lower():
                in_reactants = True
                break
        if not in_reactants:
            for prod in rxn.getListOfProducts():
                if "biomass" in prod.getSpecies().lower():
                    return rxn
    return None


def args():
    parser = ArgumentParser("Returns cell informations")
    parser.add_argument("infile", type=str, help="SBML input file (xml)")
    parser.add_argument("--biomassid", type=str, help="ID of biomass reaction")
    parser.add_argument("--taxonid", type=str, help="Taxonomy ID")
    parser.add_argument("--standalone", action="store_true", help="Standalone mode, e.g. do not retrieve taxonomy ID on Internet (true if --taxonid is provided)")
    parser.add_argument("--compartments-outfile", type=str, help="Path to store cell compartments")
    parser.add_argument("--biomassid-outfile", type=str, help="Path to store biomass reaction ID")
    parser.add_argument("--taxonid-outfile", type=str, help="Path to store host taxonomy ID")
    params = parser.parse_args()
    return params


def get_organism_from_bigg_model(model_id):
    """Try to retrieve organism info from BiGG Models for a given model ID."""
    url = f"http://bigg.ucsd.edu/api/v2/models/{model_id}"
    try:
        response = r_get(url)
        if response.status_code == 200:
            data = response.json()
            organism = data.get("organism")
            return organism
    except Exception as e:
        print(f"Error querying BiGG: {e}")
    return None

def get_taxon_id(input_name):
    """Try BiGG model name first, then NCBI directly."""
    print(f"Trying input: {input_name}")
    
    # Try resolving as a BiGG model
    organism = get_organism_from_bigg_model(input_name)
    if organism:
        print(f"Model '{input_name}' maps to organism: {organism}")
        taxon_id = get_taxonid(organism)
        if taxon_id:
            return taxon_id

    # If not a model, try directly as an organism name
    print(f"Trying NCBI search with input: {input_name}")
    return get_taxonid(input_name)


def entry_point():
    params = args()

    # test if the file exists
    with open(params.infile):
        pass

    sbml_doc = readSBMLFromFile(params.infile)

    compartments = sbml_doc.getModel().getListOfCompartments()
    comp_str = ""
    for comp in compartments:
        comp_str += f"{comp.getId()}\t{comp.getName()}\n"
    print("Compartments:")
    for comp in compartments:
        print(f"{comp.getId()}\t{comp.getName()}".replace("\n", " | "))
    if params.compartments_outfile:
        with open(params.compartments_outfile, "w") as f:
            f.write("#ID\tNAME\n")
            f.write(comp_str)

    if params.biomassid:
        biomass_rxn = sbml_doc.getModel().getReaction(params.biomassid)
    else:
        biomass_rxn = get_biomass_rxn(sbml_doc)
    if not biomass_rxn:
        print("Warning: unable to retrieve biomass reaction")
        biomass_id = ""
    else:
        biomass_id = biomass_rxn.getId()
    print(f"Biomass reaction ID: {biomass_id}")
    if params.biomassid_outfile:
        with open(params.biomassid_outfile, "w") as f:
            f.write("#ID\n")
            f.write(f"{biomass_id}\n")

    if params.taxonid:
        taxid = params.taxonid
    elif params.standalone:
        taxid = -1
    else:
        model_id = sbml_doc.getModel().getId()
        if model_id:
            taxid = get_taxon_id(sbml_doc.getModel().getId())
        if taxid == -1:
            # Try with model name
            model_name = sbml_doc.getModel().getName()
            if model_name:
                taxid = get_taxon_id(sbml_doc.getModel().getName())
    print(f"Taxonomy ID: {taxid}")

    if params.taxonid_outfile:
        with open(params.taxonid_outfile, "w") as f:
            f.write("#ID\n")
            f.write(f"{taxid}\n")


if __name__ == "__main__":
    entry_point()
