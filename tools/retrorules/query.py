import argparse
import json
import logging
import sys
from typing import Dict, List

import requests

BASE_URL = "https://retrorules.org/api/v0.7"


def from_ec_number(ec_number: str, min_diameter: int = None) -> List:
    url = f"{BASE_URL}/ecnumber"
    params = [("input", ec_number)]
    if min_diameter:
        params.append(("minDiameter", str(min_diameter)))
    return url, params


def from_substrate(substrate: str, min_diameter: int = None) -> List:
    url = f"{BASE_URL}/substrate"
    params = [("input", substrate)]
    if min_diameter:
        params.append(("minDiameter", str(min_diameter)))
    return url, params


def from_reaction_id(
    reaction_id: str, repository: str, min_diameter: int = None
) -> List:
    url = f"{BASE_URL}/reactionid"
    params = [("input", reaction_id)]
    params.append(("repo", repository))
    if min_diameter:
        params.append(("minDiameter", str(min_diameter)))
    return url, params


def from_inchi(inchi: str, min_diameter: int = None) -> List:
    url = f"{BASE_URL}/inchi"
    params = [("input", inchi)]
    if min_diameter:
        params.append(("minDiameter", str(min_diameter)))
    return url, params


def from_repository(repository: str, min_diameter: int = None) -> List:
    url = f"{BASE_URL}/repo"
    params = [("input", repository)]
    if min_diameter:
        params.append(("minDiameter", str(min_diameter)))
    return url, params


def from_smarts_id(smarts_ids: List[str], min_diameter: int = None) -> List:
    url = f"{BASE_URL}/smartsid"
    params = [("input", smart_id) for smart_id in smarts_ids]
    if min_diameter:
        params.append(("minDiameter", str(min_diameter)))
    return url, params


def query(url: str, params: Dict) -> Dict:
    response = requests.get(url, params=params)
    response.raise_for_status()
    return response.json()


def write_json(path: str, data: Dict):
    with open(path, "w") as fd:
        json.dump(data, fd, indent=4)


def main():
    parser = argparse.ArgumentParser(
        description="Query RetroRules API via command-line endpoints."
    )
    subparsers = parser.add_subparsers(dest="command")

    # Subcommand: EC number
    parser_ecn = subparsers.add_parser("ec-number", help="From EC number")
    parser_ecn.add_argument(
        "--input-ec-number-str", required=True, help="EC number such as 1.1.1.1"
    )
    parser_ecn.add_argument("--input-min-diameter-int", type=int, help="Min diameter")
    parser_ecn.add_argument(
        "--output-data-json", required=True, help="Output results, JSON format"
    )

    # Subcommand: Substrate
    parser_sub = subparsers.add_parser("substrate", help="From substrate")
    parser_sub.add_argument(
        "--input-substrate-str", required=True, help="Substrate label"
    )
    parser_sub.add_argument("--input-min-diameter-int", type=int, help="Min diameter")
    parser_sub.add_argument(
        "--output-data-json", required=True, help="Output results, JSON format"
    )

    # Subcommand: Reaction ID
    parser_rea = subparsers.add_parser("reaction-id", help="From Reaction ID")
    parser_rea.add_argument(
        "--input-reaction-id-str", required=True, help="Reaction ID"
    )
    parser_rea.add_argument("--input-repository-str", required=True, help="Repository")
    parser_rea.add_argument("--input-min-diameter-int", type=int, help="Min diameter")
    parser_rea.add_argument(
        "--output-data-json", required=True, help="Output results, JSON format"
    )

    # Subcommand: InChI
    parser_inc = subparsers.add_parser("inchi", help="From InChI")
    parser_inc.add_argument("--input-inchi-str", required=True, help="InChI")
    parser_inc.add_argument("--input-min-diameter-int", type=int, help="Min diameter")
    parser_inc.add_argument(
        "--output-data-json", required=True, help="Output results, JSON format"
    )

    # Subcommand: Repository
    parser_rep = subparsers.add_parser("repository", help="From Repository")
    parser_rep.add_argument("--input-repository-str", required=True, help="InChI")
    parser_rep.add_argument("--input-min-diameter-int", type=int, help="Min diameter")
    parser_rep.add_argument(
        "--output-data-json", required=True, help="Output results, JSON format"
    )

    # Subcommand: Smarts ID
    parser_sma = subparsers.add_parser("smarts-id", help="From Smarts ID")
    parser_sma.add_argument(
        "--input-smarts-id-str", nargs="+", required=True, help="Smarts ID"
    )
    parser_sma.add_argument("--input-min-diameter-int", type=int, help="Min diameter")
    parser_sma.add_argument(
        "--output-data-json", required=True, help="Output results, JSON format"
    )

    logging.info("Query RetroRules - start")
    args = parser.parse_args()

    try:
        logging.info("Build arguments")
        url, params = "", {}
        if args.command == "ec-number":
            url, params = from_ec_number(
                ec_number=args.input_ec_number_str,
                min_diameter=args.input_min_diameter_int,
            )
        elif args.command == "substrate":
            url, params = from_substrate(
                substrate=args.input_substrate_str,
                min_diameter=args.input_min_diameter_int,
            )
        elif args.command == "reaction-id":
            url, params = from_reaction_id(
                reaction_id=args.input_reaction_id_str,
                repository=args.input_repository_str,
                min_diameter=args.input_min_diameter_int,
            )
        elif args.command == "inchi":
            url, params = from_inchi(
                inchi=args.input_inchi_str,
                min_diameter=args.input_min_diameter_int,
            )
        elif args.command == "repository":
            url, params = from_repository(
                repository=args.input_repository_str,
                min_diameter=args.input_min_diameter_int,
            )
        elif args.command == "smarts-id":
            url, params = from_smarts_id(
                smarts_ids=args.input_smarts_id_str,
                min_diameter=args.input_min_diameter_int,
            )
        else:
            parser.print_help()
            sys.exit(1)

        logging.info("Query API")
        data = query(url=url, params=params)

        logging.info("Write data")
        write_json(path=args.output_data_json, data=data)
    except requests.HTTPError as e:
        logging.error(f"HTTP error: {e.response.status_code} - {e.response.text}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

    logging.info("Query RetroRules - end")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%d-%m-%Y %H:%M",
    )
    main()
