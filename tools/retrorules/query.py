import argparse
import json
import logging
import sys
from typing import Dict, Tuple

import requests

BASE_URL = "https://retrorules.org/api"


def from_templates(
    smarts_str: str,
    template_ids_str: str,
    reaction_ids_str: str,
    datasets_str: str,
    chemical_domain_str: str,
    ec_number_str: str,
    min_radius_int: int,
    valid_str: str,
    dedup_str: str,
    limit_int: int,
    offset_int: int,
) -> Tuple:
    url = f"{BASE_URL}/templates"
    params = []
    if smarts_str:
        params.append(("q", smarts_str))
    if template_ids_str:
        params.append(("template_ids", ",".join(template_ids_str)))
    if reaction_ids_str:
        params.append(("reaction_ids", ",".join(reaction_ids_str)))
    if datasets_str and datasets_str != "any":
        params.append(("datasets", datasets_str))
    if chemical_domain_str and chemical_domain_str != "any":
        params.append(("chemical_domain", chemical_domain_str))
    if ec_number_str:
        params.append(("ec", ec_number_str))
    if min_radius_int is not None:
        params.append(("min_radius", str(min_radius_int)))
    if valid_str and valid_str != "any":
        params.append(("valid", valid_str))
    if dedup_str and dedup_str != "any":
        params.append(("dedup", dedup_str))
    if limit_int:
        params.append(("limit", str(limit_int)))
    if offset_int:
        params.append(("offset", str(offset_int)))
    return url, params


def from_templates_summary(template_id_str: str) -> Tuple:
    url = f"{BASE_URL}/templates/{template_id_str}/summary"
    params = {}
    return url, params


def from_templates_sources(template_id_str: str) -> Tuple:
    url = f"{BASE_URL}/templates/{template_id_str}/sources"
    params = {}
    return url, params


def from_templates_count(
    smarts_str: str,
    template_ids_str: str,
    reaction_ids_str: str,
    datasets_str: str,
    chemical_domain_str: str,
    ec_number_str: str,
    min_radius_int: int,
    valid_str: str,
    dedup_str: str,
) -> Tuple:
    url = f"{BASE_URL}/templates_count"
    params = []
    if smarts_str:
        params.append(("q", smarts_str))
    if template_ids_str:
        params.append(("template_ids", ",".join(template_ids_str)))
    if reaction_ids_str:
        params.append(("reaction_ids", ",".join(reaction_ids_str)))
    if datasets_str and datasets_str != "any":
        params.append(("datasets", datasets_str))
    if chemical_domain_str and chemical_domain_str != "any":
        params.append(("chemical_domain", chemical_domain_str))
    if ec_number_str:
        params.append(("ec", ec_number_str))
    if min_radius_int is not None:
        params.append(("min_radius", str(min_radius_int)))
    if valid_str and valid_str != "any":
        params.append(("valid", valid_str))
    if dedup_str and dedup_str != "any":
        params.append(("dedup", dedup_str))
    return url, params


def from_templates_export(
    generation_token_str: str,
    smarts_str: str,
    template_ids_str: str,
    reaction_ids_str: str,
    datasets_str: str,
    chemical_domain_str: str,
    ec_number_str: str,
    min_radius_int: int,
    valid_str: str,
    dedup_str: str,
) -> Tuple:
    url = f"{BASE_URL}/templates_export"
    params = []
    if generation_token_str:
        params.append(("gen_token", generation_token_str))
    if smarts_str:
        params.append(("q", smarts_str))
    if template_ids_str:
        params.append(("template_ids", ",".join(template_ids_str)))
    if reaction_ids_str:
        params.append(("reaction_ids", ",".join(reaction_ids_str)))
    if datasets_str and datasets_str != "any":
        params.append(("datasets", datasets_str))
    if chemical_domain_str and chemical_domain_str != "any":
        params.append(("chemical_domain", chemical_domain_str))
    if ec_number_str:
        params.append(("ec", ec_number_str))
    if min_radius_int is not None:
        params.append(("min_radius", str(min_radius_int)))
    if valid_str and valid_str != "any":
        params.append(("valid", valid_str))
    if dedup_str and dedup_str != "any":
        params.append(("dedup", dedup_str))
    return url, params


def query(url: str, params: Dict):
    response = requests.get(url, params=params)
    response.raise_for_status()
    return response


def write_json(path: str, data: Dict):
    with open(path, "w") as fd:
        json.dump(data, fd, indent=4)


def write_tab(path: str, data: str):
    with open(path, "w") as fd:
        fd.write(data)


def main():
    parser = argparse.ArgumentParser(
        description="Query RetroRules API via command-line endpoints."
    )
    subparsers = parser.add_subparsers(dest="command")

    # Subcommand: templates
    parser_tem = subparsers.add_parser("templates", help="From templates")
    parser_tem.add_argument("--input-smarts-str", help="Exact SMARTS")
    parser_tem.add_argument(
        "--input-template-ids-str",
        nargs="*",
        help="Space separated list of template IDs",
    )
    parser_tem.add_argument(
        "--input-reaction-ids-str",
        nargs="*",
        help="Space separated list of reaction IDs",
    )
    parser_tem.add_argument(
        "--input-datasets-str",
        default="any",
        choices=["any", "metanetx", "rhea", "uspto"],
        help="Select a specific database",
    )
    parser_tem.add_argument(
        "--input-chemical-domain-str",
        default="any",
        choices=["any", "biochem", "orgchem"],
    )
    parser_tem.add_argument(
        "--input-ec-number-str",
        help="EC number to filter templates",
    )
    parser_tem.add_argument(
        "--input-min-radius-int",
        type=int,
        help="Single radius of the template",
    )
    parser_tem.add_argument(
        "--input-valid-str",
        default="true",
        choices=["any", "true", "false"],
        help="By default only valid templates are returned",
    )
    parser_tem.add_argument(
        "--input-dedup-str",
        default="true",
        choices=["true", "false"],
        help="By default deduplicated templates are returned",
    )
    parser_tem.add_argument(
        "--input-limit-int",
        type=int,
        help="Limit number of returned templates",
    )
    parser_tem.add_argument(
        "--input-offset-int",
        type=int,
        help="Offset the list of returned templates",
    )
    parser_tem.add_argument(
        "--output-data-json",
        required=True,
        help="Path to output JSON file",
    )

    # Subcommand: templates-summary
    parser_tem_sum = subparsers.add_parser("templates-summary", help="From templates-summary")
    parser_tem_sum.add_argument("--input-template-id-str", required=True, help="Template ID")
    parser_tem_sum.add_argument(
        "--output-data-json",
        required=True,
        help="Path to output JSON file",
    )

    # Subcommand: templates-sources
    parser_tem_sou = subparsers.add_parser("templates-sources", help="From templates-sources")
    parser_tem_sou.add_argument("--input-template-id-str", required=True, help="Template ID")
    parser_tem_sou.add_argument(
        "--output-data-json",
        required=True,
        help="Path to output JSON file",
    )

    # Subcommand: templates-count
    parser_cou = subparsers.add_parser("templates-count", help="From templates-count")
    parser_cou.add_argument("--input-smarts-str", help="Exact SMARTS")
    parser_cou.add_argument(
        "--input-template-ids-str",
        nargs="*",
        help="Space separated list of template IDs",
    )
    parser_cou.add_argument(
        "--input-reaction-ids-str",
        nargs="*",
        help="Space separated list of reaction IDs",
    )
    parser_cou.add_argument(
        "--input-datasets-str",
        default="any",
        choices=["any", "metanetx", "rhea", "uspto"],
        help="Select a specific database",
    )
    parser_cou.add_argument(
        "--input-chemical-domain-str",
        default="any",
        choices=["any", "biochem", "orgchem"],
    )
    parser_cou.add_argument(
        "--input-ec-number-str",
        help="EC number to filter templates",
    )
    parser_cou.add_argument(
        "--input-min-radius-int",
        type=int,
        help="Single radius of the template",
    )
    parser_cou.add_argument(
        "--input-valid-str",
        default="true",
        choices=["any", "true", "false"],
        help="By default only valid templates are returned",
    )
    parser_cou.add_argument(
        "--input-dedup-str",
        default="true",
        choices=["true", "false"],
        help="By default deduplicated templates are returned",
    )
    parser_cou.add_argument(
        "--output-data-json",
        required=True,
        help="Path to output JSON file",
    )

    # Subcommand: templates-export
    parser_exp = subparsers.add_parser("templates-export", help="From templates-export")
    parser_exp.add_argument(
        "--input-generation-token-str",
        help="Generation token from RetroRules web interface",
    )
    parser_exp.add_argument("--input-smarts-str", help="Exact SMARTS")
    parser_exp.add_argument(
        "--input-template-ids-str",
        nargs="*",
        help="Space separated list of template IDs",
    )
    parser_exp.add_argument(
        "--input-reaction-ids-str",
        nargs="*",
        help="Space separated list of reaction IDs",
    )
    parser_exp.add_argument(
        "--input-datasets-str",
        default="any",
        choices=["any", "metanetx", "rhea", "uspto"],
        help="Select a specific database",
    )
    parser_exp.add_argument(
        "--input-chemical-domain-str",
        default="any",
        choices=["any", "biochem", "orgchem"],
    )
    parser_exp.add_argument(
        "--input-ec-number-str",
        help="EC number to filter templates",
    )
    parser_exp.add_argument(
        "--input-min-radius-int",
        type=int,
        help="Single radius of the template",
    )
    parser_exp.add_argument(
        "--input-valid-str",
        default="true",
        choices=["any", "true", "false"],
        help="By default only valid templates are returned",
    )
    parser_exp.add_argument(
        "--input-dedup-str",
        default="true",
        choices=["true", "false"],
        help="By default deduplicated templates are returned",
    )
    parser_exp.add_argument(
        "--output-data-json",
        help="Path to output JSON file",
    )
    parser_exp.add_argument(
        "--output-data-csv",
        help="Path to output CSV file",
    )
    parser_exp.add_argument(
        "--output-data-tsv",
        help="Path to output TSV file",
    )

    logging.info("Query RetroRules - start")
    args = parser.parse_args()

    try:
        logging.info("Build arguments")
        url, params = "", {}
        if args.command == "templates":
            url, params = from_templates(
                smarts_str=args.input_smarts_str,
                template_ids_str=args.input_template_ids_str,
                reaction_ids_str=args.input_reaction_ids_str,
                datasets_str=args.input_datasets_str,
                chemical_domain_str=args.input_chemical_domain_str,
                ec_number_str=args.input_ec_number_str,
                min_radius_int=args.input_min_radius_int,
                valid_str=args.input_valid_str,
                dedup_str=args.input_dedup_str,
                limit_int=args.input_limit_int,
                offset_int=args.input_offset_int,
            )
        elif args.command == "templates-summary":
            url, _ = from_templates_summary(
                template_id_str=args.input_template_id_str,
            )
        elif args.command == "templates-sources":
            url, _ = from_templates_sources(
                template_id_str=args.input_template_id_str,
            )
        elif args.command == "templates-count":
            url, params = from_templates_count(
                smarts_str=args.input_smarts_str,
                template_ids_str=args.input_template_ids_str,
                reaction_ids_str=args.input_reaction_ids_str,
                datasets_str=args.input_datasets_str,
                chemical_domain_str=args.input_chemical_domain_str,
                ec_number_str=args.input_ec_number_str,
                min_radius_int=args.input_min_radius_int,
                valid_str=args.input_valid_str,
                dedup_str=args.input_dedup_str,
            )
        elif args.command == "templates-export":
            url, params = from_templates_export(
                generation_token_str=args.input_generation_token_str,
                smarts_str=args.input_smarts_str,
                template_ids_str=args.input_template_ids_str,
                reaction_ids_str=args.input_reaction_ids_str,
                datasets_str=args.input_datasets_str,
                chemical_domain_str=args.input_chemical_domain_str,
                ec_number_str=args.input_ec_number_str,
                min_radius_int=args.input_min_radius_int,
                valid_str=args.input_valid_str,
                dedup_str=args.input_dedup_str,
            )
        else:
            parser.print_help()
            sys.exit(1)

        logging.info("Query API")
        response = query(url=url, params=params)

        logging.info("Write data")
        if "output_data_json" in vars(args) and args.output_data_json:
            write_json(path=args.output_data_json, data=response.json())
        if "output_data_csv" in vars(args) and args.output_data_csv:
            write_tab(path=args.output_data_csv, data=response.text)
        if "output_data_tsv" in vars(args) and args.output_data_tsv:
            write_tab(path=args.output_data_tsv, data=response.text)
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
