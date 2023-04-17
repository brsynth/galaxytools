#!/usr/bin/env python
import argparse
import os
import yaml


def lint_shed_yaml(repository: str, owner: str) -> None:
    shed_yaml = os.path.join(repository, ".shed.yml")
    if not os.path.exists(shed_yaml):
        raise ValueError("No .shed.yml file found, skipping.")
    data = {}
    try:
        with open(shed_yaml) as fh:
            data = yaml.safe_load(fh)
    except Exception as e:
        raise ValueError("Failed to parse .shed.yml file [%s]" % unicodify(e))
    name = data.get("name")
    if name is None or len(name) < 3:
        raise ValueError("Name if not found or invalid")
    sowner = data.get("owner")
    if sowner is None or sowner != owner:
        raise ValueError("Owner is not found or invalid")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input-repository-str", required=True, help="Path repository"
    )
    parser.add_argument(
        "-o", "--input-owner-str", default="tduigou", help="Owner of the repository"
    )
    args = parser.parse_args()
    lint_shed_yaml(repository=args.input_repository_str, owner=args.input_owner_str)
