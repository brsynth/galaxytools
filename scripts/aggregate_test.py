#!/usr/bin/env python
import argparse
import glob
import json
import os
from typing import List


def parse_artifacts(finputs: List[str]) -> int:
    """Return 0 if all are success 1 otherwise"""
    for finput in finputs:
        data = {}
        with open(finput) as fd:
            data = json.load(fd)
        if data["summary"]["num_errors"] + data["summary"]["num_failures"] > 0:
            raise ValueError("Error found !")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input-repository-str",
        default=os.getcwd(),
        help="Path to search recursively tool_test_output.json files",
    )
    parser.add_argument(
        "-f", "--input-tool-json", nargs="+", help="tool_test_output.json files"
    )
    args = parser.parse_args()

    finputs = []
    if args.input_tool_json:
        finputs = args.input_tool_json
    else:
        if not os.path.isdir(args.input_repository_str):
            parser.error("Input directory does not exist")
        finputs = glob.glob(
            os.path.join(args.input_repository_str, "**", "tool_test_output.json"),
            recursive=True,
        )
    if len(finputs) < 1:
        parser.error("No file provided")

    parse_artifacts(finputs)
