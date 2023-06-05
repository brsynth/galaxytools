#!/usr/bin/env python

import argparse
import json
import os
from typing import Any, Dict

from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.workflows import WorkflowClient


class GalaxyUploadWorkflow(object):
    def __init__(self, url: str, key: str):
        self.gi = GalaxyInstance(url=url, key=key)
        self.wc = WorkflowClient(galaxy_instance=self.gi)

    def upload(self, workflow: str) -> None:
        data = GalaxyUploadWorkflow.read_workflow(path=workflow)

        self.wc.update_workflow(
            workflow_id=data["uuid"],
            workflow=data["steps"],
            name=data["name"],
            annotation=data["annotation"],
            menu_entry=True,
            tags=data["tags"],
            published=True,
        )

    @classmethod
    def read_workflow(cls, path: str) -> Dict[str, Any]:
        with open(path) as fd:
            return json.load(fd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser_input = parser.add_argument_group("Input")
    parser_input.add_argument(
        "--galaxy_url",
        required=True,
        help="Url of Galaxy Instance",
    )
    parser_input.add_argument(
        "--galaxy_user_key",
        required=True,
        help="Token for Galaxy Instance",
    )
    parser_input.add_argument("--workflow", required=True, help="Worfklow file (.ga)")

    args = parser.parse_args()

    guw = GalaxyUploadWorkflow(url=args.galaxy_url, key=args.galaxy_user_key)
    guw.upload(workflow=args.workflow)
