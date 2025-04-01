#!/bin/bash

planemo test \
    --galaxy_source https://github.com/galaxyproject/galaxy \
    --galaxy_branch master \
    --galaxy_python_version 3.11 \
    --test_timeout 86400 \
    --biocontainers \
    --no_dependency_resolution \
    --no_conda_auto_init \
    --docker_extra_volume . \
    --test_output_json output.json \
    --test_output output.html \
    /opt/brsynth/galaxytools/tools/rrparser #calibrator.xml
