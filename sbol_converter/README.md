# sbol-converter -- Convert between SBOL3 and other genetic design formats.

The open-source software package sbol-converter is available here : https://github.com/SynBioDex/SBOL-utilities

## How to run sbol-converter wrapper tests

In order to execute tests on sbol-converter wrapper, you need to:

  - Connect to your galaxy instance in interactive mode:

  ```bash
    docker exec -it -u root galaxy_galaxy_1 bash
  ```
  - Copy all the contents of `test-data` folder into your own test-data directory which is located in your local galaxy instance : `/galaxy/test-data`. It contains all the input files and expected output files needed for the tests.

  - Install Planemo:
  You can see here the documentation for Planemo Installation : https://planemo.readthedocs.io/en/latest/installation.html
  Note that they recommand to install Planemo by setting up a virtual environment:

  ```bash
  python3 -m venv planemo
  . planemo/bin/activate
  pip install -U planemo
  ```

  - run the tests:

  ```bash
  planemo test --conda_channels conda-forge,bioconda --conda_debug tools/synbiocad-galaxy-wrappers/sbol-converter/wrap.xml
  ```

  You can also specify your data test directory with the following option: `--test_data DIRECTORY_DATA_TEST`

  IMPORTANT: Maybe you will need to remove CONDA from your PATH for the command `planemo test` to run correctly. To do that, you can edit this file `~/.bashrc`, comment this line `PATH="/root/anaconda3/bin:$PATH"` and save changes.

  Planemo will output an html test summary `tool_test_output.html`.