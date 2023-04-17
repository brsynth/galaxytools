# RRulesParser Galaxy wrapper

The open-source software package RRulesParser is available here : https://github.com/brsynth/RRParser

## Installing

* Create a new section in the Galaxy config file `config/tool_conf.xml` and paste the following:
```
<section id="sbc-rs" name="SynBioCAD RetroSynthesis">
    <tool file="synbiocad-galaxy-wrappers/RRulesParser/wrap.xml" />
</section>
```

* Make sure that the following entry exists under Galaxy's destination tag in `config/job_conf.xml`:
```
<destination id="local" runner="local" />
```

And that the destination of the tool is referred under the tools tag in `config/job_conf.xml`:

```
<tool id="RRulesParser" destination="local" />
```

## How to run RRulesParser wrapper tests

In order to execute tests on RRulesParser wrapper, you need to:

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

  Upgrade pip if needed.

  - run the tests:

  ```bash
  planemo test --conda_channels conda-forge tools/synbiocad-galaxy-wrappers/RRulesParser/wrap.xml --biocontainers --no_conda_auto_init
  ```

  IMPORTANT: Maybe you will need to remove CONDA from your PATH for the command `planemo test` to run correctly. To do that, you can edit this file `~/.bashrc`, comment this line `PATH="/root/anaconda3/bin:$PATH"` and save changes.

  Planemo will output an html test summary `tool_test_output.html`.

## Built With

* [Galaxy](https://galaxyproject.org) - The Galaxy project


## Authors

* Melchior du Lac
* Joan Hérisson
* Thomas Duigou

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
