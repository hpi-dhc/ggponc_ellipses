# Resolving Elliptical Compounds in German Medical Text

This repository contains the code to reproduce results from the paper:

ToDo: Citation

## Preparation

1. Get access to GGPONC following the instructions on the [project homepage](https://www.leitlinienprogramm-onkologie.de/projekte/ggponc-english/) and place the the contents of the 2.0 release (v2.0_2022_03_24) in the `data` folder
2. Install Python dependencies `pip install -r requirements.txt`
3. Fill and adjust the config file [experiment.yaml](scripts/experiment.yaml). If you want to use any of the approaches that require OpenAI services you need to specify your API key in the config

## Notebooks

In `notebooks`, we provide the following Jupyter Notebooks to reproduce the results from the paper:

- [01_Dataset.ipynb](notebooks/01_Dataset.ipynb)
    - Corpus Statistics
- [02_Baseline_Aepli.ipynb](notebooks/02_Baseline_Aepli.ipynb)
    - Rule-based baseline implemented from Noëmi Aepli and Martin Volk. 2013. [Reconstructing complete lemmas for incomplete German compounds](https://link.springer.com/chapter/10.1007/978-3-642-40722-2_1). In Language Processing and Knowledge in
      the Web, pages 1–13. Springer
- [03_Generative.ipynb](notebooks/03_Generative.ipynb)
    - Generative approach using HuggingFace transformer-based model
- [04_Zero_Shot.ipynb](notebooks/04_Zero_Shot.ipynb)
    - Zero-shot approach using ChatGPT (API Key needed)
- [05_TopK.ipynb](notebooks/05_TopK.ipynb)
    - Two-fold approach using generative model to generate the *k* most likely sentences and ChatGPT to choose the best one (API Key needed) 

## Running Generative Transformer Experiments with HuggingFace and Hydra

In `scripts`, we provide [Hydra](https://github.com/facebookresearch/hydra) configurations for the different  experiments with the best hyperparameters found through grid search.
To run such an experiment, do:
- `cd scripts`
- `python run_experiment.py -cn <experiment>.yaml cuda=<cuda devices>`
    - for instance: `python run_experiment.py -cn experiment.yaml cuda=0`

If you have installed and configured [Weights & Biases](https://wandb.ai/), it will automatically sync your runs.

To run a hyperparameter sweep, specify your desired paramters in [experiment.yaml](scripts/experiment.yaml) under *params* and pass the optiom `-m` to Hydra, e.g.:
- `python run_experiment.py -m experiment.yaml cuda=0`

## Citing our Paper

ToDo: Citation

