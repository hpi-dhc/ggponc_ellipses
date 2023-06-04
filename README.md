# Resolving Elliptical Compounds in German Medical Text

## Preparation

1. Get access to GGPONC following the instructions on the [project homepage](https://www.leitlinienprogramm-onkologie.de/projekte/ggponc-english/) and place the contents of `ellipses_2023_01_30` in the `data/ellipses` folder (or adapt the path in `experiment.yaml` to point to another directory)
2. Install Python dependencies `pip install -r requirements.txt`
3. Fill and adjust the config file [experiment.yaml](scripts/experiment.yaml). If you want to use any of the approaches that require OpenAI services you need to specify your API key in the config

## Notebooks

In `notebooks`, we provide the following Jupyter Notebooks to reproduce the results from the paper:

- [01_Dataset.ipynb](notebooks/01_Dataset.ipynb)
    - Corpus Statistics
- [02_Baseline_Aepli.ipynb](notebooks/02_Baseline_Aepli.ipynb)
    - Updated implementation of a rule-based baseline by [Aepli & Volk, 2012](https://link.springer.com/chapter/10.1007/978-3-642-40722-2_1)
    - In order to run the baseline, you also need to download the full GGPONC 2.0 corpus and place it under `data/ggponc_v2` (or adapt the path in `experiment.yaml` to point to your data folder) 
- [03_Generative.ipynb](notebooks/03_Generative.ipynb)
    - Generative approach using encoder-decoder Transformer-based model (default: mT5)
- [04_Zero_Shot.ipynb](notebooks/04_Zero_Shot.ipynb)
    - Zero-shot approach using ChatGPT (API key needed)
- [05_TopK.ipynb](notebooks/05_TopK.ipynb)
    - Multiple-choice approaches using generative model to generate the `k` most likely sentences and ChatGPT to choose the best one (API Key needed) 

## Running Generative Transformer Experiments with HuggingFace and Hydra

In `scripts`, we provide [Hydra](https://github.com/facebookresearch/hydra) configurations for the different  experiments with the best hyperparameters found through grid search.
To run such an experiment, do:
- `cd scripts`
- `python run_experiment.py -cn <experiment>.yaml cuda=<cuda devices>`
    - for instance: `python run_experiment.py -cn experiment.yaml cuda=0`

If you have installed and configured [Weights & Biases](https://wandb.ai/), it will automatically sync your runs.

To run a hyperparameter sweep, specify your desired paramters in [experiment.yaml](scripts/experiment.yaml) under `params` and pass the optiom `-m` to Hydra, e.g.:
- `python run_experiment.py -m experiment.yaml cuda=0`

## Citation

Accepted at BioNLP at ACL'23

TODO: Add Citation

