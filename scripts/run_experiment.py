import hydra
from omegaconf import DictConfig, OmegaConf
from hydra.core.hydra_config import HydraConfig
from hydra.utils import to_absolute_path

import transformers
from transformers import Text2TextGenerationPipeline
import logging
import os
import sys
import wandb
from pathlib import Path

from dataset import load_data, get_dataloader
from scripts.transformers_util import get_training_args, get_trainer, get_tokenizer
from evaluation import error_analysis, get_scores, encode_decode

sys.path.append('scripts')

log = logging.getLogger(__name__)

def run(config, run_name, sweep_name):
    log.info(f'Fixing random seed {config.random_seed}')
    transformers.trainer_utils.set_seed(config.random_seed)
    
    transformers.logging.disable_default_handler()
    
    log.info(OmegaConf.to_yaml(config))
    log.info('Running in: ' + os.getcwd())

    with wandb.init(project=config.wandb_project, name=run_name, reinit=True) as run:

        wandb.log(OmegaConf.to_container(config))
        wandb.log({'hydra_sweep' : sweep_name})
        wandb.log({'experiment_dir': os.getcwd()})

        training_args = get_training_args(config, report_to="wandb")
        tokenizer = get_tokenizer(config)

        train_df, val_df, test_df = load_data(
            to_absolute_path(config.data.cnf_tsv_path), 
            to_absolute_path(config.data.controls_tsv_path) if config.data.controls_tsv_path else None,
            sample_frac=config.get("sample", None))

        wandb.log({"n_ellipses_train": (~train_df.controls).sum()})
        wandb.log({"n_ellipses_dev": (~val_df.controls).sum()})
        wandb.log({"n_ellipses_test": (~test_df.controls).sum()})

        wandb.log({"n_controls_train" : train_df.controls.sum()})
        wandb.log({"n_controls_dev" : val_df.controls.sum()})
        wandb.log({"n_controls_test" : test_df.controls.sum()})

        train_dataset, val_dataset, test_dataset = get_dataloader(train_df, val_df, test_df, tokenizer)

        trainer = get_trainer(config, tokenizer, training_args, train_dataset, val_dataset)

        trainer.train()
        
        wandb.log({'best_cp' : trainer.state.best_model_checkpoint})

        pipeline = Text2TextGenerationPipeline(model=trainer.model, tokenizer=tokenizer, max_length=config.generation_max_length, device=0)

        def get_errors(sample):
            out = pipeline(list(sample.raw_sentence))
            gen = [o['generated_text'] for o in out]
            errors = error_analysis(gen, encode_decode(sample.full_resolution, tokenizer), encode_decode(sample.raw_sentence, tokenizer))
            return errors

        log.info("Running error analysis on dev set")
        errors_valid = get_errors(val_df)               
        valid_scores = get_scores(errors_valid, "eval")
        wandb.log(valid_scores)

        log.info("Running error analysis on test set")
        errors_test = get_errors(test_df)               
        test_scores = get_scores(errors_test, "test")
        wandb.log(test_scores)

        return valid_scores['eval/exact_match']


@hydra.main(config_path='..', config_name='experiment.yaml', version_base="1.2")
def main(config: DictConfig):
    hydra_conf = HydraConfig.get()
       
    run_name = Path(os.getcwd()).name
    
    if 'id' in hydra_conf.job:
        # In a sweep
        sweep_name = hydra_conf.sweep.dir
        log.info("Running in a parameter sweep")
    else:
        sweep_name = config.date_run
        
    log.info(f"Grouping by {sweep_name}")

    metric = run(config, run_name, sweep_name)

    return metric
    
if __name__ == "__main__":
    main()
