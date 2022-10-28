from typing import Dict, Union, Any, Optional, List, Tuple, Callable

import torch
from torch import nn
from torch.utils.data import Dataset
from transformers import Seq2SeqTrainer, PreTrainedModel, TrainingArguments, DataCollator, PreTrainedTokenizerBase, \
    EvalPrediction, TrainerCallback
from transformers.deepspeed import is_deepspeed_zero3_enabled


class Seq2SeqTrainerCustom(Seq2SeqTrainer):

    def __init__(
            self,
            length_multiplier: float = 1.5,
            model: Union[PreTrainedModel, nn.Module] = None,
            args: TrainingArguments = None,
            data_collator: Optional[DataCollator] = None,
            train_dataset: Optional[Dataset] = None,
            eval_dataset: Optional[Dataset] = None,
            tokenizer: Optional[PreTrainedTokenizerBase] = None,
            model_init: Callable[[], PreTrainedModel] = None,
            compute_metrics: Optional[Callable[[EvalPrediction], Dict]] = None,
            callbacks: Optional[List[TrainerCallback]] = None,
            optimizers: Tuple[torch.optim.Optimizer, torch.optim.lr_scheduler.LambdaLR] = (None, None),
            preprocess_logits_for_metrics: Callable[[torch.Tensor, torch.Tensor], torch.Tensor] = None,
    ):
        self.length_multiplier = length_multiplier
        super().__init__(
            model=model,
            args=args,
            data_collator=data_collator,
            train_dataset=train_dataset,
            eval_dataset=eval_dataset,
            tokenizer=tokenizer,
            model_init=model_init,
            compute_metrics=compute_metrics,
            callbacks=callbacks,
            optimizers=optimizers,
            preprocess_logits_for_metrics=preprocess_logits_for_metrics
        )

    def prediction_step(
        self,
        model: nn.Module,
        inputs: Dict[str, Union[torch.Tensor, Any]],
        prediction_loss_only: bool,
        ignore_keys: Optional[List[str]] = None,
    ) -> Tuple[Optional[float], Optional[torch.Tensor], Optional[torch.Tensor]]:

        if not self.args.predict_with_generate or prediction_loss_only:
            return super().prediction_step(
                model, inputs, prediction_loss_only=prediction_loss_only, ignore_keys=ignore_keys
            )

        has_labels = "labels" in inputs
        inputs = self._prepare_inputs(inputs)

        # XXX: adapt synced_gpus for fairscale as well
        gen_kwargs = self._gen_kwargs.copy()
        gen_kwargs["max_length"] = (
            gen_kwargs["max_length"] if gen_kwargs.get("max_length") is not None else self.model.config.max_length
        )
        gen_kwargs["num_beams"] = (
            gen_kwargs["num_beams"] if gen_kwargs.get("num_beams") is not None else self.model.config.num_beams
        )
        default_synced_gpus = True if is_deepspeed_zero3_enabled() else False
        gen_kwargs["synced_gpus"] = (
            gen_kwargs["synced_gpus"] if gen_kwargs.get("synced_gpus") is not None else default_synced_gpus
        )

        if "attention_mask" in inputs:
            gen_kwargs["attention_mask"] = inputs.get("attention_mask", None)
        if "global_attention_mask" in inputs:
            gen_kwargs["global_attention_mask"] = inputs.get("global_attention_mask", None)

        # prepare generation inputs
        # some encoder-decoder models can have varying encoder's and thus
        # varying model input names
        if hasattr(self.model, "encoder") and self.model.encoder.main_input_name != self.model.main_input_name:
            generation_inputs = inputs[self.model.encoder.main_input_name]
        else:
            generation_inputs = inputs[self.model.main_input_name]

        """generated_tokens = self.model.generate(
            generation_inputs,
            **gen_kwargs,
        )"""

        tokens = []
        for input_sentence in generation_inputs:
            length = input_sentence.shape[0]
            gen = self.model.generate(
                input_sentence,
                min_length=length,
                max_length=self.length_multiplier*length,
                **gen_kwargs)

            tokens.append(
                self._pad_tensors_to_max_len(gen, gen_kwargs["max_length"])
            )

        generated_tokens = torch.stack(tokens)


        # in case the batch is shorter than max length, the output should be padded
        if generated_tokens.shape[-1] < gen_kwargs["max_length"]:
            generated_tokens = self._pad_tensors_to_max_len(generated_tokens, gen_kwargs["max_length"])

        with torch.no_grad():
            with self.compute_loss_context_manager():
                outputs = model(**inputs)
            if has_labels:
                if self.label_smoother is not None:
                    loss = self.label_smoother(outputs, inputs["labels"]).mean().detach()
                else:
                    loss = (outputs["loss"] if isinstance(outputs, dict) else outputs[0]).mean().detach()
            else:
                loss = None

        if self.args.prediction_loss_only:
            return (loss, None, None)

        if has_labels:
            labels = inputs["labels"]
            if labels.shape[-1] < gen_kwargs["max_length"]:
                labels = self._pad_tensors_to_max_len(labels, gen_kwargs["max_length"])
        else:
            labels = None

        return (loss, generated_tokens, labels)