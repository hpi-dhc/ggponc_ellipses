import numpy as np
from evaluate import load
import nltk

def postprocess_text(preds, labels):
    preds = [pred.strip() for pred in preds]
    labels = [[label.strip()] for label in labels]
    return preds, labels


class Metrics:

    def __init__(self, metric_names, tokenizer):
        metric_functions = {
            'exact_match': self.compute_exact_match,
            'google_bleu': self.compute_bleu,
            'meteor': self.compute_meteor,
            'edit_distance': self.compute_edit_distance
        }

        self.chosen_metrics = {name: load(name) for name in ['exact_match', 'google_bleu', 'meteor']}
        self.chosen_metric_functions = {name: metric_functions[name] for name in metric_names}
        self.tokenizer = tokenizer

    def compute_metrics(self, eval_preds):
        predictions, labels = eval_preds

        if isinstance(predictions, tuple):
            predictions = predictions[0]

        labels = np.where(labels != -100, labels, self.tokenizer.pad_token_id)

        decoded_preds = self.tokenizer.batch_decode(predictions, skip_special_tokens=True)
        decoded_labels = self.tokenizer.batch_decode(labels, skip_special_tokens=True)

        return {name: self.chosen_metric_functions[name](decoded_preds, decoded_labels) for name in self.chosen_metric_functions.keys()}

    def compute_exact_match(self, preds, refs):
        return self.chosen_metrics["exact_match"].compute(predictions=preds, references=refs)['exact_match']

    def compute_bleu(self, preds, refs):
        preds, refs = postprocess_text(preds, refs)
        return self.chosen_metrics["google_bleu"].compute(predictions=preds, references=refs)['google_bleu']

    def compute_meteor(self, preds, refs):
        #preds, refs = postprocess_text(preds, refs)
        return self.chosen_metrics["meteor"].compute(predictions=preds, references=refs)['meteor']

    def compute_edit_distance(self, preds, refs):
        return np.mean([nltk.edit_distance(pred, ref) for pred, ref in zip(preds, refs)])
        