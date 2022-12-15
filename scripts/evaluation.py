import numpy as np
from evaluate import load
import nltk
from collections import Counter
import difflib
import pandas as pd

def encode_decode(series, tokenizer):
    return [tokenizer.decode(s, skip_special_tokens=True) for s in tokenizer(series.tolist())['input_ids']]

def error_analysis(predictions, gt_resolutions, original_sentences):
    d = difflib.Differ()
    res = []

    for pred_gen, true, sent in zip(predictions, gt_resolutions, original_sentences):
        entry = {'pred' : pred_gen, 'ground_truth' : true, 'original' : sent}    
        if pred_gen == true:
            entry['error_type'] = 'tp'
        elif pred_gen == sent:
            entry['error_type'] = 'fn'
        elif pred_gen != sent and sent == true:
            entry['error_type'] = 'fp'
        else:
            op_codes = difflib.SequenceMatcher(None, true, pred_gen).get_opcodes()
            counts = Counter([o[0] for o in op_codes])
            del counts["equal"]
            if len(counts) > 1:
                entry['error_type'] = 'complex'
            else:
                entry['error_type'] = list(counts.keys())[0]
        res.append(entry)
    return pd.DataFrame(res)

def relative_edit_distance(p, g, o):
    ed = nltk.edit_distance
    d = ed(p,g)
    k = ed(p,o)
    l = ed(o,g)
    if d == 0:
        return 1
    return 1 - (d / (k + l))

def get_scores(error_analysis_results, key):
    m = Metrics(['exact_match', 'google_bleu'], None)
    res = {}
    error_counts =  error_analysis_results.error_type.value_counts()
    for k in ['tp', 'fn', 'fp', 'replace', 'insert', 'delete', 'complex']:
        v = error_counts.loc[k] if k in error_counts.index else 0
        res[k] = v / len(error_analysis_results)
        res[f"{k}_abs"] = v
    relative_edit_score = error_analysis_results.apply(lambda r: relative_edit_distance(r['pred'], r['ground_truth'], r['original']), axis=1)
    res["edit_distance_rel"] = relative_edit_score.mean()
    res["exact_match"] = m.compute_exact_match(error_analysis_results['pred'], error_analysis_results['ground_truth'])
    res["gleu"] = m.compute_bleu(error_analysis_results['pred'], error_analysis_results['ground_truth'])
    #res["edit_distance_abs"] = m.compute_edit_distance(error_analysis_results['pred'], error_analysis_results['ground_truth'])
    return { f"{key}/{k}":v for k, v in res.items() } 
    
def postprocess_text(preds, labels):
    preds = [pred.strip() for pred in preds]
    labels = [[label.strip()] for label in labels]
    return preds, labels

class Metrics:

    def __init__(self, metric_names, tokenizer):
        metric_functions = {
            'exact_match': self.compute_exact_match,
            'google_bleu': self.compute_bleu,
            'edit_distance': self.compute_edit_distance
        }

        self.chosen_metrics = {name: load(name) for name in ['exact_match', 'google_bleu']}
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

    def compute_edit_distance(self, preds, refs):
        return np.mean([nltk.edit_distance(pred, ref) for pred, ref in zip(preds, refs)])
        