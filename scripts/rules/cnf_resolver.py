from typing import List
from .rule_modules import RuleModule, Suffix, Prefix

class RulebasedCNFResolver:

    rule_modules: List[RuleModule]

    def __init__(self, num_nn: int):
        self.rule_modules = [Suffix(num_nn), Prefix(num_nn)]

    def predict_all(self, sentences: List[str]):
        predictions = []
        for sentence in sentences:
            predictions.append(self.predict_single(sentence))
        return predictions

    def predict_single(self, sentence: str):
        for module in self.rule_modules:
            sentence = module.predict(sentence)
        return sentence