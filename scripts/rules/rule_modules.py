from abc import ABC, abstractmethod
from typing import Tuple
from spacy.matcher import Matcher
import spacy
from charsplit import Splitter
import re

class RuleModule(ABC):

    def __init__(self, num_nn: int):
        self.current_sentence = None

        self.doc = None
        self.splitter = Splitter()
        self.nlp = spacy.load('de_core_news_lg')

        self.num_nn = num_nn

    def predict(self, sentence: str):
        self.current_sentence = sentence
        for num in range(self.num_nn, 0, -1):
            matches = self.identify_elliptic_phrases(num)
            for match_id, start, end in matches:
                span = self.doc[start:end]
                fixed_phrase = self.fix_incomplete_part(span.text, (start, end))
                self.current_sentence = self.current_sentence.replace(span.text, fixed_phrase)
        return self.current_sentence

    def identify_elliptic_phrases(self, num: int):
        self.doc = self.nlp(self.current_sentence)
        matcher = Matcher(self.nlp.vocab)
        pattern = self.generate_pattern(num)
        matcher.add('PossibleEllipses', pattern)
        matches = matcher(self.doc)
        return matches

    @abstractmethod
    def fix_incomplete_part(self, phrase: str, indices):
        pass

    @abstractmethod
    def generate_pattern(self, num_nn: int):
        pass


class Suffix(RuleModule):

    def __init__(self, num_nn: int):
        super().__init__(num_nn)

        self.matcher_nn = Matcher(self.nlp.vocab)
        self.matcher_nn.add('WordToSplit', [[{'TAG': {'IN': ['NN', 'NE', 'VAFIN', 'FM']}}]])

        self.matcher_trunc = Matcher(self.nlp.vocab)
        self.matcher_trunc.add('WordToAdd', [[{'TAG': 'TRUNC'}]])

    def generate_pattern(self, num_nn: int):
        pattern = []
        for i in range(0, num_nn):
            pattern.append({'TAG': 'TRUNC'})
            pattern.append({'TAG': {'IN': ['ADV', 'ADJA', 'ADJD']}, 'OP': '*'})
            pattern.append({'TAG': {'IN': ['KON', '$,', '$(']}, 'OP': '+'})  # more exhaustive tracking of punctuation
        pattern.append({'TAG': {'IN': ['ADV', 'ADJA', 'ADJD']}, 'OP': '*'})
        pattern.append({'TAG': {'IN': ['NN', 'NE', 'NNE', 'XY']}})  # proper nouns
        return [pattern]

    def fix_incomplete_part(self, phrase: str, indices: Tuple[int, int]):
        word = self.identify_complete_part(indices)
        if word == 'error':
            return phrase
        words = self.identify_incomplete_part(indices)
        fixed_words = [split.replace('-', word.lower()) for split in words]
        for old, new in zip(words, fixed_words):
            phrase = phrase.replace(old, new)
        return phrase

    def identify_complete_part(self, indices: Tuple[int, int]):
        matches = self.matcher_nn(self.doc[indices[0]:indices[1]])
        words = []
        for match_id, start, end in matches:
            span = self.doc[indices[0]:indices[1]][start:end]
            words.append(span.text)
        if len(words) != 1:
            return 'error'
        split = self.splitter.split_compound(words[0])[0][1:]
        if len(split) > 2:
            pass  # maybe consider as error
        return split[1]

    def identify_incomplete_part(self, indices: Tuple[int, int]):
        matches = self.matcher_trunc(self.doc[indices[0]:indices[1]])
        words = []
        for match_id, start, end in matches:
            span = self.doc[indices[0]:indices[1]][start:end]
            words.append(span.text)
        return words


class Prefix(RuleModule):
    def __init__(self, num_nn):
        super().__init__(num_nn)

        self.matcher_nn = Matcher(self.nlp.vocab)
        self.matcher_nn.add('WordToSplit', [[{'TAG': {'IN': ['NN', 'NE', 'VAFIN', 'FM']}}]])

    def generate_pattern(self, num_nn: int):
        pattern = []
        for i in range(0, num_nn):
            pattern.append({'TAG': {'IN': ['NN', 'NE', 'NNE', 'XY']}})
            pattern.append({'TAG': {'IN': ['KON', '$,', '$(']}})  # more exhaustive tracking of punctuation
        pattern.append({'TEXT': {'IN': ['-', '–']}, 'OP': '?'})
        pattern.append({'TAG': {'IN': ['NN', 'NE', 'NNE', 'XY', 'VVFIN']}})
        return [pattern]

    def fix_incomplete_part(self, phrase: str, indices: Tuple[int, int]):
        split_words = self.identify_incomplete_part(phrase)
        if not split_words:
            return phrase
        complete_words = self.identify_complete_part(indices)
        fixed_words = [split.replace('-', complete_words).replace('–', complete_words) for split in split_words]
        for old, new in zip(split_words, fixed_words):
            phrase = phrase.replace(old, new)
        return phrase

    def identify_complete_part(self, indices: Tuple[int, int]):
        matches = self.matcher_nn(self.doc[indices[0]:indices[1]])
        words = []
        for match_id, start, end in matches:
            span = self.doc[indices[0]:indices[1]][start:end]
            words.append(span.text)
        if len(words) != 1:
            pass  # maybe consider as error
        split = self.splitter.split_compound(words[0])[0][1:]
        if len(split) > 2:
            pass  # maybe consider as error
        return split[0]

    def identify_incomplete_part(self, phrase):
        return re.findall(r' [-–][a-zA-ZäöüÄÖÜß]+', phrase)