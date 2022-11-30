# -*- coding: utf-8 -*-

# Author: Noëmi Aepli, Institute of Computational Linguistics, University of Zurich

"""Script to fix TRUNC-lemmas of elliptically coordinated compounds in T+B corpus

The python program (used as a script) looks for "TRUNC + compound"-patterns, analyses
them and substitutes TRUNC-lemma with the intended compound.

If there is no solution to be found or the case is not decidable, nothing happens 
(i.e. input = output).

TRUNC stands for truncation and is the POS-tag for the first part of a compound 
according to STTS (Stuttgard Tübingen TagSet)

Dependencies:
gertwolList: file with Gertwol-analysed lemmas of compounds
wordFrequencyList: file with all words and corresponding frequencies of T+B corpus

Usage:
python ellipticCompounds.py FileIn FileOut gertwolList wordFrequencyList

Classes:
EllipticCompound()

"""

from xml.etree import cElementTree
from collections import defaultdict
import re
import codecs
import operator
import sys


# TRUNC = word with suspended hyphen (normal case: 1st word, reversed case: 2nd word)
# compound (normal case: 2nd word, reversed case: 1st word)

class EllipticCompound():

    def __init__(self, xmlIn, xmlOut, gertwol, freqDict):
        self.gertwolDict = self.createGertwolDict(gertwol)
        self.document = cElementTree.parse(xmlIn)
        self.words = self.document.findall('.//w')
        self.frequencyDict = self.createFrequencyDict(freqDict)
        self.undecDict = defaultdict(int)
        self.lemmaSuff = {}
        self.count = 0
        self.compList = []
        self.taggedList = []
        with codecs.open('../data/baseline_files/compoundListTagged.txt', 'r', 'utf-8') as fileTagged:  # list of tagged words
            for lines in fileTagged:
                self.taggedList.append(lines.split())

    # create dictionary from words in gertwolList
    # format list/dict: WORD \t GERTWOL-VERSION [e.g. Sportkletterroute \t Sport#klett~er#route]
    def createGertwolDict(self, gertwol):
        gertwolList = codecs.open(gertwol, "r", encoding="utf-8")
        gertwolDict = {}
        for line in gertwolList:
            line = re.split("\t", line)
            gertwolDict[line[0]] = line[1]
        return gertwolDict

    # for each word in corpus store frequency
    # format: WORD: FREQUENCY [e.g. Schneegrenze: 2]
    def createFrequencyDict(self, freqDict):
        frequencyList = codecs.open(freqDict, "r", encoding="utf-8")
        frequencyDict = defaultdict(int)
        for line in frequencyList:
            line = re.split("\t", line)
            frequencyDict[line[0]] = line[1]
        return frequencyDict

    # find different patterns of elliptic compounds
    def findPattern(self):
        current = 0
        for word in self.words:
            if len(self.words) > (current + 2):

                # (1) TRUNC $, TRUNC $, TRUNC KON NN/NE/ADJA
                if (word.get('pos') == "TRUNC" and word.text[-1] == "-" and self.words[current + 1].get(
                        'pos') == "$," and self.words[current + 2].get('pos') == "TRUNC" and
                    self.words[current + 2].text[-1] == "-" and self.words[current + 3].get('pos') == "$," and
                    self.words[current + 4].get('pos') == "TRUNC" and self.words[current + 4].text[-1] == "-") and (
                        self.words[current + 5].get('pos') == "KON") and (
                        self.words[current + 6].get('pos') == "NN" or self.words[current + 6].get('pos') == "NE" or
                        self.words[current + 6].get('pos') == "ADJA") and self.words[current + 6].text[0] != "-":
                    truncWord1 = self.words[current].text.strip()
                    truncWord2 = self.words[current + 2].text.strip()
                    truncWord3 = self.words[current + 4].text.strip()
                    newCompoundLemma = self.substituteCompoundLemma(self.words[current + 6], current + 6)
                    segments1 = self.analyseCompoundLemma(newCompoundLemma, truncWord1)
                    segments2 = self.analyseCompoundLemma(newCompoundLemma, truncWord2)
                    segments3 = self.analyseCompoundLemma(newCompoundLemma, truncWord3)
                    newTruncLemma1 = self.findTruncLemma(truncWord1, segments1, current)
                    newTruncLemma2 = self.findTruncLemma(truncWord2, segments2, current + 2)
                    newTruncLemma3 = self.findTruncLemma(truncWord3, segments3, current + 4)
                    self.substituteLemma(newTruncLemma1, current)
                    self.substituteLemma(newTruncLemma2, current + 2)
                    self.substituteLemma(newTruncLemma3, current + 4)

                # (2) TRUNC $, TRUNC KON NN/NE/ADJA
                elif (word.get('pos') == "TRUNC" and word.text[-1] == "-" and self.words[current + 1].get(
                        'pos') == "$," and self.words[current + 2].get('pos') == "TRUNC" and
                      self.words[current + 2].text[-1] == "-") and (self.words[current + 3].get('pos') == "KON") and (
                        self.words[current + 4].get('pos') == "NN" or self.words[current + 4].get('pos') == "NE" or
                        self.words[current + 4].get('pos') == "ADJA") and self.words[current + 4].text[0] != "-":
                    truncWord1 = self.words[current].text.strip()
                    truncWord2 = self.words[current + 2].text.strip()
                    newCompoundLemma = self.substituteCompoundLemma(self.words[current + 4], current + 4)
                    segments1 = self.analyseCompoundLemma(newCompoundLemma, truncWord1)
                    segments2 = self.analyseCompoundLemma(newCompoundLemma, truncWord2)
                    newTruncLemma1 = self.findTruncLemma(truncWord1, segments1, current)
                    newTruncLemma2 = self.findTruncLemma(truncWord2, segments2, current + 2)
                    self.substituteLemma(newTruncLemma1, current)
                    self.substituteLemma(newTruncLemma2, current + 2)


                # TRUNC KON
                elif word.get('pos') == "TRUNC" and word.text[-1] == "-" and self.words[current + 1].get(
                        'pos') == "KON":  # (3) trunc kon -WORD
                    if self.words[current + 2].text[0] == "-":
                        trunc1 = self.words[current].text.strip()
                        trunc2 = self.words[current + 2].text.strip()

                        # word 1 = trunc1 + compound1[1]
                        # word 2 = compound2[0] + trunc2
                        newCompoundLemma1 = self.substituteCompoundLemma(self.words[current + 2], current + 2)
                        newCompoundLemma2 = self.substituteCompoundLemma(self.words[current], current)
                        segments1 = self.analyseCompoundLemma(newCompoundLemma1, trunc1)
                        segments2 = self.rvsdAnalyseCompoundLemma(newCompoundLemma2, trunc2)
                        newTruncLemma1 = self.findTruncLemma(trunc1, segments1, current + 2)
                        self.checkSegment(newTruncLemma1, current + 2, 2)
                        newTruncLemma2 = self.rvsdFindTruncLemma(trunc2, segments2, current)
                        self.rvsdCheckSegment(newTruncLemma2, current + 2, 2)
                        self.substituteLemma(newTruncLemma1, current)
                        self.substituteLemma(newTruncLemma2, current + 2)

                    # (4) trunc kon APPR ART NN/NE
                    elif self.words[current + 2].get('pos') == "APPR" and self.words[current + 3].get(
                            'pos') == "ART" and (
                            self.words[current + 4].get('pos') == "NN" or self.words[current + 4].get('pos') == "NE"):
                        decideDict = defaultdict(int)
                        truncWord = self.words[current].text.strip()

                        # TRUNC + NN[1]
                        newCompoundLemma1 = self.substituteCompoundLemma(self.words[current + 4], current + 4)
                        segments1 = self.analyseCompoundLemma(newCompoundLemma1, truncWord)
                        trunc1 = self.findTruncLemma(truncWord, segments1, current + 4)
                        if trunc1[-1] != "-":  # if segments != [] -> no solution
                            newTruncLemma1 = re.sub(r'[\\, ~, #, |]', "", trunc1)
                            decideDict[(newTruncLemma1, trunc1)] = self.frequencyDict[newTruncLemma1]

                        # TRUNC + APPR
                        newCompoundLemma2 = self.substituteCompoundLemma(self.words[current + 2], current + 2)
                        segments2 = self.adjAnalyseCompoundLemma(newCompoundLemma2, truncWord)
                        trunc2 = self.findTruncLemma(truncWord, segments2, current + 2)
                        if trunc2[-1] != "-":  # if segments != [] -> no solution
                            newTruncLemma2 = re.sub(r'[\\, ~, #, |]', "", trunc2)
                            decideDict[(newTruncLemma2, trunc2)] = self.frequencyDict[newTruncLemma2]

                        try:
                            if max(decideDict.items(), key=operator.itemgetter(1))[1] != 0:
                                best = max(decideDict.items(), key=operator.itemgetter(1))[0]
                                resultList = [max(decideDict.items(), key=operator.itemgetter(1))[0]]
                                resultTruncLemma = resultList[0][1]
                            else:  # frequency of most frequent possibility == 0
                                resultTruncLemma = "*undecidable*"
                                self.undecDict["findPattern_4"] += 1
                        except:  # no dict because word smaller than 6 characters -> no segments
                            resultTruncLemma = "*undecidable*"

                        if resultTruncLemma == "*undecidable*":
                            for entries in self.taggedList:
                                if re.sub(r'[\\, ~, #, |, +]', "", trunc1.lower()) in entries[
                                    1].lower() and '+' in trunc1:
                                    self.checkSegment(trunc1, current + 4, 4)
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc1))
                                    break
                                elif re.sub(r'[\\, ~, #, |, +]', "", trunc2.lower()) in entries[
                                    1].lower() and '+' in trunc2:
                                    self.checkSegment(trunc2, current + 2, 2)
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc2))
                                    break
                            if '+' not in self.words[current].get('lemma'):
                                if '+' in trunc1 and '+' in trunc2:
                                    tr1 = re.sub(r'[\\, ~, |, +, #]', "", trunc1.split('+')[1])
                                    tr2 = re.sub(r'[\\, ~, |, +, #]', "", trunc2.split('+')[1])
                                    for entries in self.taggedList:
                                        if tr1.lower() == entries[1].lower() or re.sub(r'[\\, ~, |, +, #]', "",
                                                                                       trunc1.lower()) == entries[
                                            1].lower():
                                            trunc2 = ''
                                            break
                                        elif tr2.lower() == entries[1].lower() or re.sub(r'[\\, ~, |, +, #]', "",
                                                                                         trunc2.lower()) == entries[
                                            1].lower():
                                            trunc1 = ''
                                            break
                                if '+' in trunc1:
                                    if trunc1.split('+')[0].endswith('-'):
                                        trunc1 = trunc1.split('+')[0][:-1] + '+' + trunc1.split('+')[1]
                                    if trunc1.split('+')[1].startswith('-'):
                                        trunc1 = trunc1.split('+')[0] + '+' + trunc1.split('+')[1][1:]
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc1))
                                    if '#' not in self.words[current + 4].get('lemma'):
                                        tempInd = self.words[current + 4].get('lemma').rfind(
                                            re.sub(r'[\\, ~, |]', "", trunc1.split('+')[1][:3]))
                                        if tempInd != -1:
                                            if self.words[current + 4].get('lemma').endswith('-'):
                                                self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                        self.words[current + 4].get(
                                                                                            'lemma')[
                                                                                        :(tempInd - 1)] + '#' +
                                                                                        trunc1.split('+')[1]))
                                            else:
                                                self.words[current + 4].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 4].get(
                                                                                                'lemma')[
                                                                                            :tempInd] + '#' +
                                                                                            trunc1.split('+')[1]))
                                elif '+' in trunc2:
                                    if trunc2.split('+')[0].endswith('-'):
                                        trunc2 = trunc2.split('+')[0][:-1] + '+' + trunc2.split('+')[1]
                                    if trunc2.split('+')[1].startswith('-'):
                                        trunc2 = trunc2.split('+')[0] + '+' + trunc2.split('+')[1][1:]
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc2))
                                    if '#' not in self.words[current + 2].get('lemma'):
                                        tempInd = self.words[current + 2].get('lemma').rfind(
                                            re.sub(r'[\\, ~, |]', "", trunc2.split('+')[1][:3]))
                                        if tempInd != -1:
                                            if self.words[current + 2].get('lemma').endswith('-'):
                                                self.words[current + 2].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 2].get(
                                                                                                'lemma')[
                                                                                            :(tempInd - 1)] + '#' +
                                                                                            trunc2.split('+')[1]))
                                            else:
                                                self.words[current + 2].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 2].get(
                                                                                                'lemma')[
                                                                                            :tempInd] + '#' +
                                                                                            trunc2.split('+')[1]))
                        # self.generateList(trunc1, re.sub(r'[\\, ~, #, |, +]',"", trunc1))
                        # self.generateList(trunc2, re.sub(r'[\\, ~, #, |, +]',"", trunc2))
                        self.substituteLemma(resultTruncLemma, current)

                    # (5) trunc kon ADJA/ADJD/ADV/VVFIN/VVIZU/VVINF/VVPP/CARD/APPR !NN/NE
                    elif (self.words[current + 2].get('pos') == "ADJA" or self.words[current + 2].get(
                            'pos') == "ADJD" or self.words[current + 2].get('pos') == "ADV" or self.words[
                              current + 2].get('pos') == "CARD" or self.words[current + 2].get('pos') == "VVFIN" or
                          self.words[current + 2].get('pos') == "VVPP" or self.words[current + 2].get(
                                    'pos') == "VVIZU" or self.words[current + 2].get('pos') == "VVINF" or self.words[
                              current + 2].get('pos') == "APPR") and (
                            self.words[current + 3].get('pos') != "NN" and self.words[current + 3].get('pos') != "NE"):
                        truncWord = self.words[current].text.strip()
                        newCompoundLemma = self.substituteCompoundLemma(self.words[current + 2], current + 2)
                        segments = self.adjAnalyseCompoundLemma(newCompoundLemma, truncWord)
                        newTruncLemma = self.findTruncLemma(truncWord, segments, current + 2)
                        self.checkSegment(newTruncLemma, current + 2, 2)
                        self.substituteLemma(newTruncLemma, current)

                    # (6) trunc kon ADJA/CARD/ADV NN/NE
                    elif (self.words[current + 3].get('pos') == "NN" or self.words[current + 3].get(
                            'pos') == "NE") and (
                            self.words[current + 2].get('pos') == "ADJA" or self.words[current + 2].get(
                            'pos') == "CARD" or self.words[current + 2].get('pos') == "ADV"):
                        decideDict = defaultdict(int)
                        truncWord = self.words[current].text.strip()

                        # [1] TRUNC-NN
                        compoundLemma1 = self.words[current + 3].text
                        trunc1 = self.mergeTruncNn(truncWord, compoundLemma1, current)
                        newTruncLemma1 = re.sub(r'[\\, ~, #, |]', "", trunc1)
                        decideDict[(newTruncLemma1, trunc1)] = self.frequencyDict[newTruncLemma1]

                        # [2] TRUNC + ADJA[1]
                        newCompoundLemma2 = self.substituteCompoundLemma(self.words[current + 2], current + 2)
                        segments2 = self.adjAnalyseCompoundLemma(newCompoundLemma2, truncWord)
                        trunc2 = self.findTruncLemma(truncWord, segments2, current + 2)
                        if trunc2[-1] != "-":  # if segments != [] -> no solution
                            newTruncLemma2 = re.sub(r'[\\, ~, #, |]', "", trunc2)
                            decideDict[(newTruncLemma2, trunc2)] = self.frequencyDict[newTruncLemma2]

                        # [3] TRUNC + NN[1]
                        newCompoundLemma3 = self.substituteCompoundLemma(self.words[current + 3], current + 3)
                        segments3 = self.analyseCompoundLemma(newCompoundLemma3, truncWord)
                        trunc3 = self.findTruncLemma(truncWord, segments3, current + 3)
                        if trunc3[-1] != "-":  # if segments != [] -> no solution
                            newTruncLemma3 = re.sub(r'[\\, ~, #, |]', "", trunc3)
                            decideDict[(newTruncLemma3, trunc3)] = self.frequencyDict[newTruncLemma3]

                        try:
                            if max(decideDict.items(), key=operator.itemgetter(1))[1] != 0:
                                best = max(decideDict.items(), key=operator.itemgetter(1))[0]
                                resultList = [max(decideDict.items(), key=operator.itemgetter(1))[0]]
                                resultTruncLemma = resultList[0][1]
                            else:  # frequency of most frequent possibility == 0
                                resultTruncLemma = "*undecidable*"
                                self.undecDict["findPattern_6"] += 1
                        except:  # no dict because word smaller than 6 characters -> no segments
                            resultTruncLemma = "*undecidable*"

                        if resultTruncLemma == "*undecidable*":
                            for entries in self.taggedList:
                                if re.sub(r'[\\, ~, #, |, +]', "", trunc2.lower()) in entries[1].lower() and '+' in trunc2:
                                    self.checkSegment(trunc2, current + 2, 2)
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc2))
                                    break
                                elif re.sub(r'[\\, ~, #, |, +]', "", trunc3.lower()) in entries[
                                    1].lower() and '+' in trunc3:
                                    self.checkSegment(trunc3, current + 3, 3)
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc3))
                                    break
                            if '+' not in self.words[current].get('lemma'):
                                if '+' in trunc2 and '+' in trunc3:
                                    tr2 = re.sub(r'[\\, ~, |, +, #]', "", trunc2.split('+')[1])
                                    tr3 = re.sub(r'[\\, ~, |, +, #]', "", trunc3.split('+')[1])
                                    for entries in self.taggedList:
                                        if tr2.lower() == entries[1].lower() or re.sub(r'[\\, ~, |, +, #]', "",
                                                                                       trunc2.lower()) == entries[
                                            1].lower():
                                            trunc3 = ''
                                            break
                                        elif tr3.lower() == entries[1].lower() or re.sub(r'[\\, ~, |, +, #]', "",
                                                                                         trunc3.lower()) == entries[
                                            1].lower():
                                            trunc2 = ''
                                            break
                                if '+' in trunc2:
                                    if trunc2.split('+')[0].endswith('-'):
                                        trunc2 = trunc2.split('+')[0][:-1] + '+' + trunc2.split('+')[1]
                                    if trunc2.split('+')[1].startswith('-'):
                                        trunc2 = trunc2.split('+')[0] + '+' + trunc2.split('+')[1][1:]
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc2))
                                    if '#' not in self.words[current + 2].get('lemma'):
                                        tempInd = self.words[current + 2].get('lemma').rfind(
                                            re.sub(r'[\\, ~, |]', "", trunc2.split('+')[1][:4]))
                                        if tempInd != -1:
                                            if self.words[current + 2].get('lemma').endswith('-'):
                                                self.words[current + 2].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 2].get(
                                                                                                'lemma')[
                                                                                            :(tempInd - 1)] + '#' +
                                                                                            trunc2.split('+')[1]))
                                            else:
                                                self.words[current + 2].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 2].get(
                                                                                                'lemma')[
                                                                                            :tempInd] + '#' +
                                                                                            trunc2.split('+')[1]))
                                elif '+' in trunc3:
                                    if trunc3.split('+')[0].endswith('-'):
                                        trunc3 = trunc3.split('+')[0][:-1] + '+' + trunc3.split('+')[1]
                                    if trunc3.split('+')[1].startswith('-'):
                                        trunc3 = trunc3.split('+')[0] + '+' + trunc3.split('+')[1][1:]
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc3))
                                    if '#' not in self.words[current + 3].get('lemma'):
                                        tempInd = self.words[current + 3].get('lemma').rfind(
                                            re.sub(r'[\\, ~, |]', "", trunc3.split('+')[1][:4]))
                                        if tempInd != -1:
                                            if self.words[current + 3].get('lemma').endswith('-'):
                                                self.words[current + 3].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 3].get(
                                                                                                'lemma')[
                                                                                            :(tempInd - 1)] + '#' +
                                                                                            trunc3.split('+')[1]))
                                            else:
                                                self.words[current + 3].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 3].get(
                                                                                                'lemma')[
                                                                                            :tempInd] + '#' +
                                                                                            trunc3.split('+')[1]))
                        # self.generateList(trunc2, re.sub(r'[\\, ~, #, |, +]',"", trunc2))
                        # self.generateList(trunc3, re.sub(r'[\\, ~, #, |, +]',"", trunc3))
                        self.substituteLemma(resultTruncLemma, current)

                    # (7) trunc kon ART/APPR/APPRART NN/NE
                    elif (self.words[current + 2].get('pos') == "ART" or self.words[current + 2].get('pos') == "APPR" or
                          self.words[current + 2].get('pos') == "APPRART") and (
                            self.words[current + 3].get('pos') == "NN" or self.words[current + 3].get('pos') == "NE"):
                        truncWord = self.words[current].text.strip()
                        newCompoundLemma = self.substituteCompoundLemma(self.words[current + 3], current + 3)
                        segments = self.analyseCompoundLemma(newCompoundLemma, truncWord)
                        newTruncLemma = self.findTruncLemma(truncWord, segments, current + 3)
                        self.checkSegment(newTruncLemma, current + 3, 3)
                        self.substituteLemma(newTruncLemma, current)

                    # (8) trunc kon ART/CARD/ADV/APPR/APPRART ADJA NN/NE
                    elif (self.words[current + 2].get('pos') == "ART" or self.words[current + 2].get('pos') == "CARD" or
                          self.words[current + 2].get('pos') == "ADV" or self.words[current + 2].get('pos') == "APPR" or
                          self.words[current + 2].get('pos') == "APPRART") and (
                            self.words[current + 3].get('pos') == "ADJA") and (
                            self.words[current + 4].get('pos') == "NN" or self.words[current + 4].get('pos') == "NE"):
                        decideDict = defaultdict(int)
                        truncWord = self.words[current].text.strip()

                        # TRUNC + NN
                        newCompoundLemma1 = self.words[current + 4].text
                        trunc1 = self.mergeTruncNn(truncWord, newCompoundLemma1, current)
                        newTruncLemma1 = re.sub(r'[\\, ~, #, |]', "", trunc1)
                        decideDict[(newTruncLemma1, trunc1)] = self.frequencyDict[newTruncLemma1]

                        # TRUNC + NN[1]
                        newCompoundLemma2 = self.substituteCompoundLemma(self.words[current + 4], current + 4)
                        segments2 = self.analyseCompoundLemma(newCompoundLemma2, truncWord)
                        trunc2 = self.findTruncLemma(truncWord, segments2, current + 4)
                        if trunc2[-1] != "-":  # if segments != [] -> no solution
                            newTruncLemma2 = re.sub(r'[\\, ~, #, |]', "", trunc2)
                            decideDict[(newTruncLemma2, trunc2)] = self.frequencyDict[newTruncLemma2]

                        # TRUNC + ADJ[1]
                        newCompoundLemma3 = self.substituteCompoundLemma(self.words[current + 3], current + 3)
                        segments3 = self.adjAnalyseCompoundLemma(newCompoundLemma3, truncWord)
                        trunc3 = self.findTruncLemma(truncWord, segments3, current + 3)
                        if trunc3[-1] != "-":  # if segments != [] -> no solution
                            newTruncLemma3 = re.sub(r'[\\, ~, #, |]', "", trunc3)
                            decideDict[(newTruncLemma3, trunc3)] = self.frequencyDict[newTruncLemma3]

                        try:
                            if max(decideDict.items(), key=operator.itemgetter(1))[1] != 0:
                                best = max(decideDict.items(), key=operator.itemgetter(1))[0]
                                resultList = [max(decideDict.items(), key=operator.itemgetter(1))[0]]
                                resultTruncLemma = resultList[0][1]
                            else:  # frequency of most frequent possibility == 0
                                resultTruncLemma = "*undecidable*"
                                self.undecDict["findPattern_8"] += 1
                        except:  # no dict because word smaller than 6 characters -> no segments
                            resultTruncLemma = "*undecidable*"

                        if resultTruncLemma == "*undecidable*":
                            for entries in self.taggedList:
                                if re.sub(r'[\\, ~, #, |, +]', "", trunc2.lower()) in entries[
                                    1].lower() and '+' in trunc2:
                                    self.checkSegment(trunc2, current + 4, 4)
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc2))
                                    break
                                elif re.sub(r'[\\, ~, #, |, +]', "", trunc3.lower()) in entries[
                                    1].lower() and '+' in trunc3:
                                    self.checkSegment(trunc3, current + 3, 3)
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc3))
                                    break
                            if '+' not in self.words[current].get('lemma'):
                                if '+' in trunc2 and '+' in trunc3:
                                    tr2 = re.sub(r'[\\, ~, |, +, #]', "", trunc2.split('+')[1])
                                    tr3 = re.sub(r'[\\, ~, |, +, #]', "", trunc3.split('+')[1])
                                    for entries in self.taggedList:
                                        if tr2.lower() == entries[1].lower() or re.sub(r'[\\, ~, |, +, #]', "",
                                                                                       trunc2.lower()) == entries[
                                            1].lower():
                                            trunc3 = ''
                                            break
                                        elif tr3.lower() == entries[1].lower() or re.sub(r'[\\, ~, |, +, #]', "",
                                                                                         trunc3.lower()) == entries[
                                            1].lower():
                                            trunc2 = ''
                                            break
                                if '+' in trunc2:
                                    if trunc2.split('+')[0].endswith('-'):
                                        trunc2 = trunc2.split('+')[0][:-1] + '+' + trunc2.split('+')[1]
                                    if trunc2.split('+')[1].startswith('-'):
                                        trunc2 = trunc2.split('+')[0] + '+' + trunc2.split('+')[1][1:]
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc2))
                                    if '#' not in self.words[current + 4].get('lemma'):
                                        tempInd = self.words[current + 4].get('lemma').rfind(
                                            re.sub(r'[\\, ~, |]', "", trunc2.split('+')[1][:3]))
                                        if tempInd != -1:
                                            if self.words[current + 4].get('lemma').endswith('-'):
                                                self.words[current + 4].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 4].get(
                                                                                                'lemma')[
                                                                                            :(tempInd - 1)] + '#' +
                                                                                            trunc2.split('+')[1]))
                                            else:
                                                self.words[current + 4].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 4].get(
                                                                                                'lemma')[
                                                                                            :tempInd] + '#' +
                                                                                            trunc2.split('+')[1]))
                                elif '+' in trunc3:
                                    if trunc3.split('+')[0].endswith('-'):
                                        trunc3 = trunc3.split('+')[0][:-1] + '+' + trunc3.split('+')[1]
                                    if trunc3.split('+')[1].startswith('-'):
                                        trunc3 = trunc3.split('+')[0] + '+' + trunc3.split('+')[1][1:]
                                    self.words[current].set('lemma', re.sub(r'[\\, ~, |]', "", trunc3))
                                    if '#' not in self.words[current + 3].get('lemma'):
                                        tempInd = self.words[current + 3].get('lemma').rfind(
                                            re.sub(r'[\\, ~, |]', "", trunc3.split('+')[1][:3]))
                                        if tempInd != -1:
                                            if self.words[current + 3].get('lemma').endswith('-'):
                                                self.words[current + 3].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 3].get(
                                                                                                'lemma')[
                                                                                            :(tempInd - 1)] + '#' +
                                                                                            trunc3.split('+')[1]))
                                            else:
                                                self.words[current + 3].set('lemma', re.sub(r'[\\, ~, |]', "",
                                                                                            self.words[current + 3].get(
                                                                                                'lemma')[
                                                                                            :tempInd] + '#' +
                                                                                            trunc3.split('+')[1]))
                        # self.generateList(trunc2, re.sub(r'[\\, ~, #, |, +]',"", trunc2))
                        # self.generateList(trunc3, re.sub(r'[\\, ~, #, |, +]',"", trunc3))
                        self.substituteLemma(resultTruncLemma, current)

                    # (9) trunc kon NN
                    elif self.words[current + 2].get('pos') == "NN" or self.words[current + 2].get('pos') == "NE":
                        truncWord = self.words[current].text.strip()
                        newCompoundLemma = self.substituteCompoundLemma(self.words[current + 2], current + 2)
                        segments = self.analyseCompoundLemma(newCompoundLemma, truncWord)
                        newTruncLemma = self.findTruncLemma(truncWord, segments, current + 2)
                        self.checkSegment(newTruncLemma, current + 2, 2)
                        self.substituteLemma(newTruncLemma, current)

                # (10) TRUNC APPRART NN/NE/ADJA (APPRART instead of KON)
                elif (word.get('pos') == "TRUNC" and word.text[-1] == "-" and self.words[current + 1].get(
                        'pos') == "APPRART") and (
                        self.words[current + 2].get('pos') == "NN" or self.words[current + 2].get('pos') == "NE" or
                        self.words[current + 2].get('pos') == "ADJA"):
                    truncWord = self.words[current].text.strip()
                    newCompoundLemma = self.substituteCompoundLemma(self.words[current + 2], current + 2)
                    segments = self.analyseCompoundLemma(newCompoundLemma, truncWord)
                    newTruncLemma = self.findTruncLemma(truncWord, segments, current + 2)
                    self.checkSegment(newTruncLemma, current + 2, 2)
                    self.substituteLemma(newTruncLemma, current)

                # (11) WORD KON -WORD" (word after next starts with "-")
                elif self.words[current].text[0] == "-" and len(self.words[current].text) > 3 and self.words[
                    current - 1].get('pos') == "KON" and self.words[current - 2].get('pos') != "TRUNC":
                    truncWord = self.words[current].text.strip()
                    newCompoundLemma = self.substituteCompoundLemma(self.words[current - 2], current - 2)
                    segments = self.rvsdAnalyseCompoundLemma(newCompoundLemma, truncWord)
                    newTruncLemma = self.rvsdFindTruncLemma(truncWord, segments, current - 2)
                    self.rvsdCheckSegment(newTruncLemma, current - 2, -2)
                    self.substituteLemma(newTruncLemma, current)

            current = current + 1

        return self.words

    # substitute compound with Gertwol-version of compound
    def substituteCompoundLemma(self, compoundLemma, index):
        if compoundLemma.text in self.gertwolDict:
            newCompoundLemma = self.gertwolDict[compoundLemma.text].strip()
        # self.words[index].set('lemma', newCompoundLemma)
        # if not in Gertwol & lemma not unknown -> take lemma
        elif self.words[index].get('lemma') != "unk":
            newCompoundLemma = self.words[index].get('lemma')
        # not in Gertwol & lemma unknown -> check if compound can be splitted; otherwise take word
        else:
            newCompoundLemma = compoundLemma.text
            temp = self.words[index - 2].get('lemma')
            border = 0
            beginning = ''
            ending = ''
            # Check if entire word is in taggedList.
            for entries in self.taggedList:
                p = re.compile(entries[1])
                if p.search(newCompoundLemma) != None and abs(len(newCompoundLemma) - len(entries[1])) < 3:
                    if '#' in entries[0]:
                        splitted = re.sub(r'[\\, ~, |]', "", entries[0], re.UNICODE).split('#')
                        beginning = splitted[0]
                        ending = entries[3][len(beginning):]
                    elif '+' in entries[0]:
                        splitted = re.sub(r'[\\, ~, |]', "", entries[0], re.UNICODE).split('+')
                        beginning = splitted[0]
                        ending = entries[3][len(beginning):]

            # Some specific checks in order to inrease correct splittings for unknown words (especially proper nouns).
            if self.words[index].get('lemma').endswith('tal'):
                beginning = self.words[index].get('lemma')[:-3]
                ending = 'tal'
            elif self.words[index].get('lemma').endswith('Tal'):
                beginning = self.words[index].get('lemma')[:-3]
                ending = 'Tal'
            elif newCompoundLemma.endswith('tal') or newCompoundLemma.endswith('tals') or newCompoundLemma.endswith(
                    'tales'):
                beginning = newCompoundLemma[:(newCompoundLemma.find('tal'))]
                ending = 'tal'
            elif newCompoundLemma.endswith('Tal') or newCompoundLemma.endswith('Tals') or newCompoundLemma.endswith(
                    'Tales'):
                beginning = newCompoundLemma[:(newCompoundLemma.find('Tal'))]
                ending = 'Tal'
            elif newCompoundLemma.endswith('thal') or newCompoundLemma.endswith('thals') or newCompoundLemma.endswith(
                    'thales') or newCompoundLemma.endswith('thale'):
                beginning = newCompoundLemma[:(newCompoundLemma.find('thal'))]
                ending = 'tal'
            elif newCompoundLemma.endswith('Thal') or newCompoundLemma.endswith('Thals') or newCompoundLemma.endswith(
                    'Thales') or newCompoundLemma.endswith('Thale'):
                beginning = newCompoundLemma[:(newCompoundLemma.find('Thal'))]
                ending = 'Tal'
            elif newCompoundLemma.endswith(u'thäler') or newCompoundLemma.endswith(u'thälern'):
                beginning = newCompoundLemma[:(newCompoundLemma.find(u'thäler'))]
                ending = 'tal'
            elif newCompoundLemma.endswith(u'täler') or newCompoundLemma.endswith(u'tälern'):
                beginning = newCompoundLemma[:(newCompoundLemma.find(u'täler'))]
                ending = 'tal'
            elif newCompoundLemma.endswith(u'Thäler') or newCompoundLemma.endswith(u'Thälern'):
                beginning = newCompoundLemma[:(newCompoundLemma.find(u'Thäler'))]
                ending = 'Tal'
            elif newCompoundLemma.endswith(u'Täler') or newCompoundLemma.endswith(u'Tälern'):
                beginning = newCompoundLemma[:(newCompoundLemma.find(u'Täler'))]
                ending = 'Tal'
            if newCompoundLemma.endswith('alp'):
                beginning = newCompoundLemma[:-3]
                ending = 'alp'
            if newCompoundLemma.endswith('pass'):
                beginning = newCompoundLemma[:-4]
                ending = 'pass'
            if newCompoundLemma.endswith('Pass'):
                beginning = newCompoundLemma[:-4]
                ending = 'Pass'
            if newCompoundLemma.endswith('gletscher'):
                beginning = newCompoundLemma[:-9]
                ending = 'gletscher'
            if newCompoundLemma.endswith('gebiet'):
                beginning = newCompoundLemma[:-6]
                ending = 'gebiet'
            if newCompoundLemma.endswith('Gebiet'):
                beginning = newCompoundLemma[:-6]
                ending = 'Gebiet'
            if newCompoundLemma.endswith('Massiv'):
                beginning = newCompoundLemma[:-6]
                ending = 'Massiv'

            # Check if there is an entry 'trunc+(ending of compound)' in taggedList.
            if self.words[index].get('pos') != 'TRUNC' and compoundLemma.text.find('-') == -1 and (
                    beginning == '' and ending == ''):
                if '-' in temp:
                    bef = re.sub('-', '', temp)
                elif '#' in temp:
                    bef = temp.split('#')[0]
                else:
                    bef = temp
                for entries in self.taggedList:
                    for border in range(1, len(newCompoundLemma)):
                        p = re.compile(entries[1])
                        if p.search(bef + newCompoundLemma[border:]) != None and abs(
                                len(bef + newCompoundLemma[border:]) - len(entries[1])) < 3:
                            beginning = newCompoundLemma[:border]
                            ending = entries[3][len(bef):]
                            break

                # If the tagger can't return any useful result, check possible splittings.
                if beginning == '' and ending == '':
                    for bord in range(1, len(newCompoundLemma)):
                        if (newCompoundLemma[bord:].lower() in self.gertwolDict or newCompoundLemma[
                                                                                   bord:].title() in self.gertwolDict) and (
                        len(newCompoundLemma[bord:])) > 2:
                            ending = newCompoundLemma[bord:]
                            beginning = newCompoundLemma[:bord]
                            for k, v in self.gertwolDict.items():
                                if ending in k:
                                    if '#' in v:
                                        temporary = v.split('#')
                                        if temporary[len(temporary) - 1][:2] == ending[:2] and abs(
                                                len(temporary[len(temporary) - 1]) - len(ending)) < 4:
                                            ending = temporary[len(temporary) - 1]
                                            ending = re.sub(r'[\\, ~, #, |]', "", ending, re.UNICODE)
                                            break
                            break
                    ending = ending[:(ending.find('&#10'))]

                # If no suffix is found, check if the prefix can be found in Gertwol.
                if len(ending) < 2:
                    for bord in range(3, len(newCompoundLemma)):
                        for k, v in self.gertwolDict.items():
                            begin = newCompoundLemma[:bord]
                            if begin in k:
                                if begin + '#' in v:
                                    tempo = v.split('#')
                                    if tempo[0][:2] == begin[:2]:
                                        ending = newCompoundLemma[len(begin):]
                                        beginning = newCompoundLemma[:(len(newCompoundLemma) - len(ending))]
                                        for entries in self.taggedList:
                                            p = re.compile(entries[1])
                                            if p.search(ending) != None and abs(len(ending) - len(entries[1])) < 3:
                                                ending = entries[3]
                                        break

                # Lemmatize endings.
                for entries in self.taggedList:
                    p = re.compile(entries[1])
                    if p.search(temp + ending) != None and abs(len(temp + ending) - len(entries[1])) < 3:
                        ending = entries[3][len(temp):]

            # Check tokens with hyphen(s).
            if compoundLemma.text.find('-') != -1 and (ending == '' and beginning == ''):
                t = newCompoundLemma.rfind('-')
                b = newCompoundLemma.find('-')
                without = re.sub('-', '', newCompoundLemma)
                for entries in self.taggedList:
                    p = re.compile(re.sub('-', '', entries[1]))
                    if p.search(without) != None and abs(len(without) - len(re.sub('-', '', entries[1]))) < 3:
                        if '#' in entries[0]:
                            splitted = entries[0].split('#')
                            beginning = splitted[0]
                            ending = entries[3][len(beginning):]
                if ((t + 1) < len(newCompoundLemma)) and ending == '' and beginning == '':
                    if (len(re.findall('-', newCompoundLemma)) > 1 and newCompoundLemma[t + 1].isupper()) or len(
                            re.findall('-', newCompoundLemma)) == 1:
                        beginning = newCompoundLemma[:t]
                        ending = newCompoundLemma[t:]
                        for entries in self.taggedList:
                            p = re.compile(entries[1])
                            if p.search(ending[1:]) != None and abs(len(ending[1:]) - len(entries[1])) < 3:
                                ending = '-' + entries[3]
                    elif len(re.findall('-', newCompoundLemma)) > 1 and newCompoundLemma[b + 1].isupper():
                        beginning = newCompoundLemma[:b]
                        ending = newCompoundLemma[b:]
                        for entries in self.taggedList:
                            p = re.compile(entries[1])
                            if p.search(ending[1:]) != None and abs(len(ending[1:]) - len(entries[1])) < 3:
                                ending = '-' + entries[3]

            # Return new lemma.
            if beginning != '' and ending != '' and len(beginning) > 0:
                if beginning[len(beginning) - 1] == '-':
                    beginning = beginning[:(len(beginning) - 1)]
                if ending[0] == '-':
                    ending = ending[1:]
                self.words[index].set('lemma',
                                      re.sub(r'[\\, ~, |]', "", beginning) + '#' + re.sub(r'[\\, ~, |]', "", ending))
                self.undecDict["splitCompound_A"] -= 1
            else:
                self.words[index].set('lemma', compoundLemma.text)

        # self.generateList(newCompoundLemma, self.words[index].text) # Generates a list of words which should be lemmatized by the TreeTagger.
        print(newCompoundLemma)
        return newCompoundLemma

    # Generates a list with all (unlemmatized) compounds and truncs.
    def generateList(self, word, wordTxt):
        tem = [word, wordTxt]
        if tem not in self.compList:
            self.compList.append(tem)

    # Prints all compounds and truncs to 'compoundList.txt'.
    def outputList(self):
        with open('compoundList.txt', 'a') as extList:
            for entry in self.compList:
                extList.write('\n' + entry[0].encode('utf-8') + '\t' + entry[1].encode('utf-8'))

    # NOUN: split compoundLemma at word boundary (#, -, |)
    def analyseCompoundLemma(self, newCompoundLemma, truncWord):
        if "#" in newCompoundLemma:
            segments = newCompoundLemma.split("#")
        elif "-" in newCompoundLemma[1:-1]:
            segments = []
            list = newCompoundLemma.split("-")
            for item in list[1:]:
                seg = "-" + item
                segments.append(seg)
        elif "|" in newCompoundLemma:
            segments = newCompoundLemma.split("|")
        else:  # according to Gertwol not splitable -> find word boundary without Gertwol
            segments = self.splitCompound(newCompoundLemma, truncWord)
        # segments = newCompoundLemma
        return segments

    # ADJECTIVE: split compoundLemma at word boundary (#, -, |, ~) --> only ADJ split at ~
    def adjAnalyseCompoundLemma(self, newCompoundLemma, truncWord):
        if "#" in newCompoundLemma:
            segments = newCompoundLemma.split("#")
        elif "-" in newCompoundLemma[1:-1]:
            segments = []
            list = newCompoundLemma.split("-")
            for item in list[1:]:
                seg = "-" + item
                segments.append(seg)
        elif "|" in newCompoundLemma:
            segments = newCompoundLemma.split("|")
        elif "~" in newCompoundLemma:
            segments = newCompoundLemma.split("~")
        else:  # according to Gertwol not splitable -> find word boundary without Gertwol
            segments = self.splitCompound(newCompoundLemma, truncWord)
        return segments

    # REVERSED (necessary if "-" in word - to keep it at correct place)
    # split compoundLemma at word boundary (#, -, |)
    def rvsdAnalyseCompoundLemma(self, newCompoundLemma, truncWord):
        if "#" in newCompoundLemma:
            segments = newCompoundLemma.split("#")
        elif "-" in newCompoundLemma[1:-1]:
            segments = []
            list = newCompoundLemma.split("-")
            for item in list[:-1]:
                seg = item + "-"
                segments.append(seg)
        elif "|" in newCompoundLemma:
            segments = newCompoundLemma.split("|")
        else:  # according to Gertwol not splitable -> find word boundary without Gertwol
            segments = self.rvsdSplitCompound(newCompoundLemma, truncWord)
        return segments

    # if Gertwol can't split the compound:
    # generate every possible split (but at least 3 characters in one "syllable")
    # and take the one which appears most frequent in the corpus
    def splitCompound(self, compound, trunc):
        decideDict = defaultdict(int)
        border = 3
        rest = compound[border:]
        best = ""
        while len(rest) >= 3:
            probableWord = trunc[:-1] + rest
            decideDict[tuple([probableWord, rest])] = self.frequencyDict[probableWord]
            border += 1
            rest = compound[border:]

        try:
            if max(decideDict.items(), key=operator.itemgetter(1))[1] != 0:
                best = max(decideDict.items(), key=operator.itemgetter(1))[0]
                segments = [max(decideDict.items(), key=operator.itemgetter(1))[0][1]]
            else:  # frequency of most frequent possibility == 0
                segments = []
                self.undecDict["splitCompound_A"] += 1

        except:  # no dict because word smaller than 6 characters -> no segments
            segments = []
            self.undecDict["splitCompound_B"] += 1

        return segments

    # REVERSED (necessary in order to attach missing word-part at correct place)
    # if Gertwol can't split the compound:
    # generate every possible split (but at least 3 characters in one "syllable")
    # and take the one which appears most frequent in the corpus
    def rvsdSplitCompound(self, compound, trunc):
        decideDict = defaultdict(int)
        border = len(compound) - 3
        rest = compound[:border]
        best = ""
        while len(rest) >= 3:
            probableWord = rest + trunc[1:]
            decideDict[tuple([probableWord, rest])] = self.frequencyDict[probableWord]
            border -= 1
            rest = compound[:border]

        try:
            if max(decideDict.items(), key=operator.itemgetter(1))[1] != 0:
                best = max(decideDict.items(), key=operator.itemgetter(1))[0]
                segments = [max(decideDict.items(), key=operator.itemgetter(1))[0][1]]
            else:  # frequency of most frequent possibility == 0
                segments = []
                self.undecDict["rvsdSplitCompound_A"] += 1
        except:  # no dict because word smaller than 6 characters -> no segments
            segments = []
            self.undecDict["rvsdSplitCompound_B"] += 1

        return segments

    # for TRUNC KON ADJA NN - pattern --> simply merge TRUNC & NN into 1 word
    def mergeTruncNn(self, truncWord, compound, index):
        newTruncWord = truncWord + "+" + compound
        # self.generateList(newTruncWord, re.sub(r'[\\, ~, #, |, +]',"", newTruncWord))
        return newTruncWord

    # find correct TRUNC-Lemma
    def findTruncLemma(self, truncWord, segments, index):
        if len(segments) == 0:  # no segments -> no result
            newTruncWord = truncWord
        elif len(segments) > 0 and len(segments) <= 2:  # only one possibility
            newTruncWord = truncWord[:-1] + "+" + segments[-1]
        else:  # several possibilities to form new lemma -> find result through word-frequency in corpus
            trunc1 = truncWord[:-1] + "+" + segments[-1]
            trunc2 = truncWord[:-1] + "+" + segments[-2] + segments[-1]
            newTruncWord_1 = re.sub(r'[\\, ~, #, |]', "", trunc1)
            newTruncWord_2 = re.sub(r'[\\, ~, #, |]', "", trunc2)

            if self.frequencyDict[newTruncWord_1] > self.frequencyDict[newTruncWord_2]:
                newTruncWord = trunc1
            elif self.frequencyDict[newTruncWord_1] < self.frequencyDict[newTruncWord_2]:
                newTruncWord = trunc2
            else:  # not decidable by word-frequency
                newTruncWord = "*undecidable*"
                self.undecDict["findTruncLemma"] += 1

        if '#' in self.words[index].get('lemma'):
            suff = self.words[index].get('lemma').split('#')
            newTruncWord = re.sub(r'[\\, ~, |]', "", truncWord[:-1]) + '+' + suff[1]

        if self.words[index].get('lemma').endswith('tal'):
            newTruncWord = re.sub(r'[\\, ~, |]', "", truncWord[:-1]) + '+' + 'tal'

        if self.words[index].get('lemma').endswith('Tal'):
            newTruncWord = re.sub(r'[\\, ~, |]', "", truncWord[:-1]) + '+' + 'Tal'

        # self.generateList(newTruncWord, re.sub(r'[\\, ~, #, |, +]',"", newTruncWord))
        if '+' in newTruncWord:
            if newTruncWord.split('+')[0].endswith('-'):
                newTruncWord = newTruncWord.split('+')[0][:-1] + '+' + newTruncWord.split('+')[1]
            if newTruncWord.split('+')[1].startswith('-'):
                newTruncWord = newTruncWord.split('+')[0] + '+' + newTruncWord.split('+')[1][1:]
        return newTruncWord

    # Check that all segmentation-symbols are inserted (so long as it is decidable).
    def checkSegment(self, word, index, z):
        beg = ''
        endi = ''
        wor = re.sub(r'[\\, ~, |, +, #]', "", word)
        if '#' not in self.words[index].get('lemma'):
            if self.words[index].get('lemma').endswith('tal'):
                beg = self.words[index].get('lemma')[:-3]
                endi = 'tal'
            elif self.words[index].get('lemma').endswith('Tal'):
                beg = self.words[index].get('lemma')[:-3]
                endi = 'Tal'
            elif wor.endswith('tal') or wor.endswith('tals') or wor.endswith('tales'):
                beg = self.words[index].text[:(wor.find('tal'))]
                endi = 'tal'
            elif wor.endswith('Tal') or wor.endswith('Tals') or wor.endswith('Tales'):
                beg = self.words[index].text[:(wor.find('Tal'))]
                endi = 'Tal'
            elif wor.endswith('thal') or wor.endswith('thals') or wor.endswith('thales') or wor.endswith('thale'):
                beg = self.words[index].text[:(wor.find('thal'))]
                endi = 'tal'
            elif wor.endswith('Thal') or wor.endswith('Thals') or wor.endswith('Thales') or wor.endswith('Thale'):
                beg = self.words[index].text[:(wor.find('Thal'))]
                endi = 'Tal'
            elif wor.endswith(u'thäler') or wor.endswith(u'thälern'):
                beg = self.words[index].text[:(wor.find(u'thäler'))]
                endi = 'tal'
            elif wor.endswith(u'täler') or wor.endswith(u'tälern'):
                beg = self.words[index].text[:(wor.find(u'täler'))]
                endi = 'tal'
            elif wor.endswith(u'Thäler') or wor.endswith(u'Thälern'):
                beg = self.words[index].text[:(wor.find(u'Thäler'))]
                endi = 'Tal'
            elif wor.endswith(u'Täler') or wor.endswith(u'Tälern'):
                beg = self.words[index].text[:(wor.find(u'Täler'))]
                endi = 'Tal'
        if word != "*undecidable*" and beg == '' and endi == '':
            # self.words[index-z].set('lemma', re.sub(r'[\\, ~, |]',"", word))
            if '+' in word:
                endi = re.sub(r'[\\, ~, |]', "", word.split('+')[1])
                tempInd = self.words[index].text.lower().rfind(endi[:3].lower())
                if tempInd != -1:
                    beg = self.words[index].text[:tempInd]
                else:
                    if 'a' in endi[:3]:
                        k = endi[:3].find('a')
                        n = endi[:3][:k] + u'ä' + endi[:3][(k + 1):]
                        if n in self.words[index].text:
                            ind = self.words[index].text.rfind(n)
                            beg = self.words[index].text[:ind]
                    if 'o' in endi[:3]:
                        k = endi[:3].find('o')
                        n = endi[:3][:k] + u'ö' + endi[:3][(k + 1):]
                        if n in self.words[index].text:
                            ind = self.words[index].text.rfind(n)
                            beg = self.words[index].text[:ind]
                    elif 'u' in endi[:3]:
                        k = endi[:3].find('u')
                        n = endi[:3][:k] + u'ü' + endi[:3][(k + 1):]
                        if n in self.words[index].text:
                            ind = self.words[index].text.rfind(n)
                            beg = self.words[index].text[:ind]
        if beg != '' and endi != '' and len(beg) > 0:
            if beg[len(beg) - 1] == '-':
                beg = beg[:(len(beg) - 1)]
            if endi[0] == '-':
                endi = endi[1:]
            if self.words[index].get('pos') == 'VVPP':
                if beg.endswith('ge'):
                    beg = beg[:-2]
            if self.words[index].get('pos') == 'VVIZU':
                if beg.endswith('zu'):
                    beg = beg[:-2]
            self.words[index].set('lemma', re.sub(r'[\\, ~, |]', "", beg) + '#' + endi)

    def rvsdCheckSegment(self, word, index, z):
        if word != "*undecidable*":
            # self.words[index+z].set('lemma', re.sub(r'[\\, ~, |]',"", word))
            if '+' in word:
                beg = re.sub(r'[\\, ~, |]', "", word.split('+')[0])
                endi = ''
                tempWord = re.sub(r'[\\, ~, |, +, #]', "", word)
                tempEndi = self.words[index].text[len(beg):]
                # Lemmatize the second word-segment.
                for entries in self.taggedList:
                    if entries[1].lower() == self.words[index].text.lower():
                        if len(entries[3][len(beg):]) > 1:
                            endi = entries[3][len(beg):].lower()
                            break
                if endi == '':
                    for k, v in self.gertwolDict.items():
                        if tempEndi.lower() == k.lower():
                            endi = re.sub(r'[\\, ~, |, +, #]', "", v.lower())
                            endi = endi[:(endi.find('&#10'))]
                            break
                if endi == '':
                    for entries in self.taggedList:
                        if entries[1].lower() == tempEndi.lower():
                            if len(entries[3]) > 1:
                                endi = entries[3].lower()
                                break
                # No lemma found. -> Take wordform in text.
                if endi == '':
                    endi = tempEndi
                # Return results.
                if beg != '' and endi != '' and len(beg) > 0:
                    if beg[len(beg) - 1] == '-':
                        beg = beg[:(len(beg) - 1)]
                    if endi[0] == '-':
                        endi = endi[1:]
                    self.words[index].set('lemma', re.sub(r'[\\, ~, |]', "", beg) + '#' + endi)

    # REVERSED (first word = compound, second word = "TRUNC")
    # find correct TRUNC-Lemma
    def rvsdFindTruncLemma(self, truncWord, segments, index):
        if len(segments) == 0:  # no segments -> no result
            newTruncWord = truncWord
        elif len(segments) <= 2:  # only one possibility
            newTruncWord = segments[0] + "+" + truncWord[1:]
        else:  # several possibilities to form new lemma -> find result through word-frequency in corpus
            trunc1 = segments[0] + "+" + truncWord[1:]
            trunc2 = segments[0] + "+" + segments[1] + truncWord[1:]
            newTruncWord_1 = re.sub(r'[\\, ~, #, |]', "", trunc1)
            newTruncWord_2 = re.sub(r'[\\, ~, #, |]', "", trunc2)

            if self.frequencyDict[newTruncWord_1] > self.frequencyDict[newTruncWord_2]:
                newTruncWord = trunc1
            elif self.frequencyDict[newTruncWord_1] < self.frequencyDict[newTruncWord_2]:
                newTruncWord = trunc2
            else:  # not decidable by word-frequency
                newTruncWord = "*undecidable*"
                self.undecDict["rvsdFindTruncLemma"] += 1

        if '+' in newTruncWord:
            for k, v in self.gertwolDict.items():
                if newTruncWord.split('+')[1].lower() == k.lower():
                    newTruncWord = newTruncWord.split('+')[0] + '+' + re.sub(r'[\\, ~, |, +, #]', "", v.lower())
                    newTruncWord = newTruncWord[:(newTruncWord.find('&#10;'))]
                    break

        if '+' in newTruncWord and len(newTruncWord.split('+')[0]) > 0:
            s1 = newTruncWord.split('+')[0]
            s2 = newTruncWord.split('+')[1]
            if s1[len(s1) - 1] == '-':
                s1 = s1[:(len(s1) - 1)]
            if s2[0] == '-':
                s2 = s2[1:]
            newTruncWord = s1 + '+' + s2

        if '#' in self.words[index].get('lemma'):
            pref = self.words[index].get('lemma').split('#')
            for entries in self.taggedList:
                p = re.compile(entries[1])
                if p.search(truncWord[1:]) != None and abs(len(truncWord[1:]) - len(entries[1])) < 3:
                    newTruncWord = pref[0] + '+' + entries[3]
                else:
                    if truncWord[1:] == u'Gneiße':
                        trWord = 'Gneis'
                    else:
                        trWord = truncWord[1:]
                    newTruncWord = pref[0] + '+' + trWord

        # self.generateList(newTruncWord, re.sub(r'[\\, ~, #, |, +]',"", newTruncWord))
        return newTruncWord

    # substitute original TRUNC-lemma with resulting lemma
    def substituteLemma(self, result, index):
        if result != "*undecidable*":
            newResult = re.sub(r'[\\, ~, |]', "", result)
            self.words[index].set('lemma', newResult)
            self.count += 1

    # write output to new xml-file
    def output(self, outfilename):
        with open(outfilename, "wb") as file:
            self.document.write(file, encoding="utf-8", xml_declaration=True)

    # give information about the frequency & method name where a case couldn't be decided
    # in those cases nothing happened (input = output)
    def infoUndecidables(self, infilename):
        print("{0} lemmas have been changed in {1}.".format(self.count, infilename))
        if (len(self.undecDict) == 0 or (
                self.undecDict.get("splitCompound_A") == 0 and self.undecDict.get("splitCompound_B") == None)):
            print("There were no undecidable cases.")
        else:
            print("Function name & frequency of undecidable cases which have not been changed:")
            for k, v in self.undecDict.items():
                print("{0: <20} {1}".format(k, v))


def exchangeLemmas(originalFile, outputFile, gertwolList, wordFreqs):
    fixLemma = EllipticCompound(originalFile, outputFile, gertwolList, wordFreqs)
    fixLemma.findPattern()
    fixLemma.output(outputFile)
    fixLemma.infoUndecidables(originalFile)


# fixLemma.outputList()


"""if __name__ == "__main__":
	exchangeLemmas()"""
