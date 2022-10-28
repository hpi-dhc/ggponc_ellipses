from torch.utils.data import Dataset

class EllipsesDataset(Dataset):

    def __init__(self, sentences, references, tokenizer):
        self.tokenizer = tokenizer

        self.sentences = sentences.apply(lambda x: self.tokenizer(x, return_tensors='pt', return_length=True))
        self.references = references.apply(lambda x: self.tokenizer(text_target=x, return_tensors='pt'))

        if len(sentences) == len(references):
            self.n_samples = len(sentences)
        else:
            raise ValueError('Different number of samples and labels')

    def __getitem__(self, index):
        return dict(
            input_ids=self.sentences.iloc[index]["input_ids"].flatten(),
            #length=self.sentences.iloc[index]["length"].flatten(),
            attention_mask=self.sentences.iloc[index]["attention_mask"].flatten(),
            labels=self.references.iloc[index]["input_ids"].flatten()
        )

    def __len__(self):
        return self.n_samples

def load_dataset(ellipses_file, control_file):
    pass