import pandas as pd
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

def load_data(cnf_tsv, control_tsv=None, sample_frac = None):
    dataset = pd.read_csv(cnf_tsv, sep='\t')
    dataset['controls'] = False
        
    if control_tsv:
        df_controls = pd.read_csv(control_tsv, sep='\t')
        df_controls['controls'] = True
        df_controls['full_resolution'] = df_controls.raw_sentence
        
        dataset = pd.concat([dataset, df_controls])
        
    train_df = dataset[dataset.split == 'train']
    val_df = dataset[dataset.split == 'dev']
    test_df = dataset[dataset.split == 'test']
    if sample_frac:
        train_df = train_df.sample(frac=sample_frac)
        val_df = val_df.sample(frac=sample_frac)
        test_df = test_df.sample(frac=sample_frac)
    return train_df, val_df, test_df

def get_dataloader(train_df, val_df, test_df, tokenizer):
    train_data = EllipsesDataset(train_df.raw_sentence, train_df.full_resolution, tokenizer)
    val_data = EllipsesDataset(val_df.raw_sentence, val_df.full_resolution, tokenizer)
    test_data = EllipsesDataset(test_df.raw_sentence, test_df.full_resolution, tokenizer)
    return train_data, val_data, test_data