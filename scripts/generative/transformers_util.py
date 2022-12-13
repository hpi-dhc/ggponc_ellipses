from transformers import Seq2SeqTrainingArguments, AutoModelForSeq2SeqLM, DataCollatorForSeq2Seq, AutoTokenizer, Seq2SeqTrainer
from evaluation import Metrics

def get_tokenizer(config):
    tokenizer = AutoTokenizer.from_pretrained(config.model_name)
    if config.model_name == "facebook/mbart-large-50-many-to-many-mmt":
        tokenizer.src_lang = 'de_DE'
        tokenizer.tgt_lang = 'de_DE'
    if config.model_name == "facebook/m2m100_418M":
        tokenizer.src_lang = 'de'
        tokenizer.tgt_lang = 'de'
    return tokenizer

def get_training_args(config, report_to=None):
    return Seq2SeqTrainingArguments(
        output_dir='./results',
        num_train_epochs=config.num_epochs,
        per_device_train_batch_size=config.train_batch_size,
        per_device_eval_batch_size=config.eval_batch_size,
        logging_dir='./logs',
        logging_steps=100,
        warmup_steps=config.warmup_steps,
        load_best_model_at_end=True,
        predict_with_generate=True,
        learning_rate=config.learning_rate,
        generation_max_length=config.generation_max_length,
        evaluation_strategy="epoch",
        save_strategy="epoch", # epoch
        report_to=report_to,
        fp16=config.fp16,
        metric_for_best_model="exact_match",
        save_total_limit=1,
    )

def get_trainer(config, tokenizer, training_args, train_data, val_data):
    model = AutoModelForSeq2SeqLM.from_pretrained(config.model_name)

    metrics = Metrics(config.metrics, tokenizer)
    data_collator = DataCollatorForSeq2Seq(tokenizer=tokenizer, model=model)

    return Seq2SeqTrainer(
        model=model,
        args=training_args,
        train_dataset=train_data,
        eval_dataset=val_data,
        data_collator=data_collator,
        compute_metrics=metrics.compute_metrics
    )
    

