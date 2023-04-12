from IPython.display import display, display_markdown, Markdown, HTML
from evaluation import error_analysis, get_scores, encode_decode

from difflib import HtmlDiff
from bs4 import BeautifulSoup
import pandas as pd

def show_errors(df):
    for i, r in df.reset_index().iterrows():
        display(Markdown(f"{i};{r.file};{r.sentence_id}"))
        display(Markdown(f"__Input:__"))
        display(Markdown(r.original))
        if r.error_type == "tp":
            display(Markdown("__Prediction (correct):__"))
            display(Markdown(r.pred))
        else:
            display(Markdown(f"__Error type:__ {r.error_type}"))
            html_diff =  HtmlDiff().make_file([r.ground_truth], [r.pred])

            soup = BeautifulSoup(html_diff, features="html.parser")
            for tab in soup.find_all('table')[1:]:
                tab.decompose()
            tab0 = soup.find("table")

            for tag in tab0.find_all("td"):
                if tag.has_attr('class'):
                    tag.decompose()
                else:
                    tag.name = "div"
            display(HTML(str(soup)))
            
def print_eval_row(errors, scores):
    col_order = ['exact_match', 'edit_distance_rel', 'gleu', 'tp', 'tn', 'fp', 'fn', 'insert', 'delete', 'replace', 'complex'] 
    errors_df = (errors.error_type.value_counts()).to_frame().T# / len(errors)).to_frame().T
    row = pd.concat([pd.DataFrame([scores]).rename(columns=lambda n: n.split('/')[1]), errors_df.reset_index()], axis=1)
    return row.round(3)[col_order]
    #return row
    
def calculate_errors(out, sample, tokenizer):
    gen_text = [o if type(o) == str else o['generated_text'] for o in out]
    errors = error_analysis(gen_text, encode_decode(sample.full_resolution, tokenizer), encode_decode(sample.raw_sentence, tokenizer))
    errors = pd.concat([errors, sample[['file', 'sentence_id']].reset_index()], axis=1)
    return errors