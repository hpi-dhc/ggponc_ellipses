from IPython.display import display, display_markdown, Markdown, HTML
from difflib import HtmlDiff
from bs4 import BeautifulSoup

def show_errors(df):
    for _, r in df.iterrows():
        display(Markdown(f"__Input ({r.file},{r.sentence_id}):__"))
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