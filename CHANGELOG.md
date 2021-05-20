## v1.5.0 (2021-05-20)

### Feat

- add background enrichment support

## v1.4.0 (2020-12-27)

### Feat

- update term heatmap colormap

### Fix

- update term heatmap to show by default the annotation

## v1.3.0 (2020-12-27)

### Feat

- Add UI for term enrichment method

## v1.2.0 (2020-12-18)

### Feat

- add fscore option and set fscore as default for ca.add_terms_to_features()
- add sample_term_scores() method and update rank_enrichment() and sample_enrichment() to use it
- add set_log_level() to DBBact class

### Refactor

- **DBBact**: move calculation of recall/precision etc. out of wordcloud function to get_wordcloud_stats

## v1.1.1 (2020-12-02)

### Fix

- **DBBact**: fix error in rank_enrichment

## v1.1.0 (2020-12-02)

### Feat

- **DBBact**: add rank_enrichment() to identify terms correlated with some per-feature stat

## v1.0.4 (2020-12-02)

### Fix

- remove docs for smaller install

## v1.0.3 (2020-12-02)

### Fix

- fix DBBact import in __init__.py

## v1.0.2 (2020-12-01)

### Refactor

- add github workflows for automatic release and publish

## v1.0.1 (2020-12-01)

### Refactor

- Start semantic versioning
