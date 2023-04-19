## v1.15.0 (2023-04-19)

### Fix
- Add support for random_seed in sample_enrichment() function

## v1.12.0 (2022-05-30)

### Feat

- Add the plot_term_pcoa() function for plotting dbBact terms based PCA

## v1.11.0 (2021-08-25)

### Feat

- Added sample_to_many_enrich() to identify features enriched in one set of samples (even one sample) compared to another set of samples

### Fix

- **requirements.txt**: add requests to requirements so importing dbbact-calour will work

## v1.10.1 (2021-08-12)

### Fix

- fix deepcopy in enrichment()

## v1.10.0 (2021-08-11)

### Feat

- Add support for "parentterm" term_type in enrichment()

## v1.9.0 (2021-07-25)

### Feat

- change to ontology term_id based annotations in add annotation gui

## v1.8.0 (2021-06-23)

### Feat

- Add annotationid to annotation enrichemnt analysis

### Fix

- fix load message for adding dbbact terms
- remove debug prints

## v1.7.0 (2021-05-28)

### Feat

- Fix enrichment so replicates results across runs with same random_seed

### Fix

- change print to logger.debug

## v1.6.0 (2021-05-26)

### Feat

- update label in venn diagram to show corrent number if clipped at max_size

## v1.5.1 (2021-05-20)

### Fix

- remove old code for background enrichment

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
## v1.12.0 (2022-05-30)

### Feat

- Add the plot_term_pcoa() function for plotting dbBact terms based PCA

## v1.11.0 (2021-08-25)

### Feat

- Added sample_to_many_enrich() to identify features enriched in one set of samples (even one sample) compared to another set of samples

### Fix

- **requirements.txt**: add requests to requirements so importing dbbact-calour will work

## v1.10.1 (2021-08-12)

### Fix

- fix deepcopy in enrichment()

## v1.10.0 (2021-08-11)

### Feat

- Add support for "parentterm" term_type in enrichment()

## v1.9.0 (2021-07-25)

### Feat

- change to ontology term_id based annotations in add annotation gui

## v1.8.0 (2021-06-23)

### Feat

- Add annotationid to annotation enrichemnt analysis

### Fix

- fix load message for adding dbbact terms
- remove debug prints

## v1.7.0 (2021-05-28)

### Feat

- Fix enrichment so replicates results across runs with same random_seed

### Fix

- change print to logger.debug

## v1.6.0 (2021-05-26)

### Feat

- update label in venn diagram to show corrent number if clipped at max_size

## v1.5.1 (2021-05-20)

### Fix

- remove old code for background enrichment

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
