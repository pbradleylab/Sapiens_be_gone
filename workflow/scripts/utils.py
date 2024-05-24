def get_subsample_attributes(subsample, attribute, pep):
    subsample_rows = pep.subsample_table[pep.subsample_table["subsample"] == subsample]
    return pep.get_sample(subsample_rows["sample_name"].tolist()[0])[attribute]