import pandas as pd


def get_source_data():
    df = pd.read_csv('mapping.csv')
    col1 = 'source_property'
    col2 = 'nde_property'
    result = dict(zip(df[col1], df[col2]))
    result = {k: v for k, v in result.items() if not pd.isna(k)
              and not pd.isna(v)}
    return result
