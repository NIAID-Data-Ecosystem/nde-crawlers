from biothings.utils.dataload import tab2dict


def get_source_data(source_name):
    datafile = f'hub/dataload/sources/{source_name}/{source_name}.csv'
    test = f'{source_name}.csv'
    result = tab2dict(test, cols=[0, 1], key=0, sep=',')
    # convert values that are lists to strings
    for key in result:
        if isinstance(result[key], list):
            result[key] = str(result[key])
    return result


print(get_source_data('biocontainers'))
