import pandas as pd
import synapseclient

syn = synapseclient.Synapse()
syn.login()


def load_table(table_id, verbose=False):
    """"""
    print(f"Loading {table_id}")
    table = syn.tableQuery(f"select * from {table_id}").asDataFrame()
    if verbose:
        print(table.head())
        # [Molecular, DataType, genelists, MSE, numFeatures, numSamples, corVal, compound, method, meanCor]
        print(table.genelists)
    return table


def load_file(file_id, delimiter='\t'):
    table = pd.read_csv(syn.get(file_id).path, delimiter=delimiter)
    return table


def load_excel(file_id):
    table = pd.read_excel(
        syn.get(file_id).path,
        engine='openpyxl'
    )
    return table


if __name__ == '__main__':
    t_id = 'syn26469964'
    # load_table(table_id=t_id)
    print(syn.get('syn26718014'))
    load_file(file_id='syn26718014')

    # t = load_excel('syn26532699')
