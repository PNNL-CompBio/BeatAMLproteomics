import synapseclient
import synapseutils

syn = synapseclient.Synapse()
syn.login()
files = synapseutils.syncFromSynapse(
    syn, 'syn51400754 ',
    path='nmf_trials',
    downloadFile=True
)

# syn, 'syn51400754 ',
#  path='nmf_trials',