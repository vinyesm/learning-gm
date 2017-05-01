import synapseclient
import synapseutils

syn = synapseclient.Synapse()
syn.login('vinyesm', 'M02p4sse?')
#files = synapseutils.syncFromSynapse(syn, 'syn2820442')

#files = synapseutils.syncFromSynapse(syn,'syn2820442','/home/marina/Marina/learning-gm/DREAM5/sub_challenge1')
#files = synapseutils.syncFromSynapse(syn,'syn2867578','/home/marina/Marina/learning-gm/DREAM5/sub_challenge2')
files = synapseutils.syncFromSynapse(syn, 'syn2787211', '/home/marina/Marina/learning-gm/DREAM5/Network-Inference')

#entity = syn.get('syn2820442')
#print(entity)