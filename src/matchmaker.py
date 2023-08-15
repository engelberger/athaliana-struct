from chimerax.core.commands import run
import argparse


# parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument('--dir', type=str, help='The directory where the models are located')
parser.add_argument('--file_1', type=str, help='The name of the first file')
parser.add_argument('--file_2', type=str, help='The name of the second file')
# An argument for output html file
parser.add_argument('--output', type=str, help='The name of the output html file')
args = parser.parse_args()


run(session, f" open {args.file_1} {args.file_2}")
run(session, 'mmaker #2 to #1 showAlignment true verbose true')
run(session, f'log save {args.dir}/{args.output}')
run(session, 'exit')
# Create a list of the models in the session
models_names_list = []

for m in session.models:
    models_names_list.append(m.name)
    

#count = 1
## Create for loop from that loops over the elements in the list of models in the session
#for m in session.models:
#    # run show models
#    run(session, 'show models')
#    # Create a empty dataframe to store the interface residues
#    interface_residues_df = pd.DataFrame()
#    #interface_residues = run(session, f'select #1/B:228,301').residues
#    run(session, f'interfaces select #{count}/{chain_id_0} contacting #{count}/{chain_id_1} bothSides false')
#    interface_residues = run(session, 'select sel').residues
#    print(interface_residues)
#    # Create a list with the names of the residues in the interface
#    interface_residues_names_list = []
#    for r in interface_residues.names:
#        interface_residues_names_list.append(r)
#    # Create a list with the numbers of the residues in the interface
#    interface_residues_numbers_list = []
#    for n in interface_residues.numbers:
#        interface_residues_numbers_list.append(n)
#    # Create a list with the chain ids of the residues in the interface
#    interface_residues_chain_ids_list = []
#    for c in interface_residues.chain_ids:
#        interface_residues_chain_ids_list.append(c)
#    # Print the list of the names of the residues in the interface and the length of the list
#    print(interface_residues_names_list)
#    print(len(interface_residues_names_list))
#    # Print the list of the numbers of the residues in the interface and the length of the list
#    print(interface_residues_numbers_list)
#    print(len(interface_residues_numbers_list))
#    # Print the list of the chain ids of the residues in the interface and the length of the list
#    print(interface_residues_chain_ids_list)
#    print(len(interface_residues_chain_ids_list))
#    
#    
#    # Add the lists to a dataframe as columns
#    interface_residues_df['residue_names'] = interface_residues_names_list
#    interface_residues_df['residue_numbers'] = interface_residues_numbers_list
#    interface_residues_df['residue_chain_ids'] = interface_residues_chain_ids_list
#    model_name = m.name
#    # Split the model name by _unrelaxed and preserve the first part of the name
#    model_name_split = model_name.split('_unrelaxed')
#    # Split the name of the model by _whole and preserve the first part of the name
#    model_name_split_whole = model_name.split('_whole_')[1]
#    # Int the model_name_split_whole replace the string '_unrelaxed_rank_1' with empty string ''
#    model_name_split_whole_empty = model_name_split_whole.replace('_unrelaxed_rank_1', '')
#
#    interface_residues_df.to_csv(f'{dir}/{model_name_split_whole_empty}_interface_residues.csv',index=False)
#    # add 1 to the count
#    count += 1
#    # Run hide models
#    run(session, f'hide ~#{count} models')
#    run(session, f'save {dir}/figures/{model_name_split_whole_empty}.png transparentBackground true width 3000 height 3000 supersample 3')
## Save the dataframe to a csv file in the dir
##