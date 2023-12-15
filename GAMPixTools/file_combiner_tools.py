import os
import numpy as np
import h5py
import argparse

#combine processed h5 files coming from the dune_tools main function.
#combines all the files in directory into a single h5 file with the same keys and the data of all the files
#outputs into output_filename in directory
def combine_output_files(directory, output_filename):
    #find file paths
    path_list = os.listdir(directory)
    path_list = [directory + '/' + path for path in path_list if path[-3:] == '.h5']
    
    #read files
    h5file_list = [h5py.File(path, 'r') for path in path_list]
    dictionary_list = [{**h5file} for h5file in h5file_list]
    
    #combine into a single dictionary
    output_dict = {}
    for key in dictionary_list[0].keys():
        
        #check if value should be concatenated or if its a 'datatype' value that should be the same for all files (like pdgMask)
        if np.array(dictionary_list[0][key][()]).ndim == 0: #do not concatenate, assert that all files have the same value
            value = dictionary_list[0][key][()]
            for dictionary in dictionary_list:
                assert dictionary[key][()] == value, key + ' is not the same in all files.'
            output_dict[key] = value
            continue      
        else: #concatenate
            values = np.concatenate([dictionary[key][()] for dictionary in dictionary_list])
            output_dict[key] = values
            
    #write to file
    with h5py.File(directory +'/'+ output_filename + '.h5', 'w') as f:
        f.update(output_dict)


def main(args):
    input_directory = args.input_directory
    output_filename = args.output_filename
    
    print('combining files')
    combine_output_files(input_directory, output_filename)

        
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('input_directory', type = str,
                        help = 'input directory containing h5 files to combine')
    parser.add_argument('-o', '--output_filename', type = str, default = 'combined_output',
                        help = 'name for outputfile after combining h5 files. location will be same directory as files are in')
    args = parser.parse_args()
    
    main(args)
    
