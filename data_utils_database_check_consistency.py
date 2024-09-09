#!/usr/bin/env python3
import sys
import os
import argparse

# Add the directory containing data_utils to sys.path
sys.path.insert(0, os.path.expanduser('~/lib'))


def main():
    parser = argparse.ArgumentParser(description='Check the format of the data_utils_database.txt file.')
    parser.add_argument('-f', '--database_file', type=str, required=True,
                        help='Path to the data_utils_database.txt file')
    args = parser.parse_args()
    d_handles_map={
        'magnuson2018_30348759':'CANCER_2018_30348759_MAGNUSON',
        'luoma2020_32603654':'CANCER_2020_32603654_LUOMA',
        'chen2024_38981439':'COAD_2024_38981439_CHEN'}

    # Create a dictionary with the list as keys and None as values

    s_data_handles=set()
    with open(args.database_file) as fp:
        for x in fp:
            data_handle, rest =x.rstrip('\n\r').split('.',1)
            key, value=rest.split('=',1)
            if data_handle in d_handles_map:
                data_handle_new=d_handles_map[data_handle]
                print('{}.{}={}'.format(data_handle_new, key, value))
                
if __name__ == "__main__":
    main()
