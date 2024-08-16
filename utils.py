# util for natural sorting
import re

def cmp_natural(text):
    '''
    compare function for sorting sample or gene names with (string+number )

    USAGE: alist.sort(key=cmp_natural) sorts in human order
    '''
    def atoi(text): 
        return int(text) if text.isdigit() else text
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]
