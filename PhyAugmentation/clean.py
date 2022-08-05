import os
import pickle

def remove():
     ### Clean the directory
    skip_list = ['build.py', 'hdock', 'tcr_unknown.pickle', 'pep_database.pickle', 'main.py',
                 'clean.py', 'results.txt', 'peptide_build.py', 'hdock_example.py', 'README.md']
    for item in os.listdir():
        if item in skip_list:
            pass
        else:
            os.remove(item)

if __name__ == '__main__':
    # peps = pickle.load( open('pep_database.pickle', 'rb' ) )
    # for pep in peps:
    #     print(pep)
    # print(len(peps))
    remove()
