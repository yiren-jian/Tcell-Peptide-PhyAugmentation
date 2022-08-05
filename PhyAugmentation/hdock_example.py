import subprocess
import pickle

def main():
    p = subprocess.Popen("./hdock TCR.pdb peptide.pdb > hdock.out", stdout=subprocess.PIPE, shell=True)
    p_status = p.wait()

    with open('hdock.out') as f:
        lines = f.readlines()
    energy = lines[-2].split(' ')[-1].replace('\n','')
    energys.append(float(energy))
    print(energys)

    with open('energys.pickle', 'wb') as handle:
        pickle.dump(energys, handle)

if __name__ == '__main__':
    main()
