

""" Prints the top x percent of VP indexes from a file containing all the indexes from a model """
import sys
from collections import defaultdict

def main():
    file = sys.argv[1]
    content = open(file).readlines()
    content = [line.strip() for line in content]
    content.sort()
    d = defaultdict(int)
    for ind in content:
        d[ind] += 1
    l = [(val,key) for key,val in d.items()]
    l.sort()
    l.reverse()
    percent = .90
    total = float(len(content))
    ss = 0.0
    for x in l:
        ss += x[0]/total
        #print(x,ss,x[0]/total)
        for i in range(x[0]):
            print(x[1])
        if(ss > percent): break
        


if __name__ == "__main__":
    main()
