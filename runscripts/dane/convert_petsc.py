import scipy.io, PetscBinaryIO
import sys

def convert(file_in, file_out):
    A = scipy.io.mmread(file_in)
    A = A.tocsr()
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(file_out,'w'), A)

if __name__=='__main__':
    args = sys.argv[1:]
    print(args)
    convert(args[0], args[1])
