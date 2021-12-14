import scipy.sparse
from scipy.io import savemat, loadmat
import numpy as np

def genSDPdata(m, nlist):
    
    ncones = len(nlist)
    SDPdata = {}
    SDPOriginal = {}
    for i in range(ncones):
        print("Generating cone {0}".format(i + 1))
        A, pA, mtypeinBlock = genBlockData(m + 1, nlist[i])
        print("Cone statistic:"
              " {0} sparse {1} dense {2} "
              "rank one {3} zero".format(mtypeinBlock[0],
                                         mtypeinBlock[1],
                                         mtypeinBlock[2],
                                         mtypeinBlock[3]))
        SDPdata[i] = pA
        SDPOriginal[i] = A
        
    return SDPdata, SDPOriginal

def genBlockData(m, dim):
    """
    Generate a block of data, containing m + 1 matrices of dimension dim
    
    """
    blockoriginal = {}
    blockpacked = []
    mtypeinblock = [0, 0, 0, 0]
    
    for i in range(m):
        A, mattype = genSomeMat(dim)
        Apacked = np.asarray(A.todense().T[np.triu_indices(dim, 0)])[0]
        blockoriginal[i] = A
        blockpacked.append(Apacked)
        mtypeinblock[mattype] += 1
        
    return blockoriginal, scipy.sparse.csr_matrix(blockpacked), mtypeinblock

def genSomeMat(dim, mattype=None):
    """
    Randomly generate a matrix from 
    
    - Sparse
    - Dense 
    - Rank 1
    - Zero
    
    """
    if mattype is None:
        mattype = np.random.randint(4)
    
    if mattype == 0:
        # Sparse matrix
        A = scipy.sparse.random(dim, dim, density=0.6, format="csr")
        return A + A.T, mattype
    elif mattype == 1:
        A = scipy.sparse.random(dim, dim, density=1.0, format="csr")
        return A.dot(A.T), mattype
    elif mattype == 2:
        x = scipy.sparse.random(dim, 1, 0.8, format="csr")
        return x.dot(x.T), mattype
    else:
        return scipy.sparse.random(dim, dim, density=0.0, format="csr"), mattype

def writeSDPData(SDPdata, coneDims):

    m = SDPdata[0].shape[0]
    b = np.random.randn(m)
    b_str = ", ".join([str(b[i]) for i in np.arange(m)])
    
    with open("SDPdata.h", "w") as f:
        f.write("int m = {0};\n".format(m - 1))
        f.write("int    ncones = {0};\n".format(len(SDPdata)))
        f.write("double b[{0:d}] = {{{1:s}}};\n".format(b.size, b_str))

        for key, val in SDPdata.items():

            print("Writing SDP cone {0}".format(key + 1))
            A = val
            Ap = A.indptr
            Ai = A.indices
            Ax = A.data
            Ap_str = ", ".join([str(Ap[i]) for i in np.arange(Ap.size)])
            Ai_str = ", ".join([str(Ai[i]) for i in np.arange(Ai.size)])
            Ax_str = ", ".join([str(Ax[i]) for i in np.arange(Ax.size)])
            f.write("int  cdim{0} = {1};\n".format(key + 1, coneDims[key]))
            f.write("int    Ap_{0}[{1:d}] = {{{2:s}}};\n".format(key + 1, Ap.size, Ap_str))
            f.write("int    Ai_{0}[{1:d}] = {{{2:s}}};\n".format(key + 1, Ai.size, Ai_str))
            f.write("double Ax_{0}[{1:d}] = {{{2:s}}};\n".format(key + 1, Ax.size, Ax_str))
        
    f.close()

        
if __name__ == "__main__":

    np.random.seed(24)
    # Number of constraints
    m = 10
    # Dimension of SDP cones
    nlist = [3, 3, 4]
    SDPdata, SDPOriginal = genSDPdata(m, nlist)
    writeSDPData(SDPdata, nlist)

    savemat("SDPdata.mat", 
            {"SDPdata{0}".format(key): 
            {"A_{0}".format(i): SDPOriginal[key][i] for i in range(m + 1)} 
            for key in SDPOriginal.keys()})

    print("Data is already written.")
