import pdb
import numpy as np

def db2(S, MAX_LAYER = -1):
    #forward discrete wavelet transform db2 (1 vanishing moment)
    #S: array full of data
    #T: number of elements evaluated

    N = len(S)
    T = N
    
    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N
        #make sure that array can be resized
        S.resize(p2, refcheck=False)
        print S
        T = p2
    
    layers = 0
    step = 2
    
    while step < 2*T:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break

        start = step/2
        
        #predict 1
        S[start::step] -= S[::step]

        #update 1
        S[::step] += S[start::step]/2
         
        step <<= 1
        layers += 1
    
    return S

def idb2(W, MAX_LAYER = -1):
    N = len(W)

    layers = 0
    step = N

    #magic
    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >=MAX_LAYER >= 0): step = 1 << MAX_LAYER
    
    while step > 1:
        start = step/2
        
        #update' 1
        W[::step] -= W[start::step]/2
        
        #predict' 1
        W[start::step] += W[::step]

        step >>= 1
        layers += 1

    return W


def db4(S, MAX_LAYER = -1):
    #forward discrete wavelet transform    
    #S: array full of data
    #T: number of elements evaluated

    #split array, first half are even elements and second half are odd element
    N = len(S)
    T = N
    
    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N        S.resize(p2, refcheck=False)
        print S
        T = p2
    

    s3 = np.sqrt(3)
    s2 = np.sqrt(2)
    
    layers = 0
    step = 2
    
    while step < 2*T:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break

        start = step/2
        
        #update
        S[::step]  += s3 * S[start::step]
        
        #predict
        S[(start+step) :: step] += - s3/4*S[step::step] - (s3-2)/4*S[:-step:step] 
        S[start] += - s3/4*S[0] - (s3-2)/4*S[-step]
        
        #update 2
        S[:-step:step]  = S[:-step:step]  - S[(start+step) :: step]
        S[-step] = S[-step] - S[start]
        
        #normalize
        S[::step]  = (s3-1)/s2 * S[::step]  
        S[start::step] = (s3+1)/s2*S[start::step]
        
        step <<= 1
        layers += 1
        
    return S
    
def idb4(W, MAX_LAYER = -1):
    #inverse discrete wavelet transform    
    #same as forward, except all additions and subtractions are flipped
    #W: array full of wavelet coefficients
    #MAX_LAYER: total level of coefficients in transform
    #  It is possible for a transform to be partial, such as data with size 2^16 to go 8 levels deep instead of 16
 
    N = len(W)
    s3 = np.sqrt(3)
    s2 = np.sqrt(2)

    layers = 0
    step = N
    
    #magic
    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >= MAX_LAYER >= 0): step = 1 << MAX_LAYER
    
    while step > 1:
        start = step/2

        #normalize'
        W[::step] = (s3+1)/s2 * W[::step]
        W[start::step] = (s3-1)/s2 * W[start::step]

        #update 2'
        W[:-step:step] = W[:-step:step] + W[(start+step) :: step]
        W[-step] = W[-step] + W[start]
        
        #predict'        
        W[(start+step) :: step] += + s3/4*W[step::step]  + (s3-2)/4*W[:-step:step]
        W[start] += + s3/4*W[0] + (s3-2)/4*W[-step]
        
        #update 1'
        W[::step] += -s3 * W[start::step]
        
        step >>= 1
            
    return W


def layer_to_interlace(old):
    T = len(old)
    new = np.zeros(T)
    step = 2
    while T > 1:
        half = T/2
        start = step/2
        new[start::step] = old[half:T]

        start <<= 1
        T >>= 1
    new[0] = old[0]
    return new

def interlace_to_layer(old):
    T = len(old)
    new = np.zeros(T)
    step = 2
    while T > 1:
        half = T/2
        start = step/2
        new[half:T] = old[start::step]

        step <<= 1
        T >>= 1
    new[0] = old[0]
    return new

def waveletprint(W, MAX_LAYER = -1):
    N = len(W)
    logn = (N).bit_length() - 1
    
    initstep = N
    if logn > MAX_LAYER >= 0:
        initstep  >>= (logn - MAX_LAYER)
        
    step = initstep
    scale = 0
    print initstep

    lowdone = False
    if step == 1:
        print "\nScale 2^0:", W
        
    while step > 1:
        #stuff
        if step==initstep and not lowdone:
            start = 0
            lowdone = True
            sl = W[::step]
        else:
            start = step >> 1
            sl = W[start::step]
            step >>= 1

        st = "\nScale 2^-" + str(scale) + ":"
        print st, sl
        scale += 1

if __name__ == "__main__":
    s = np.array([32.0, 10.0, 20.0, 38.0, 37.0, 28.0, 38.0, 34.0, 18.0, 24.0, 18.0, 9.0, 23.0, 24.0, 28.0, 34.0])
    #s = np.arange(8, dtype = np.float64)
    #s = np.array([9.0, 7.0, 3.0, 5.0])
    #s = np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0])
    #s = np.array([1.0, 3.0, -2.0, 1.5, -0.5, 2.0, 0.0, 1.0])
    ml = -1
    #print s
    ds = db2(np.copy(s), MAX_LAYER = ml)
    #print ds
    #print interlace_to_layer(ds)
    ids = idb2(np.copy(ds), MAX_LAYER = ml)
    #print ids

waveletprint(ds, MAX_LAYER = ml)
