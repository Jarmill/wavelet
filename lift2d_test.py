import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image
import pdb

def db2_2d(S, MAX_LAYER = -1):
    """forward discrete wavelet transform db2/haar (1 vanishing moment)"""
    #S: array full of data
    #T: number of elements evaluated

    N = min(S.shape)
    T = N
    S = S.astype(np.double)
    """
    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N
        #make sure that array can be resized
        S.resize(p2, refcheck=False)
        print S
        T = p2
    """
    print "haar"
    layers = 0
    step = 2
    vert = False
    
    while step < 2*T:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break
        start = step/2

        #predict 1
        S[start::step, :] -= S[::step, :]

        #update 1
        S[::step, :] += S[start::step, :]/2.0
        
        if vert:
            step <<= 1
            layers += 1
            vert = False
            S = S.T
        else:
            vert = True
            S = S.T

    return S

def idb2_2d(W, MAX_LAYER = -1):
    """inverse discrete wavelet transform db2/haar (1 vanishing moment)"""
    W = W.T
    N = W.shape[0]

    step = N
    vert = False
    
    #magic
    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >=MAX_LAYER >= 0): step = 1 << MAX_LAYER

    while step > 1:
        start = step/2
        
        #update' 1
        W[::step, :] -= W[start::step, :]/2.0

        #predict' 1
        W[start::step, :] += W[::step, :]

        if vert:
            step >>= 1
            vert = False
            W = W.T
        else:
            vert = True
            W = W.T


    return (W.T).astype(np.uint8)

def db4_2d(S, MAX_LAYER = -1):
    """forward discrete wavelet transform db4 2d (2 vanishing moments)"""
    #S: array full of data
    #T: number of elements evaluated

    #split array, first half are even elements and second half are odd element
    N = min(S.shape)
    T = N
    vert = False
    S = S.astype(np.double)
    """
    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N        S.resize(p2, refcheck=False)
        print S
        T = p2
    """

    s3 = np.sqrt(3)
    s2 = np.sqrt(2)

    layers = 0
    step = 2

    while step < 2*T:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break

        start = step/2

        #update
        S[::step, :]  += s3 * S[start::step, :]

        #predict
        S[(start+step) :: step, :] += - s3/4*S[step::step, :] - (s3-2)/4*S[:-step:step, :]
        S[start, :] += - s3/4*S[0, :] - (s3-2)/4*S[-step, :]

        #update 2
        S[:-step:step, :]  = S[:-step:step, :]  - S[(start+step) :: step, :]
        S[-step, :] = S[-step, :] - S[start, :]

        #normalize
        S[::step, :]  = (s3-1)/s2 * S[::step, :]
        S[start::step, :] = (s3+1)/s2*S[start::step, :]

        if vert:
            step <<= 1
            layers += 1
            vert = False
            S = S.T
        else:
            vert = True
            S = S.T

    return S

def idb4_2d(W, MAX_LAYER = -1):
    """inverse discrete wavelet transform db4 2d (2 vanishing moments)"""
    #same as forward, except all additions and subtractions are flipped
    #W: array full of wavelet coefficients
    #MAX_LAYER: total level of coefficients in transform
    #  It is possible for a transform to be partial, such as data with size 2^16 to go 8 levels deep instead of 16

    N = min(W.shape)
    W = W.T
    s3 = np.sqrt(3)
    s2 = np.sqrt(2)

    layers = 0
    step = N
    vert = False

    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >= MAX_LAYER >= 0): step = 1 << MAX_LAYER

    while step > 1:
        start = step/2

        #normalize'
        W[::step, :] = (s3+1)/s2 * W[::step, :]
        W[start::step, :] = (s3-1)/s2 * W[start::step, :]

        #update 2'
        W[:-step:step, :] = W[:-step:step] + W[(start+step) :: step, :]
        W[-step, :] = W[-step, :] + W[start, :]

        #predict'        
        W[(start+step) :: step, :] += + s3/4*W[step::step, :]  + (s3-2)/4*W[:-step:step, :]
        W[start, :] += + s3/4*W[0, :] + (s3-2)/4*W[-step, :]

        #update 1'
        W[::step, :] += -s3 * W[start::step, :]

        if vert:
            step >>= 1
            vert = False
            W = W.T
        else:
            vert = True
            W = W.T


    return (W.T).astype(np.uint8)

def cdf97_2d(S, MAX_LAYER = -1):
    """CDF9/7 29 forward transformation, used for JPEG2000"""

    N = min(S.shape)
    S = S.astype(np.double)
    T = N
    """
    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N        S.resize(p2, refcheck=False)
        print S
        T = p2
    """
    #coefficients
    alpha = -1.5861343420693648
    beta = -0.0529801185718856
    gamma = 0.8829110755411875
    delta = 0.4435068520511142
    kodd = 1/1.1496043988602418
    keven = 1.1496043988602418

    layers = 0
    step = 2
    vert = False

    while step < 2*T:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break

        start = step/2

        #update 1
        #u1(z) = alpha(1 + z^-1)
        #odd[n] += alpha * (even[n] + even[n+1])
        S[start:-start:step, :]  += alpha * (S[:-step:step, :] + S[step::step, :])
        S[-start, :] += alpha * (S[-step, :] + S[0, :])
        
        #predict 1
        #p1(z) = beta(1 + z)
        #even[n] += beta * (odd[n] + odd[n-1])

        S[step::step, :] += beta * (S[step+start::step, :] + S[start:-start:step, :])
        S[0, :] += beta * (S[start, :] + S[-start, :])
        
        #update 2
        #u2(z) = gamma*(1 + z^-1)
        #odd[n] += gamma * (even[n] + even[n+1]

        S[start:-start:step, :] += gamma * (S[step::step, :] + S[:-step:step, :])
        S[-start, :] += gamma * (S[-step, :] +  S[0, :])

        #predict 2
        #p2(z) = delta * (1 + z)
        #even[n] += delta * (odd[n] + odd[n-1])
        S[step::step, :] += delta * (S[start+step::step, :] + S[start:-start:step, :])
        S[0, :] += delta * (S[start, :] + S[-start, :])

        #normalize
        #even[n] *= keven
        #odd[n]  *= kodd
        S[::step, :]  *= keven
        S[start::step, :] *= kodd
        
        if vert:
            step <<= 1
            layers += 1
            vert = False
            S = S.T
        else:
            vert = True
            S = S.T

    return S

def icdf97_2d(W, MAX_LAYER = -1):
    """CDF9/7 2d inverse transformation, used for JPEG2000"""
    #coefficients
    alpha = -1.5861343420693648
    beta = -0.0529801185718856
    gamma = 0.8829110755411875
    delta = 0.4435068520511142
    kodd = 1/1.1496043988602418
    keven = 1.1496043988602418
    
    N = min(W.shape)
    W = W.T
    step = N
    vert = False


    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >= MAX_LAYER >= 0): step = 1 << MAX_LAYER
    while step > 1:
        start = step/2
        
        #normalize'
        #even[n] /= keven
        #odd[n]  /= kodd
        W[::step, :]  /= keven
        W[start::step, :] /= kodd


        #predict 2'
        #p2(z) = delta * (1 + z)
        #even[n] -= delta * (odd[n] + odd[n-1])
        W[step::step, :] -= delta * (W[start+step::step] + W[start:-start:step])
        W[0, :] -= delta * (W[start] + W[-start])

        #update 2'
        #u2(z) = gamma*(1 + z^-1)
        #odd[n] -= gamma * (even[n] + even[n+1]

        W[start:-start:step, :] -= gamma * (W[step::step, :] + W[:-step:step, :])
        W[-start, :] -= gamma * (W[-step, :] +  W[0, :])
        
        #predict 1'
        #p1(z) = beta(1 + z)
        #even[n] -= beta * (odd[n] + odd[n-1])

        W[step::step, :] -= beta * (W[step+start::step, :] + W[start:-start:step, :])
        W[0, :] -= beta * (W[start, :] + W[-start, :])
        
        #update 1
        #u1(z) = alpha(1 + z^-1)
        #odd[n] -= alpha * (even[n] + even[n+1])
        W[start:-start:step, :]  -= alpha * (W[:-step:step, :] + W[step::step, :])
        W[-start, :] -= alpha * (W[-step, :] + W[0, :])


        if vert:
            step >>= 1
            vert = False
            W = W.T
        else:
            vert = True
            W = W.T


    return (W.T).astype(np.uint8)

def legall53_2d(S, MAX_LAYER = -1):
    """Legall 5/3 forward transformation, used for JPEG2000"""
    """
    N = len(S)
    T = N

    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N        S.resize(p2, refcheck=False)
        T = p2
    """
    N = min(S.shape)
    T = N
    layers = 0
    step = 2
    vert = False
    S = S.astype(np.double)

    while step < 2*T:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break

        start = step/2

        #update 1
        #u1(z) = (1 + z^-1)/2
        S[start:-start:step, :] -= np.floor((S[step::step, :] + S[:-step:step, :])/2)
        S[-start, :] -= np.floor((S[-step, :] + S[0, :])/2)

        #predict 1
        #p1(z) = (1 + z)/4
        S[step::step, :] += np.floor((S[step+start::step, :] + S[start:-start:step, :])/4)
        S[0] += np.floor((S[start, :] + S[-start, :])/4)
        if vert:
            step <<= 1
            layers += 1
            vert = False
            S = S.T
        else:
            vert = True
            S = S.T

    return S

def ilegall53_2d(W, MAX_LAYER = -1):
    """Legall 5/3 inverse transformation, used for JPEG2000"""    

    N = min(W.shape)
    step = N
    vert = False

    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >= MAX_LAYER >= 0): step = 1 << MAX_LAYER

    W = W.T
    
    while step > 1:
        start = step/2

        #predict 1'
        #p1(z) = (1 + z)/4
        W[step::step, :] -= np.floor((W[step+start::step, :] + W[start:-start:step, :])/4)
        W[0, :] -= np.floor((W[start, :] + W[-start, :])/4)

        #update 1'
        #u1(z) = (1 + z^-1)/2
        W[start:-start:step, :] += np.floor((W[step::step, :] + W[:-step:step, :])/2)
        W[-start, :] += np.floor((W[-step, :] + W[0, :])/2)

        if vert:
            step >>= 1
            vert = False
            W = W.T
        else:
            vert = True
            W = W.T


    return (W.T).astype(np.uint8)

"""
def mipmap(S, MAX_LAYER = -1):
    N = S.shape[0]
    m = np.zeros(N, 3*N/2)
"""    

def wavelet_reorder_2d(W, MAX_LAYER = -1):
    """print non-interlaced wavelet coefficients"""
    c = np.zeros_like(W)
    N = min(W.shape)
    layers = 1
    step = 2
    start = step/2

    endh = W.shape[0]
    midh = endh/2

    endv = W.shape[1]
    midv = endv/2
    
    while step < 2*N:
        if (MAX_LAYER >= 1 and layers > MAX_LAYER): break        
        c[midh:endh, midv:endv] = W[start::step, start::step]   #HH
        c[:midh, midv:endv] = W[start::step, ::step]            #LH
        c[midh:endh, :midv] = W[::step, start::step]            #HL

        endh /= 2
        midh /= 2
        endv /= 2
        midv /= 2
        step *= 2
        start *=2
        layers += 1
    c[:endh, :endv] = W[::start, ::start]                       #LL

    
    return c

def dwt2d(S, MAX_LAYER = 3, wavelet = "legall53"):
    """wrapper for 2d forward discrete wavelet transform"""
    if wavelet == "legall53":
        return legall53_2d(S, MAX_LAYER)
    elif wavelet == "cdf97":
        return cdf97_2d(S, MAX_LAYER)
    elif wavelet == "db2" or wavelet == "haar":
        return db2_2d(S, MAX_LAYER)
    elif wavelet == "db4":
        return db4_2d(S, MAX_LAYER)
    else:
        return legall53_2d(S, MAX_LAYER)

def idwt2d(W, MAX_LAYER = 3, wavelet = "legall53"):
    """wrapper for 2d forward discrete wavelet transform"""
    if wavelet == "legall53":
        return ilegall53_2d(W, MAX_LAYER)
    elif wavelet == "cdf97":
        return icdf97_2d(W, MAX_LAYER)
    elif wavelet == "db2" or wavelet == "haar":
        return idb2_2d(W, MAX_LAYER)
    elif wavelet == "db4":
        return idb4_2d(W, MAX_LAYER)
    else:
        return icdf97_2d(W, MAX_LAYER)

def psnr(truth, estimate):
    #returns Peak Signal to Noise Ratio in decibels
    pdb.set_trace()
    N = truth.size
    mean = np.mean(truth)
    #num = sum([abs(truth[i] - mean)**2 for i in range(N)])
    num = np.sum(truth ** 2)
    denom = np.sum((truth-estimate)**2)
    snr = 10 * np.log10(num/denom)

    return snr

def imtest(image = "dock.jpg", wavelet = "legall53", ml = 2):
    #img = mpimg.imread("..\\img\\"+ image +".jpg")
    #img = img[:, :, 1]
    img = Image.open("..\img\\" + image).convert("L")
    fig, ax = plt.subplots(2, 2)
    ax[0, 0].imshow(img, cmap = plt.get_cmap('gray'))
    ax[0, 0].set_title("Original")

    #fig2 = plt.figure()
    #fig2 = plt.hist(img.ravel(), bins = 256, range=(0, 256))

    iW = dwt2d(np.copy(img), MAX_LAYER = ml, wavelet = wavelet)
    ax[1, 0].imshow(iW, cmap = plt.get_cmap('gray'))
    ax[1, 0].set_title("Wavelet")

    irW = wavelet_reorder_2d(iW, MAX_LAYER = ml)
    ax[1, 1].imshow(irW, cmap = plt.get_cmap('gray'))
    ax[1, 1].set_title("Deinterlaced")

    tW = idwt2d(np.copy(iW), MAX_LAYER = ml, wavelet = wavelet)
    ax[0, 1].imshow(tW, cmap = plt.get_cmap('gray'))
    ax[0, 1].set_title("Recovered")
    
    plt.show()

if __name__ == "__main__":
    imtest(image = "flower.jpg", wavelet = "haar", ml = 3)
    #imtest(image = "circle_stripe.png", wavelet = "haar", ml = 2)
