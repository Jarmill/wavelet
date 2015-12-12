import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pdb

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
    N = S.shape[0]
    T = N
    layers = 0
    step = 2
    vert = False

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

    N = W.shape[0]
    step = N
    vert = True

    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >= MAX_LAYER >= 0): step = 1 << MAX_LAYER
    
    while step > 1:
        start = step/2

        #predict 1'
        #p1(z) = (1 + z)/4
        W[:, step::step] -= np.floor((W[:, step+start::step] + W[:, start:-start:step])/4)
        W[:, 0] -= np.floor((W[:, start] + W[:, -start])/4)

        #update 1'
        #u1(z) = (1 + z^-1)/2
        W[:, start:-start:step] += np.floor((W[:, step::step] + W[:, :-step:step])/2)
        W[:, -start] += np.floor((W[:, -step] + W[:, 0])/2)

        if vert:
            step <<= 1
            vert = False
            W = W.T
        else:
            vert = True
            W = W.T

        step >>= 1

    return W

def wavelet_reorder_2d(W, MAX_LAYER = -1):
    """print non-interlaced wavelet coefficients"""
    c = np.zeros_like(W)
    N = W.shape[0]
    layers = 1
    step = 2
    start = step/2

    end = N
    mid = end/2
    
    while step < 2*N:
        if (MAX_LAYER >= 1 and layers > MAX_LAYER): break        
        c[mid:end, mid:end] = W[start::step, start::step]   #HH
        c[:mid, mid:end] = W[start::step, ::step]           #LH
        c[mid:end, :mid] = W[::step, start::step]                 #HL

        end /= 2
        mid/= 2
        step *= 2
        start *=2
        layers += 1
    c[:end, :end] = W[::start, ::start]

    
    return c

def dwt2d(S, MAX_LAYER = 3, wavelet = "legall53_2d"):
    """wrapper for 2d forward discrete wavelet transform"""
    if wavelet == "legall53_2d":
        return legall53_2d(S, MAX_LAYER)
    else:
        return legall53_2d(S, MAX_LAYER)

def idwt2d(S, MAX_LAYER = 3, wavelet = "legall53_2d"):
    """wrapper for 2d forward discrete wavelet transform"""
    if wavelet == "legall53_2d":
        return ilegall53_2d(S, MAX_LAYER)
    else:
        return ilegall53_2d(S, MAX_LAYER)

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

def imtest(image = "dock", wv = "legall53_2d", ml = 2):
    img = mpimg.imread("..\\img\\"+ image +".jpg")
    img = img[:, :, 1]
    fig, ax = plt.subplots(2, 2)
    ax[0, 0].imshow(img, cmap = plt.get_cmap('gray'))
    ax[0, 0].set_title("Original")

    #fig2 = plt.figure()
    #fig2 = plt.hist(img.ravel(), bins = 256, range=(0, 256))

    iW = dwt2d(np.copy(img), MAX_LAYER = ml, wavelet = "legall53_2d")
    ax[1, 0].imshow(iW, cmap = plt.get_cmap('gray'))
    ax[1, 0].set_title("Wavelet")

    irW = wavelet_reorder_2d(iW, MAX_LAYER = ml)
    ax[1, 1].imshow(irW, cmap = plt.get_cmap('gray'))
    ax[1, 1].set_title("Deinterlaced")

    tW = idwt2d(np.copy(iW), MAX_LAYER = ml, wavelet = "legall53_2d")
    ax[0, 1].imshow(tW, cmap = plt.get_cmap('gray'))
    ax[0, 1].set_title("Recovered")
    
    plt.show()


if __name__ == "__main__":
    imtest(image = "dock", ml = 2)
