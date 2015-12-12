import pdb
import numpy as np
import matplotlib.pyplot as plt


# --- wavelet transforms ---
def db2(S, MAX_LAYER = -1):
    """forward discrete wavelet transform db2/haar (1 vanishing moment)"""
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
    """inverse discrete wavelet transform db2/haar (1 vanishing moment)"""
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
    """forward discrete wavelet transform db4 (2 vanishing moments)"""
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
    """inverse discrete wavelet transform db4 (2 vanishing moments)"""
    #same as forward, except all additions and subtractions are flipped
    #W: array full of wavelet coefficients
    #MAX_LAYER: total level of coefficients in transform
    #  It is possible for a transform to be partial, such as data with size 2^16 to go 8 levels deep instead of 16

    N = len(W)
    s3 = np.sqrt(3)
    s2 = np.sqrt(2)

    layers = 0
    step = N

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

def db6(S, MAX_LAYER = -1):
    """forward discrete wavelet transform db6 (3 vanishing moments)"""
    #S: array full of data
    #T: number of elements evaluated

    N = len(S)
    T = N

    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N        S.resize(p2, refcheck=False)
        print S
        T = p2

    #coefficients
    alpha = -0.4122865950
    beta = -1.5651362796
    beta2 = 0.3523876576
    gamma = 0.0284590896
    gamma2 = 0.4921518449
    delta = -0.3896203900
    zeta = 1.9182029462

    layers = 0
    step = 2

    #periodic extension used, boundaries do warp

    while step < 2*T:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break

        start = step/2

        #update 1
        #u1(z) = alpha
        #even[n] += -alpha * odd[n]
        S[::step]  += -alpha * S[start::step]
        
        #predict 1
        #p1(z) = beta * z + beta'
        #even[n] += -beta * odd[n+1] - beta' * odd[n]

        #S[start:-step:step] += -beta * S[step::step] - beta2 * S[:-step:step]
        #S[-step] += -beta * S[0] - beta2 * S[-step]
        S[:-step:step] += -beta * S[start+step::step] - beta2 * S[start:-step:step]
        S[-step] += -beta * S[start] - beta2 * S[-start]

        
        #update 2
        #u2(z) = gamma + gamma'z^-1
        #odd[n] += -gamma * even[n] - gamma' * even[n-1]

        S[start+step::step] += -gamma * S[step::step] - gamma2 * S[:-step:step]
        S[start] += -gamma * S[0] - gamma2 * S[-step]
        
        
        #predict 2
        #p2(z) = delta
        #even[n] += -delta * odd[n]
        S[::step] += -delta * S[start::step]
        
        
        #normalize
        #even[n] *= 1/zeta
        #odd[n]  *= zeta
        S[::step]  *= 1/zeta
        S[start::step] *= zeta
        
        step <<= 1
        layers += 1

    return S

def idb6(W, MAX_LAYER = -1):
    """inverse discrete wavelet transform db6 (3 vanishing moments)"""
    #same as forward, except all additions and subtractions are flipped
    #W: array full of wavelet coefficients
    #MAX_LAYER: total level of coefficients in transform
    #  It is possible for a transform to be partial, such as data with size 2^16 to go 8 levels deep instead of 16

    N = len(W)
    s3 = np.sqrt(3)
    s2 = np.sqrt(2)

    layers = 0
    step = N

    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >= MAX_LAYER >= 0): step = 1 << MAX_LAYER

    #coefficients
    alpha = -0.4122865950
    beta = -1.5651362796
    beta2 = 0.3523876576
    gamma = 0.0284590896
    gamma2 = 0.4921518449
    delta = -0.3896203900
    zeta = 1.9182029462

    while step > 1:
        start = step/2
        
        #normalize'
        #even[n] *= zeta
        #odd[n]  *= 1/zeta
        W[::step]  *= zeta
        W[start::step] *= 1/zeta
        
        
        #predict 2'
        #p2'(z) = -delta
        #even[n] += delta * odd[n]
        W[::step] += delta * W[start::step]
        
        #update 2'
        #u2'(z) = -gamma + -gamma'z^-1
        #odd[n] += gamma * even[n] + even' * odd[n-1]

        W[start+step::step] += gamma * W[step::step] + gamma2 * W[:-step:step]
        W[start] += gamma * W[0] + gamma2 * W[-step]
        
        #predict 1'
        #p1'(z) = -beta * z - beta'
        #even[n] += beta * odd[n+1] + beta' * odd[n]

        W[:-step:step] += beta * W[start+step::step] + beta2 * W[start:-step:step]
        W[-step] += beta * W[start] + beta2 * W[-start]
        
        #update 1'
        #u1'(z) = -alpha
        #odd[n] += alpha * even[n]
        W[::step]  += alpha * W[start::step]
        
        step >>= 1

    return W

def cdf97(S, MAX_LAYER = -1):
    """CDF9/7 forward transformation, used for JPEG2000"""

    N = len(S)
    T = N

    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N        S.resize(p2, refcheck=False)
        print S
        T = p2

    #coefficients
    alpha = -1.5861343420693648
    beta = -0.0529801185718856
    gamma = 0.8829110755411875
    delta = 0.4435068520511142
    kodd = 1/1.1496043988602418
    keven = 1.1496043988602418

    layers = 0
    step = 2

    while step < 2*T:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break

        start = step/2

        #update 1
        #u1(z) = alpha(1 + z^-1)
        #odd[n] += alpha * (even[n] + even[n+1])
        S[start:-start:step]  += alpha * (S[:-step:step] + S[step::step])
        S[-start] += alpha * (S[-step] + S[0])
        
        #predict 1
        #p1(z) = beta(1 + z)
        #even[n] += beta * (odd[n] + odd[n-1])

        S[step::step] += beta * (S[step+start::step] + S[start:-start:step])
        S[0] += beta * (S[start] + S[-start])
        
        #update 2
        #u2(z) = gamma*(1 + z^-1)
        #odd[n] += gamma * (even[n] + even[n+1]

        S[start:-start:step] += gamma * (S[step::step] + S[:-step:step])
        S[-start] += gamma * (S[-step] +  S[0])

        #predict 2
        #p2(z) = delta * (1 + z)
        #even[n] += delta * (odd[n] + odd[n-1])
        S[step::step] += delta * (S[start+step::step] + S[start:-start:step])
        S[0] += delta * (S[start] + S[-start])

        #normalize
        #even[n] *= keven
        #odd[n]  *= kodd
        S[::step]  *= keven
        S[start::step] *= kodd
        
        step <<= 1
        layers += 1

    return S

def icdf97(W, MAX_LAYER = -1):
    """CDF9/7 inverse transformation, used for JPEG2000"""
    #coefficients
    alpha = -1.5861343420693648
    beta = -0.0529801185718856
    gamma = 0.8829110755411875
    delta = 0.4435068520511142
    kodd = 1/1.1496043988602418
    keven = 1.1496043988602418
    
    N = len(W)
    step = N


    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >= MAX_LAYER >= 0): step = 1 << MAX_LAYER
    while step > 1:
        start = step/2
        
        #normalize'
        #even[n] /= keven
        #odd[n]  /= kodd
        W[::step]  /= keven
        W[start::step] /= kodd


        #predict 2'
        #p2(z) = delta * (1 + z)
        #even[n] -= delta * (odd[n] + odd[n-1])
        W[step::step] -= delta * (W[start+step::step] + W[start:-start:step])
        W[0] -= delta * (W[start] + W[-start])

        #update 2'
        #u2(z) = gamma*(1 + z^-1)
        #odd[n] -= gamma * (even[n] + even[n+1]

        W[start:-start:step] -= gamma * (W[step::step] + W[:-step:step])
        W[-start] -= gamma * (W[-step] +  W[0])
        
        #predict 1'
        #p1(z) = beta(1 + z)
        #even[n] -= beta * (odd[n] + odd[n-1])

        W[step::step] -= beta * (W[step+start::step] + W[start:-start:step])
        W[0] -= beta * (W[start] + W[-start])
        
        #update 1
        #u1(z) = alpha(1 + z^-1)
        #odd[n] -= alpha * (even[n] + even[n+1])
        W[-start] -= alpha * (W[-step] + W[0])
        W[start:-start:step]  -= alpha * (W[:-step:step] + W[step::step])


        step >>= 1

    return W

def legall53(S, MAX_LAYER = -1):
    """Legall 5/3 forward transformation, used for JPEG2000"""

    N = len(S)
    T = N

    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N        S.resize(p2, refcheck=False)
        T = p2

    layers = 0
    step = 2

    while step < 2*T:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break

        start = step/2

        #update 1
        #u1(z) = (1 + z^-1)/2
        S[start:-start:step] -= (S[step::step] + S[:-step:step])/2
        S[-start] -= (S[-step] + S[0])/2

        #predict 1
        #p1(z) = (1 + z)/4
        S[step::step] += (S[step+start::step] + S[start:-start:step])/4
        S[0] += (S[start] + S[-start])/4

        step <<= 1
        layers += 1

    return S

def ilegall53(W, MAX_LAYER = -1):
    """Legall 5/3 inverse transformation, used for JPEG2000"""    
    layers = 0
    N = len(W)
    step = N


    #if a max number of layers is specified, lower step so the high freq are processed
    if (((N).bit_length() - 1) >= MAX_LAYER >= 0): step = 1 << MAX_LAYER

    while step > 1:
        start = step/2

        #predict 1'
        #p1(z) = (1 + z)/4
        W[step::step] -= (W[step+start::step] + W[start:-start:step])/4
        W[0] -= (W[start] + W[-start])/4

        #update 1'
        #u1(z) = (1 + z^-1)/2
        W[start:-start:step] += (W[step::step] + W[:-step:step])/2
        W[-start] += (W[-step] + W[0])

        step >>= 1

    return W

def dwt(S, MAX_LAYER = -1, wavelet = "db4"):
    """wrapper for forward discrete wavelet transform"""
    if wavelet == "db2" or wavelet == "haar":
        return db2(S,  MAX_LAYER = MAX_LAYER)
    elif wavelet == "db4":
        return db4(S, MAX_LAYER = MAX_LAYER)
    elif wavelet == "db6":
        return db6(S, MAX_LAYER = MAX_LAYER)
    elif wavelet == "cdf97":
        return cdf97(S, MAX_LAYER = MAX_LAYER)
    elif wavelet == "legall53":
        return legall53(S, MAX_LAYER = MAX_LAYER)
    else:
        return db4(S, MAX_LAYER = MAX_LAYER)

def idwt(W, MAX_LAYER = -1, wavelet = "db4"):
    """wrapper for inverse discrete wavelet transform"""
    if wavelet == "db2" or wavelet == "haar":
        return idb2(W,  MAX_LAYER = MAX_LAYER)
    elif wavelet == "db4":
        return idb4(W, MAX_LAYER = MAX_LAYER)
    elif wavelet == "db6":
        return idb6(W, MAX_LAYER = MAX_LAYER)
    elif wavelet == "cdf97":
        return icdf97(W, MAX_LAYER = MAX_LAYER)
    elif wavelet == "legall53":
        return ilegall53(W, MAX_LAYER = MAX_LAYER)
    else:
        return idb4(W, MAX_LAYER = MAX_LAYER)

# --- coefficient ordering ---
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

# --- thresholding and denoising ---
def threshold(W, LAM = -1,  MAX_LAYER = -1, N_ELEMENTS = 0, shrink = "hard", policy = "none"):
    #W: array full of wavelet coefficients
    #Possible flags, choice between:
    #   LAM: set threshold value
    #   MAX_LAYER: how many layers in transform
    #   N_ELEMENTS: preserve number of wavelet coefficients
    #   SHRINK: "hard", "soft", "semisoft"
    #   POLICY: none, universal, SURE

    #if THRESHOLD is an integer, THRESHOLD is applied on every level
    #if it is an array, then THRESHOLD is level dependent

    #Warning: Horrible code. Need to completely rehaul.

    shrink = shrink.lower()
    policy = policy.lower()
    N = len(W)
    logN = (N).bit_length() - 1

    shrinkage = False

    if LAM > 0:
        THRESHOLD = LAM

    elif policy != "none":

        start = 1
        step = 2

        if policy == "universal":
            #universal thresholding is the upper bound
            #doesnt matter whether hard, soft, or semisoft are used
            median = np.median(np.abs(W))
            lnN = np.log(N)
            THRESHOLD = np.abs(median/0.6745*np.sqrt(2.0*lnN))


        elif policy == "universal_level" and MAX_LAYER > 0:
            layer = 0
            THRESHOLD = []
            stepmax = step << MAX_LAYER
            currN = N >> 1
            while step < stepmax:
                ln_currN = np.log(currN)
                median = np.median(np.abs(W[start::step]))
                thresh = np.abs(median/0.6745 * np.sqrt(2 * ln_currN))
                THRESHOLD.append(thresh)

                start <<= 1
                step  <<= 1
                currN >>= 1


        elif policy == "neighcoeff" and MAX_LAYER > 0:
            stepmax = step << MAX_LAYER
            currN = N >> 1
            L = 3
            THRESHOLD = "Neighcoeff Override!"
            while step < stepmax:
                lam = 2/3.0 * np.log(currN)
                #use median absolute deviation

                #WARNING
                #BAD CODE HERE
                roller = np.empty((currN, 3))
                roller[:, 0] = np.roll(W[start::step], -1)
                roller[:, 1] = W[start::step]
                roller[:, 2] = np.roll(W[start::step], 1)

                s2 = np.sum(roller**2, axis = 1)

                sigma = np.median(np.abs(roller), axis = 1)/0.6745
                #sigma = np.median(np.abs(np.diff(tw) - np.median(np.diff(tw)))) / np.sqrt(2) / 0.6745

                beta = 1 - lam * L * sigma**2 / s2
                beta[beta < 0] = 0

                W[start::step] *= beta

                start <<= 1
                step  <<= 1
                currN >>= 1


        """        
        elif policy == "sure" and MAX_LAYER > 0:
            #sure stuff goes here
            #requires soft thresholding
            THRESHOLD = []
            shrink = "soft"
            layer = 0
            while layer < MAX_LAYER-1 and i > 0:
                #pdb.set_trace()
                w_target = W[i:j]
                sigma = np.median(np.abs(w_target))/0.6745
                s = lambda t: sure(w_target, t, sigma)
                lni = np.log(i)
                t0 = sigma*np.sqrt(2.0*lni)
                #threshold minimizes SURE
                #start with universal/2
                t = opt.minimize(fun = s, x0 = t0/2, bounds = [[0, t0]])
                nu = sure_sparsity(w_target)
                if nu <= 1: t = t0
                THRESHOLD.append(t)
                i >>= 1
                j >>= 1
        """

    elif N_ELEMENTS > 0:
        #THRESHOLD = sorted([abs(i) for i in W], reverse=True)[N_ELEMENTS]
        THRESHOLD = np.abs(W[np.argpartition(-np.abs(W), N_ELEMENTS)[N_ELEMENTS]])

        shrink = "hard"


    else:
        return W


    if type(THRESHOLD) != list and type(THRESHOLD) != np.ndarray:
       # print "Constant Threshold:", THRESHOLD
        W[np.abs(W) < THRESHOLD] = 0.0
        if shrink == "soft":
            W[np.abs(W) >= THRESHOLD] -= THRESHOLD * np.sign(W[np.abs(W) >= THRESHOLD])
           # print "soft threshold"

        elif shrink == "semisoft":
            W[np.abs(W) >= THRESHOLD] = np.sqrt(W[np.abs(W) >= THRESHOLD]**2 - THRESHOLD**2) * np.sign(W[np.abs(W) >= THRESHOLD])


    elif type(THRESHOLD) == list:
        #level dependent thresholding is done from the highest frequency to the lowest frequency

        start = 1
        step  = 2

        if MAX_LAYER < 0: MAX_LAYER = logN
        i = 0

        while i < MAX_LAYER:
            t = THRESHOLD[i]

            W[start::step][np.abs(W[start::step]) < t] = 0.0

            if shrink == "soft":
                W[start::step][np.abs(W[start::step]) >= t] -= t * np.sign(W[start::step][abs(W[start::step]) >= t])
                #print "soft threshold"

            elif shrink == "semisoft":
                W[start::step][abs(W[start::step]) >= t] = np.sqrt(W[start::step][abs(W[start::step]) >= t]**2 - t**2) * np.sign(W[start::step][abs(W[start::step]) >= t])


            i += 1
            start <<= 1
            step  <<= 1

    else:
        print THRESHOLD

    return W

# --- printing and plotting ---
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


# --- test functions ---
def doppler(t):
    return np.sqrt(t*(t+1)) * np.sin(2*np.pi*1.05/(t+0.05))

def sinusoid(t):
    return np.sin(2 * np.pi * 4 * t)

def blocks(t):
    K = lambda c: (1 + np.sign(c))/2.0
    x = np.array([[.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81]]).T
    h = np.array([[4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2]]).T
    return np.sum(h*K(t/2 - x), axis=0)

def bumps(t):
    K = lambda c : (1. + np.abs(c)) ** -4.
    x = np.array([[.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81]]).T
    h = np.array([[4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 2.1, 4.2]]).T
    w = np.array([[.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005]]).T
    return np.sum(h*K((x-t/2)/w), axis=0)

def heavisine(t):
    return 4 * np.sin(2*np.pi*t) - np.sign(t/2 - 0.3) - np.sign(0.72 - t/2)


# --- plots and accuracy ---
def psnr(truth, estimate):
    #returns Peak Signal to Noise Ratio in decibels
    N = len(truth)
    mean = np.mean(truth)
    #num = sum([abs(truth[i] - mean)**2 for i in range(N)])
    num = np.sum(truth ** 2)
    denom = np.sum((truth-estimate)**2)
    snr = 10 * np.log10(num/denom)

    return snr

def approx_plot(t, W, MAX_LAYER = -1, wavelet = "db4", f = []):
    #sequence of approximations plotted, each with 1/2 as much data as the previous
    N = len(W)
    logN = (N).bit_length() - 1

    #if a max number of layers is specified, lower step so the high freq are processed
    if (logN >= MAX_LAYER >= 0):
        step = 1 << MAX_LAYER
    else:
        step = N
        MAX_LAYER = logN


    if wavelet == "db2" or wavelet == "haar":
        factor = 1
    elif wavelet == "db4":
        factor = 1/np.sqrt(2)
    elif wavelet == "db6":
        factor = 1/np.sqrt(2)
    elif wavelet == "cdf97":
        factor = 1/np.sqrt(2)
    elif wavelet == "legall53":
        factor = 1

    curr_factor = 1.0
    step = 1
    fig = plt.figure()
    ax = fig.add_subplot(111)

    prev = 1

    if f != []:
        ax.plot(t, f, label = "truth")
    
    for level in range(MAX_LAYER, 0, -1):
        #use approximation
        ids = idwt(np.copy(W[::step]), level, wavelet)
        ids *= curr_factor
        #print level, np.max(ids), prev/np.max(ids)
        prev = np.max(ids)
        ts = t[step/2::step]
        
        name = "level" + str(level)
        ax.plot(ts, ids, label = name)

        curr_factor *= factor
        step <<= 1

    ax.legend(loc = "lower right")
    ax.set_title(wavelet + " Approximation")
    return fig


def coeff_pyramid(t, W, MAX_LAYER = -1, wavelet = "db4"):
    #2 plots, a coefficient stem plot, and a plot showing the bases (inverse of each layer of coefficients)
    N = len(W)
    logN = (N).bit_length() - 1

    #if a max number of layers is specified, lower step so the high freq are processed
    if (logN >= MAX_LAYER >= 0):
        step = 1 << MAX_LAYER
    else:
        step = N
        MAX_LAYER = logN

    #wavelet coefficients
    fig, axs = plt.subplots(MAX_LAYER+1, 1, sharex = True)
    fig.suptitle(wavelet + " wavelet coefficients", fontsize = 20)
    
    #basis functions
    fig2, axs2 = plt.subplots(MAX_LAYER+1, 1, sharex = True, sharey = True)
    fig2.suptitle(wavelet + " wavelet basis functions", fontsize = 20)

    i = 0

    T = N >> (MAX_LAYER)
    Fs = 1 / (t[1]-t[0])

    currplot = 0
    max_freq = Fs / (1 << (MAX_LAYER - 1))

    while step > 1:
        start = step >> 1
        layer_t = t[start::step]

        w_empty = np.zeros(N)

        if i == 0:
            title_str = "Scaling Coeff (0 - {:.2f} Hz)".format(Fs /(1 << (MAX_LAYER - 1)))
            layer_w = W[::step]
            axs[i].set_title(title_str)
            axs2[i].set_title(title_str)
            w_empty[::step] = W[::step]

        else:
            title_str = "Wavelet Coeff level {:x} ({:.2f} - {:.2f} Hz)".format(i-1, max_freq, max_freq*2)
            axs[i].set_title(title_str)
            axs2[i].set_title(title_str)

            max_freq *= 2
            layer_w = W[start::step]
            w_empty[start::step] = W[start::step]

            step >>= 1

        axs[i].stem(layer_t, layer_w, basefmt = "black", markerfmt = " ")

        axs2[i].axhline(y = 0, color = 'k')
        idw = idwt(np.copy(w_empty), MAX_LAYER, wavelet)
        axs2[i].plot(t, idw)


        i += 1

    return [fig, fig2]

if __name__ == "__main__":
    #s = np.array([32.0, 10.0, 20.0, 38.0, 37.0, 28.0, 38.0, 34.0, 18.0, 24.0, 18.0, 9.0, 23.0, 24.0, 28.0, 34.0])
    N = 1024
    x = np.linspace(0, 2, N)
    s = doppler(x)

    
    ml = 7
    wv = "legall53"
    #print s
    ds = dwt(np.copy(s), MAX_LAYER = ml, wavelet = wv)
    #waveletprint(ds, MAX_LAYER = ml)
    ids = idwt(np.copy(ds), MAX_LAYER = ml, wavelet = wv)
    #print ids
    approx_plot(x, np.copy(ds), MAX_LAYER = ml, wavelet = wv, f = s)
    coeff_pyramid(x, ds, MAX_LAYER = ml, wavelet = wv)
    
    """
    x = np.linspace(0, 2, N)
    sigma = 0.30
    noise = np.random.normal(loc = 0.0, scale = sigma, size = N)
    #s = heavisine(x)
    rs = s + noise
    ml = 1
    wavelet_type = "cdf97"
    #print s
    ds = dwt(np.copy(s), MAX_LAYER = ml, wavelet = wavelet_type)
    drs = dwt(np.copy(rs), MAX_LAYER = ml, wavelet = wavelet_type)
    #print ds
    #print interlace_to_layer(ds)
    ids = idwt(np.copy(ds), MAX_LAYER = ml, wavelet = wavelet_type)
    #print ids

    #waveletprint(ds, MAX_LAYER = ml)
    #coeff_pyramid(x, ds, MAX_LAYER = ml, wavelet = wavelet_type)
    #approx_plot(x, ds, MAX_LAYER = ml, wavelet = wavelet_type, f = s)
    """

    """
    #coeff_pyramid(x, drs, MAX_LAYER = ml, wavelet = wavelet_type)

    coeff_univ_level = threshold(np.copy(drs), policy = "universal_level", MAX_LAYER = 4)
    coeff_univ = threshold(np.copy(drs), policy = "universal")
    coeff_neighcoeff = threshold(np.copy(drs), policy = "neighcoeff", MAX_LAYER = 4)
    #coeff_pyramid(x, drst, MAX_LAYER = ml, wavelet = wavelet_type)
    #approx_plot(x, ds2, MAX_LAYER = ml, wavelet = wavelet_type)

    univ_level = idwt(np.copy(coeff_univ_level), MAX_LAYER = ml, wavelet = wavelet_type)
    univ = idwt(np.copy(coeff_univ), MAX_LAYER = ml, wavelet = wavelet_type)
    neighcoeff = idwt(np.copy(coeff_neighcoeff), MAX_LAYER = ml, wavelet = wavelet_type)
    psnr_truth = psnr(s, ids)
    print "True PSNR:\t", psnr_truth
    psnr_noise = psnr(s, rs)
    print "Random PSNR:\t", psnr_noise
    psnr_univ = psnr(s, univ)
    print "Universal PSNR:\t", psnr_univ
    psnr_univ_level = psnr(s, univ_level)
    print "UnivLevel PSNR:\t", psnr_univ_level
    psnr_neighcoeff = psnr(s, neighcoeff)
    print "Neighcoeff PSNR:\t", psnr_neighcoeff

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(r"Denoising (AWGN, $\sigma$={:.2f})".format(sigma))
    ax.plot(x, s, label = "Ground Truth ({:.2f} dB)".format(psnr_truth))
    ax.plot(x, rs, label = "Noisy Data ({:.2f} dB)".format(psnr_noise))
    ax.plot(x, univ, label = "Univ Threshold ({:.2f} dB)".format(psnr_univ))
    ax.plot(x, univ_level, label = "Univ Level Threshold ({:.2f} dB)".format(psnr_univ_level))
    ax.plot(x, neighcoeff, label = "Neighcoeff Threshold ({:.2f} dB)".format(psnr_neighcoeff))
    ax.legend(loc = "lower right")
    """

    plt.show()
    
