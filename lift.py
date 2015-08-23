#all operations use db4 wavelet, lifting scheme
#algorithm from http://www.bearcave.com/misl/misl_tech/wavelets/daubechies/
import math
import pdb
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np

def db2(S, MAX_LAYER = -1):
    #forward discrete wavelet transform db2 (1 vanishing moment)
    #S: array full of data
    #T: number of elements evaluated

    #split array, first half are even elements and second half are odd element
    N = len(S)
    T = N
    
    #if N is not a power of 2
    if N & (N-1) != 0 and N > 0:
        #find lowest power of 2 greater than N
        p2 = 1 << (N).bit_length()
        #pad with zeros up to N
        #S += [0.0]*(p2-N)
        #make sure that array can be resized
        S.resize(p2, refcheck=False)
        print S
        T = p2
    
    layers = 0
    
    while T > 1:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break
        
        #split
        half = T/2
        Spart = np.empty(T)
        Spart[:half] = S[:T:2]
        Spart[half:] = S[1:T:2]
        S[:T] = Spart        

        #predict 1
        S[half:T] -= S[:half]
        
        #update 1
        S[:half] += S[half:T]/2              
         
        T >>= 1
        layers += 1
    
    return S

def idb2(W, MAX_LAYER = -1):
    N = len(W)
    T = 2
    
    layer = 0
    
    if MAX_LAYER >= 0 and MAX_LAYER < (N).bit_length():
        #print MAX_LAYER, N
        T = N >> (MAX_LAYER-1)
    
    while T <= N:
    
        half = T/2
                
        #update' 1
        W[:half] -= W[half:T]/2
        
        #predict' 1
        W[half:T] += W[:half]

        #merge
        half = T/2
        Wpart = np.empty(T)
        Wpart[::2] = W[:T/2]
        Wpart[1::2] = W[T/2:T]
        W[:T] = Wpart
        
        T <<= 1
    
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
        #pad with zeros up to N
        #S += [0.0]*(p2-N)
        #make sure that array can be resized
        S.resize(p2, refcheck=False)
        print S
        T = p2
    
    layers = 0

    s3 = math.sqrt(3)
    s2 = math.sqrt(2)
    
    while T > 1:
        if (MAX_LAYER >= 0) and (layers >= MAX_LAYER): break
        
        #split
        half = T/2
        Spart = np.empty(T)
        Spart[:half] = S[:T:2]
        Spart[half:] = S[1:T:2]
        S[:T] = Spart
        #S = S[:T:2] + S[1:T:2] + S[T:]
        
        
        #update
        S[:half] += s3 * S[half:T]
        
        #predict
        S[(half+1):T] += - s3/4*S[1:(half)] - (s3-2)/4*S[:(half-1)]
        S[half] += - s3/4*S[0] - (s3-2)/4*S[half-1]
        
        #update 2
        S[:(half-1)] = S[:(half-1)] - S[(half+1):T]
        S[half-1] = S[half-1] - S[half]
        
        #normalize
        S[:half] = (s3-1)/s2 * S[:half]
        S[half:T] = (s3+1)/s2*S[half:T]
        
        layers += 1
        
        T /= 2
        
    return S
    
def idb4(W, MAX_LAYER = -1):
    #inverse discrete wavelet transform    
    #same as forward, except all additions and subtractions are flipped
    #W: array full of wavelet coefficients
    #MAX_LAYER: total level of coefficients in transform
    #  It is possible for a transform to be partial, such as data with size 2^16 to go 8 levels deep instead of 16
    #RECOVER_LAYER: How many layers desired for recovery, new array will be truncated when returned
    
    N = len(W)
    s3 = math.sqrt(3)
    s2 = math.sqrt(2)

    T = 2
    
    layer = 0
    
    if MAX_LAYER >= 0 and MAX_LAYER < (N).bit_length():
        #print MAX_LAYER, N
        T = N >> (MAX_LAYER-1)
    
    while T <= N:
        half = T/2

        #normalize'
        W[:half] = (s3+1)/s2 * W[:half]
        W[half:T] = (s3-1)/s2*W[half:T]   
        
        #update 2'
        W[:(half-1)] = W[:(half-1)] + W[(half+1):T]
        W[half-1] = W[half-1] + W[half]
        
        #predict'        
        W[(half+1):T] += + s3/4*W[1:(half)] + (s3-2)/4*W[:(half-1)]
        W[half] += + s3/4*W[0] + (s3-2)/4*W[half-1]
        
        #update 1'
        W[:half] += -s3 * W[half:T]

        #merge odd and even parts together
        #[0,2,4,1,3,5] to [0,1,2,3,4,5]
        Wpart = np.empty(T)
        Wpart[::2] = W[:T/2]
        Wpart[1::2] = W[T/2:T]
        W[:T] = Wpart
        
        T <<= 1
            
    return W
def dwt(S, MAX_LAYER = -1, wavelet = "db4"):
    if wavelet == "db2" or wavelet == "haar":
        return db2(S,  MAX_LAYER = MAX_LAYER)
    elif wavelet == "db4":
        return db4(S, MAX_LAYER = MAX_LAYER)
    #elif wavelet == "db6":
    #    return db6(S, MAX_LAYER = MAX_LAYER)
    else:
        return db4(S, MAX_LAYER = MAX_LAYER)
    
def idwt(W, MAX_LAYER = -1, wavelet = "db4"):
    if wavelet == "db2" or wavelet == "haar":
        return idb2(W,  MAX_LAYER = MAX_LAYER)
    elif wavelet == "db4":
        return idb4(W, MAX_LAYER = MAX_LAYER)
    elif wavelet == "db6":
        return idb6(W, MAX_LAYER = MAX_LAYER)
    else:
        return idb4(W, MAX_LAYER = MAX_LAYER)    
    
def waveletprint(W, MAX_LAYER = -1):
    N = len(W)
    i = 0
    sc = 0
    max_i = 0
    logn = (N).bit_length()-1
    if logn > MAX_LAYER  >= 0:
        sc = logn-MAX_LAYER
        max_i = N >> (MAX_LAYER+1)
    j = 2*max_i if max_i != 0 else 1
    
    while i < N:
        st = "\nScale 2^-"+str(sc)+":"
        print st, W[i:j]
        i = j
        j = j << 1
        sc += 1
    print "\n"
    
def sure(W, THRESHOLD, sigma):
    #Stein's Unbiased Risk Estimator
    #must be minimized
    risk = 0
    var = sigma*sigma
    for x in W:
        if x <= THRESHOLD:
            risk += x*x - var
        else:
            risk += var + THRESHOLD*THRESHOLD
    return risk

def sure_sparsity(W):
    #measure of sparsity of coefficient layer
    N = float(len(W))
    num = np.sum(W**2 - 1)
    denom = math.log(N, 2)**1.5
    if denom == 0: denom = 1e-64
    nu = num/(denom * math.sqrt(N))
    
    return nu

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
        
        i = N/2
        j = N
        median = np.median(abs(W))
        lni = (i).bit_length() - 1
        univ_threshold = median/0.6745*math.sqrt(2.0*lni)
        #print "\nMEDIAN:", median, "\nTHRESHOLD:", THRESHOLD
        
        if policy == "universal":
            #universal thresholding is the upper bound
            #doesnt matter whether hard, soft, or semisoft are used
            median = np.median(abs(W[i:]))
            #THRESHOLD = univ_threshold
            THRESHOLD = median/0.6745*math.sqrt(2.0*lni)
        
        elif policy == "universal_level" and MAX_LAYER > 0:
            layer = 0
            THRESHOLD = []
            while layer < MAX_LAYER-1 and i > 0:
                lni = (i).bit_length() - 1
                median = np.median(abs(W[i:j]))
                THRESHOLD.append(median/0.6745*math.sqrt(2.0*lni))
                i >>= 1
                j >>= 1
                
        elif policy == "sure" and MAX_LAYER > 0:
            #sure stuff goes here
            #requires soft thresholding
            THRESHOLD = []
            shrink = "soft"
            layer = 0
            while layer < MAX_LAYER-1 and i > 0:
                #pdb.set_trace()
                w_target = W[i:j]
                sigma = np.median(abs(w_target))/0.6745
                s = lambda t: sure(w_target, t, sigma)
                lni = math.log(i)
                t0 = sigma*math.sqrt(2.0*lni)
                #threshold minimizes SURE
                #start with universal/2
                t = opt.minimize(fun = s, x0 = t0/2, bounds = [[0, t0]])
                nu = sure_sparsity(w_target)
                if nu <= 1: t = t0
                THRESHOLD.append(t)
                i >>= 1
                j >>= 1
                
        
        elif policy == "neighcoeff" and MAX_LAYER > 0:
            while layer < MAX_LAYER - 1 and i > 0:
                pass
                                
        
    elif N_ELEMENTS > 0:
        #THRESHOLD = sorted([abs(i) for i in W], reverse=True)[N_ELEMENTS]
        THRESHOLD = W[np.argpartition(-abs(W), N_ELEMENTS)]
        
        shrink = "hard"
        
        
    else:
        return W
    
    
    if type(THRESHOLD) != list and type(THRESHOLD) != np.ndarray:
       # print "Constant Threshold:", THRESHOLD
        W[np.abs(W) < THRESHOLD] = 0.0
        if shrink == "soft":
            W[abs(W) >= THRESHOLD] -= THRESHOLD * np.sign(W[abs(W) >= THRESHOLD])
           # print "soft threshold"
        
        elif shrink == "semisoft":
            W[abs(W) >= THRESHOLD] = np.sqrt(W[abs(W) >= THRESHOLD]**2 - THRESHOLD**2) * np.sign(W[abs(W) >= THRESHOLD])
        
    elif type(THRESHOLD) == list:
        #TODO
        i = N/2
        j = N
        c = 0
        if MAX_LAYER < 0: MAX_LAYER = logN
        while c < MAX_LAYER and c < len(THRESHOLD) and i > 0:
            #repeated indexing is compiler optimized away
            t = THRESHOLD[c]
            #print t
            W[i:j][np.abs(W[i:j]) < t] = 0.0
            if shrink == "soft":
                W[i:j][abs(W[i:j]) >= t] -= t * np.sign(W[i:j][abs(W[i:j]) >= t])
                #print "soft threshold"
        
            elif shrink == "semisoft":
                W[i:j][abs(W[i:j]) >= t] = np.sqrt(W[i:j][abs(W[i:j]) >= t]**2 - t**2) * np.sign(W[i:j][abs(W[i:j]) >= t])

            i >>= 1
            j >>= 1
            
            c += 1
    
    return W
    
def psnr(truth, estimate):
    #returns Peak Signal to Noise Ratio in decibels
    N = len(truth)
    mean = np.mean(truth)
    #num = sum([abs(truth[i] - mean)**2 for i in range(N)])
    num = np.sum(truth ** 2)
    denom = np.sum((truth-estimate)**2)
    snr = 10 * np.log10(num/denom)
    
    return snr
    
def main():
    print "wavelets!"

def doppler(t):
    return np.sqrt(t*(t+1)) * np.sin(2*math.pi*1.05/(t+0.05))

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

def approx_plot(t, W, MAX_LAYER = -1, wv = "db4"):
    N = len(W)
    if MAX_LAYER == -1: 
        MAX_LAYER = (N).bit_length() - 1
    if wv == "db2":
        factor = 1
    else:
        factor = 1/math.sqrt(2)
    curr_factor = 1.0
    minlevel = N
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for level in range(MAX_LAYER, 0, -1):
        #use approximation
        ids = idwt(np.copy(W[:minlevel]), level, wv)
        ids *= curr_factor
        #print curr_factor
        space = int(N/minlevel)
        ts = t[space/2::space]
        name = "truth" if level == MAX_LAYER else "level" + str(level)
        ax.plot(ts, ids, label = name)
        
        curr_factor *= factor
        minlevel /= 2
 
    ax.legend(loc = "lower right")
    return fig

def coeff_plot(W, MAX_LAYER = -1):
    N = len(W)
    logN = (N).bit_length() - 1
    if MAX_LAYER == -1: MAX_LAYER = (N).bit_length() - 1
    
    fig = plt.figure()
    
    ax = fig.add_subplot(111)
    ax.plot(np.arange(N), W)
    wmax = np.max(W)
    wmin = np.min(W)
    divisions = 1 << np.arange(logN - MAX_LAYER, logN)
    ax.vlines(divisions, wmin, wmax)
    ax.set_xlim([0, N-1])
    ax.set_ylim([wmin,wmax])
    
    return fig
 
def coeff_pyramid(t, W, MAX_LAYER = -1, wavelet = "db4"):
    N = len(W)
    logN = (N).bit_length() - 1
    if MAX_LAYER == -1: MAX_LAYER = (N).bit_length() - 1
    
    fig, axs = plt.subplots(MAX_LAYER+1, 1, sharex = True)
    fig2, axs2 = plt.subplots(MAX_LAYER+1, 1, sharex = True, sharey = True)
    
    T = N >> (MAX_LAYER)
    Fs = 1 / (t[1]-t[0])
    
    currplot = 0
    max_freq = Fs / (1 << (MAX_LAYER - 1))
    
    for i in range(MAX_LAYER+1):
        #print T, 2*T
        spacing = N/T
        layer_t = t[spacing/2::spacing]
        #print layer_t.shape, layer_w.shape
        
        w_empty = np.zeros(N)
        
        if currplot == 0:
            title_str = "Scaling Coeff (0 - {:.2f} Hz)".format(Fs /(1 << (MAX_LAYER - 1)))
            layer_w = W[0:T]
            axs[i].set_title(title_str)
            axs2[i].set_title(title_str)
            w_empty[0:T] = W[0:T]
        
        else:
            #ax.set_title("Wavelet Coeff level " + str(currplot - 1) + " (" + str(max_freq) + " - " + str(max_freq*2) + " Hz)")
            title_str = "Wavelet Coeff level {:x} ({:.2f} - {:.2f} Hz)".format(currplot-1, max_freq, max_freq*2)
            axs[i].set_title(title_str)
            axs2[i].set_title(title_str)
            max_freq *= 2
            layer_w = W[T:2*T]
            w_empty[T:2*T] = W[T:2*T]
            T <<= 1
            
        currplot += 1
        
        
        axs[i].stem(layer_t, layer_w, basefmt = "black", markerfmt = " ")
        idw = idwt(np.copy(w_empty), MAX_LAYER, wavelet)
        axs2[i].plot(t, idw)
    return fig
 
def approx_test(wv = "db4"):
    N = 1024
    logN = (N).bit_length()-1
    t = np.linspace(0, 2, N)
    s = bumps(t)
    ML = 7
    ds = dwt(np.copy(s), ML, wv)

    ap = approx_plot(t, ds, ML, wv)
    cpyr = coeff_pyramid(t, ds, ML, wv)
    #cp = coeff_plot(ds, MAX_LAYER = ML)

def cascade_test():
    s2 = math.sqrt(2)
    s3 = math.sqrt(3)
    #hk = np.array([0.6830127, 1.1830127, 0.3169873, -0.1830127])/s2
    #hk = [1/s2, 1/s2]
    #hk = [(1+s3)/(4*s2), (3+s3)/(4*s2), (3-s3)/(4*s2), (1-s3)/(4*s2)]
    import scipy.signal
    hk = scipy.signal.daub(2)
    print hk
    (x, phi, psi) = scipy.signal.cascade(hk, 10)
    fphi = np.fft.fft(phi)
    fpsi = np.fft.fft(psi)
    
    
    f = np.linspace(-2/(x[1]-x[0]), 2/(x[1]-x[0]), len(x))
    #print f
    fig, (ax, ax2)= plt.subplots(1, 2)
    ax.plot(x, phi)
    ax.plot(x, psi)
    
    ax2.plot(f, np.abs(fphi))
    ax2.plot(f, np.abs(fpsi))
    plt.show()
    
def lift_test(wv = "db4"):
    #s = [1.0, 3.0, -2.0, 1.5, -0.5, 2.0]
    #s = np.array([32.0, 10.0, 20.0, 38.0, 37.0, 28.0, 38.0, 34.0, 18.0, 24.0, 18.0, 9.0, 23.0, 24.0, 28.0, 34.0])
    N = 2048
    t = np.linspace(0, 2, N)
    #s =  [math.sqrt(x/(2.0*N)*(x/(2.0*N)+1)) * math.sin(2*math.pi*1.05/(x/(2.0*N)+0.05))  for x in range(N)]
    s = heavisine(t)
    sigma = 0.50
    noise = np.random.normal(loc = 0.0, scale = sigma, size = N)
    rs = s + noise
    
    ML = 8
    
    ds = dwt(np.copy(s), MAX_LAYER = ML, wavelet = wv)
    rds = dwt(np.copy(rs), MAX_LAYER = ML, wavelet = wv)
    ids = idwt(np.copy(ds), MAX_LAYER = ML, wavelet = wv)
    #print "\nTruth"
    #waveletprint(ds)
    #waveletprint(threshold(ds, NUM_LAYERS = 1))
    #print "\nRandom:"
    #waveletprint(rds, MAX_LAYER = ML)
    
    #soft threshold
    srds = threshold(np.copy(rds), MAX_LAYER = ML, policy = "universal", shrink = "soft")
    isrds = idwt(np.copy(srds), MAX_LAYER = ML, wavelet = wv)
    #print "\nSoft:"
    #waveletprint(srds, MAX_LAYER = ML)
    
    #hard threshold
    hrds = threshold(np.copy(rds), MAX_LAYER = ML, policy = "universal", shrink = "hard")
    ihrds = idwt(np.copy(hrds), MAX_LAYER = ML, wavelet = wv)
    #print "\nHard:"
    #waveletprint(hrds, MAX_LAYER = ML)
    
    #hard universal threshold (level dependent)
    ulrds = threshold(np.copy(rds), MAX_LAYER = int(ML/2), policy = "universal_level", shrink = "hard")
    iulrds = idwt(np.copy(ulrds), MAX_LAYER = ML, wavelet = wv)
    
    #sure threshold
    surds = threshold(np.copy(rds), MAX_LAYER = int(ML/2), policy = "sure", shrink = "soft")
    isurds = idwt(np.copy(surds), MAX_LAYER = ML, wavelet = wv)
    
    tsnr = psnr(s, ids)
    #print "True PSNR:", tsnr, "dB"
    rsnr = psnr(s, rs)
    #print "Random PSNR:", rsnr, "dB"
    #nsnr = psnr(s, inrds)
    #print "Top 40 PSNR:", nsnr, "dB"
    ssnr = psnr(s, isrds)
    #print "Soft Universal Threshold PSNR:", ssnr, "dB"
    hsnr = psnr(s, ihrds)
    #print "Hard Universal Threshold PSNR:", hsnr, "dB"
    ulsnr = psnr(s, iulrds)
    #print "Hard Universal Level-Dependent PSNR:", ulsnr, "dB"
    susnr = psnr(s, isurds)
    #print "SURE Soft PSNR:", susnr, "dB"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t, s, label = "Ground Truth ({:.2f} dB)".format(tsnr))
    ax.plot(t, rs, label = "Noisy Data ({:.2f} dB)".format(rsnr))
    ax.plot(t, isrds, label = "Soft Threshold ({:.2f} dB)".format(ssnr))
    ax.plot(t, ihrds, label = "Hard Threshold ({:.2f} dB)".format(hsnr))
    ax.plot(t, iulrds, label = "Hard Universal Level ({:.2f} dB)".format(ulsnr))
    ax.plot(t, isurds, label = "SURE Soft ({:.2f} dB)".format(susnr))
    #ax.plot(t, inrds, label = "Top 40")
    ax.legend(loc="lower right")
    plt.show()
    
    #dwt and idwt are namespace sensitive, operate in place
    #use np.copy if copies of arrays are desired.
    
def haar_test():
    s = np.array([32.0, 10.0, 20.0, 38.0, 37.0, 28.0, 38.0, 34.0, 18.0, 24.0, 18.0, 9.0, 23.0, 24.0, 28.0, 34.0])
    #s = np.array([9, 7, -1, 1])
    N = len(s)
    #t = np.linspace(0, 2, N)
    print s
    ds = dwt(np.copy(s), wavelet = "haar")
    waveletprint(ds)
    ids = idwt(np.copy(ds), wavelet = "haar")
    print ids

    print "avg", np.average(s)
    
if __name__ == "__main__":
    #main()
    lift_test("db4")
    #approx_test("db2")
    approx_test("db4")
    #impulse_response_test()
    #cascade_test()
    #haar_test()
    plt.show()
    
