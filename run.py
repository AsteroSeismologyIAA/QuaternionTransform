import numpy as np
import math
from tqdm import tqdm
import multiprocessing as mp


def fNinner(N,J,am,ji,dpi,rci,cci,rmi,rci1,rci2,rci3):
    
    for I in range(0, N):
        #th1 = dmod(I * 1.0 * (J + ji) / n, 1.0)  # Assuming dmod function is defined
        # I+1 and J+1 due to Fortran Indexes
        
        th1 = ((I+1) * 1.0 * ((J+1) + ji) / N)%1  # Fractional part
        gr = th1 * dpi
        c1 = math.cos(gr)  # Assuming dcos function is defined
        s1 = math.sin(gr)  # Assuming dsin function is defined

        ff = math.exp(math.asinh(th1 / 2.0))  # Assuming exp and asinh functions are defined
        wo = dpi * (ff - 1.0 / ff) / 2
        wp = dpi * (ff + 1.0 / ff) / 2

        z1 = math.cos(wp) * math.cos(wo)
        z2 = math.sin(wp) * math.cos(wo)
        z3 = math.cos(wp) * math.sin(wo)
        z4 = math.sin(wp) * math.sin(wo)

        rci += am[I] * c1
        cci += am[I] * s1
        rmi += am[I] * z1
        rci2 += am[I] * z2
        rci3 += am[I] * z3
        rci1 += am[I] * z4      
        
    #return rci, cci, rmi, rci2, rci3, rci1
    return rci,cci,rmi,rci1,rci2,rci3

print("Number of processors: ", mp.cpu_count())



def quaternion(input_file,anuino,anufio,rti,tpas,output_file):

    dpi = np.arctan(1.0) * 8.0
    wo = wp = th = gr = st = rci1 = rci2 = rci3 = rmi = 0.0

    AM = np.zeros(8388608, dtype=np.float64)
    
    print("file, initial frequency, final frequency, days, sampling(sec)",(input_file,anuino,anufio,rti,tpas))


    n = rti * 86400.0 / (1 * tpas)
    nw = math.log(n * 1.0) / math.log(2.0)


    am = np.zeros(int(n))
    AMEDI = 0


    with open(input_file, 'r') as file:
        for I in range(0, int(n)):
            line = file.readline().split()
            if not line:
                pass
            else:
                am[I] = float(line[0])
                AMEDI += am[I]

    print (AMEDI, am,len(am),int(n))



    N = int(n)
    AMEDI = AMEDI / N

    sig = 0.0

    for I in range(0, N):
        sig += (am[I] - AMEDI)**2

    sig = math.sqrt(sig / N)

    #print(sig,AMEDI,N)

    fny = 1000000.0 / (2 * tpas)
    delnu = 1.0e6 / (n * tpas)
    ji = math.floor(anuino / delnu)
    anuin = ji * delnu
    jf = math.floor(anufio / delnu)
    NP = jf - ji

    #print('Frequencies, data, days, Ny')
    #print(NP, N, N * tpas / 86400.0, fny)

    F = anuin
    BMAX = 0.0

    it = 0
    print("Iterations: ",NP,"x",N,"=",NP*N)
    for J in tqdm(range(0, NP),desc="NP",leave=False):
        rci = 0.0
        cci = 0.0
        rmi = 0.0
        rci1 = 0.0
        rci2 = 0.0
        rci3 = 0.0

        rci,cci,rmi,rci1,rci2,rci3 = fNinner(N,J,am,ji,dpi,rci,cci,rmi,rci1,rci2,rci3)
        
        rci = rci * 2.0 / n
        cci = cci * 2.0 / n
        rmi = rmi * 2.0 / n
        rci2 = rci2 * 2.0 / n
        rci3 = rci3 * 2.0 / n
        rci1 = rci1 * 2.0 / n

        b = rci**2 + cci**2
        Bg = rmi**2 + (rci2**2 + rci1**2 + rci3**2)

        F = anuin + (J - 1) * delnu
        if Bg > BMAX:
            FMAX = F
            BMAX = Bg
            jmax = J

        FDAT[J] = float(F)
        TRADAT[J] = float(b)
        TRADATg[J] = float(Bg)
        print (F,b,Bg)
    
    with open(output_file, 'w') as file_out:
        for I in range(0, NP):
            file_out.write(f"{FDAT[I]} {TRADAT[I]} {TRADATg[I]}\n")

    return FDAT,TRADAT,TRADATg


TRADAT = np.zeros(8388608, dtype=np.float64)
FDAT = np.zeros(8388608, dtype=np.float64)
TRADATg = np.zeros(8388608, dtype=np.float64)

if __name__== "__main__":

    output_file = 'ps.dgt'
    input_file = "../golfsel1.dat"
    anuino = 0
    anufio = 25000
    rti = 16
    tpas = 20

    FDAT,TRADAT,TRADATg = quaternion(input_file,anuino,anufio,rti,tpas,output_file)
    