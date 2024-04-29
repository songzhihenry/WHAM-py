"""
# Created on Wed Aug  4 16:08:37 2021

# Implementation of the weight histogram analysis method The estimation of of
# the error bar, or the statistical uncertainty is based on
#
#  "Use of the weighted histogram analysis method for the analysis of
#  simulated and parallel tempering simulations", Chodera JD, Swope WC,
#  Pitera JW, Seok C, and Dill KA, JCTC, 2006
#
# Feng Ding <fding@unc.edu> Copyright 2011-2012

# @author: Zhiyuan Song <zhiyuas@g.clemson.edu>
"""
import argparse
import pandas as pd
import numpy as np
import math
#Import matplotlib.pyplot as plt
parser = argparse.ArgumentParser()
parser.add_argument('taskfile',help='Task file needed')
args = parser.parse_args()

#Read task file
df = pd.DataFrame([ line.strip('\n').split() for line in open('{}'.format(args.taskfile)) if not line.startswith('#') and not line.startswith('\n')])

#get key words
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
key_words_nu = []
key_words_str = []
for i in range(len(df)):
    if df.iloc[i,1] == 'PARAM':
        continue
    if is_number(df.iloc[i,1]):
        key_words_nu.append(df.iloc[i,0])
    else:
        key_words_str.append(df.iloc[i,0])

#Global params, values below will be changed after reading a task file
EPS = 1.0E-6;
NUM_BINS = 100;
NUM_SNAPS = 0;
MAX_STEPS = 100;
MAX_PARAM_STEPS = 40;
TOLERANCE = 0.0001;
MIN_TEMP = 270;
MAX_TEMP = 370;
TEMP_BINS = 100;
COMPUTE_DCV = 0;
COMPUTE_DPARAM = 0;

#Get params
for i in key_words_nu:
    if '.' in list(str(df[df[0] == i][1])):
        exec('{} = float(df[df[0] == i][1])'.format(i))
    else:
        exec('{} = int(df[df[0] == i][1])'.format(i))
for i in key_words_str:
    exec('{} = df.iloc[np.flatnonzero(df[0] == i)[0]].iat[1]'.format(i))
G_MK = np.zeros((NUM_BINS,NUM_REPLICA))
NE_ML = np.zeros((NUM_BINS,NUM_REPLICA))
HE_M = np.zeros(NUM_BINS)
H_M = np.zeros(NUM_BINS)

#Check validation
if DELTA_CV == 'YES':
    COMPUTE_DCV = 1
if DELTA_PARAM == 'YES':
    COMPUTE_DPARAM = 1

#Read the temperatures and check replica numbers
TEMPS = np.loadtxt('{}'.format(TEMPERATURES))
if len(TEMPS) != NUM_REPLICA:
    print ("Number of replica and temp mismatch")
    exit()
if TEMPS[0] != MIN_TEMP:
    print ("Min temp mismatch")
    exit()
if TEMPS[-1] != MAX_TEMP:
    print("Nax temp mismatch")

#Input for one dimensional PMF
try:
    PMF1D
except:
    print("No computation on PMF1D")
else:
    #temp_rw_stru = int(df.iloc[np.flatnonzero(df[0] == 'RESISTRU')[0]].iat[2])
    pmf1d_list = df.iloc[np.flatnonzero(df[0] == 'PMF1D')[0]].iat[1].split(',')
    pmf1d_bins = [int(x) for x in df.iloc[np.flatnonzero(df[0] == 'PMF1D')[0]].iat[2].split(',')]
    pmf1d_temps = [int(x) for x in df.iloc[np.flatnonzero(df[0] == 'PMF1D')[0]].iat[3].split(',')]
    if len(pmf1d_bins) != len(pmf1d_list):
        if len(pmf1d_bins) == 1:
            for i in range(len(pmf1d_list)-len(pmf1d_bins)):
                pmf1d_bins.append(pmf1d_bins[0])
        else:
            print("PMF1D param bins error input")
            exit()
    if len(pmf1d_temps) != len(pmf1d_list):
        if len(pmf1d_temps) == 1:
            for i in range(len(pmf1d_list)-len(pmf1d_temps)):
                pmf1d_temps.append(pmf1d_temps[0])
        else:
            print("PMF1D param temps error input")
            exit()

#Input for two dimensional PMF
try:
    PMF2D
except:
    print("No computation on PMF2D")
else:
    pmf2d_list = [x.split(',') for x in df.iloc[np.flatnonzero(df[0] == 'PMF2D')[0]].iat[1].split('/')]
    pmf2d_bins = [x.split(',') for x in df.iloc[np.flatnonzero(df[0] == 'PMF2D')[0]].iat[2].split('/')]
    pmf2d_bins = [list(map(int, x)) for x in pmf2d_bins]
    pmf2d_temps = [int(x) for x in df.iloc[np.flatnonzero(df[0] == 'PMF2D')[0]].iat[3].split('/')]
    if len(pmf2d_bins) != len(pmf2d_list):
        print("PMF2D param bins error input")
        exit()
    if len(pmf2d_temps) != len(pmf2d_list):
        print("PMF2D param temps error input")
        exit()

#Read param file list and input params
def read_flist():
    global FILE_LIST, NUM_REPLICA
    files_path = [line.strip('\n') for line in open('{}'.format(FILE_LIST))]
    param_trj = []
    for i in range(len(files_path)):
        param_tmp = []
        for j in range(NUM_REPLICA):
            param_tmp.append(np.loadtxt('{}.{}'.format(j,files_path[i])))
        param_tmp = np.array(param_tmp)    
        param_trj.append(param_tmp)
    return param_trj,files_path

# Convert strings in dssp files to numbers
def dssp_convert(filename):
    global RESISTRU
    A = open('{}.dssp'.format(filename))
    data = []
    lines = A.readlines()
    for line in lines:
        data.append(line.strip('\n'))
    A.close()
    data_t = np.zeros((4,len(data)))
    data_RH = np.zeros((len(data[0]),len(data)))
    data_RB = np.zeros((len(data[0]),len(data)))
    data_RT = np.zeros((len(data[0]),len(data)))
    data_RU = np.zeros((len(data[0]),len(data)))
    for i in range(len(data)):
        hang = list(data[i])
        #count
        helix = hang.count('H') + hang.count('G') + hang.count('I') 
        beta = hang.count('E') + hang.count('B')
        turn = hang.count('T') + hang.count('S')
        #unstr = len(hang) - helix - beta - turn
        unstr = hang.count('C')
        #percentage
        data_t[0][i] = helix/(len(hang))*100
        data_t[1][i] = beta/(len(hang))*100
        data_t[2][i] = turn/(len(hang))*100
        data_t[3][i] = unstr/(len(hang))*100
        if RESISTRU == 'YES':
            for j in range(len(hang)):
                if hang[j] == 'H' or hang[j] == 'G' or hang[j] == 'I':
                    data_RH[j][i] = 1
                if hang[j] == 'E' or hang[j] == 'B':
                    data_RB[j][i] = 1
                if hang[j] == 'T' or hang[j] == 'S':
                    data_RT[j][i] = 1
                if hang[j] == 'C':
                    data_RU[j][i] = 1
    if RESISTRU == 'YES':
        return data_t,data_RH,data_RB,data_RT,data_RU
    else:
        return data_t

def change_first_two_index(tri_matrix):
    list_and_array = []
    tri_matrix = np.array(tri_matrix).transpose(1,0,2)
    for i in range(len(tri_matrix)):
        list_and_array.append(tri_matrix[i])
    return list_and_array
# Read the file that include the extensions of large matrix 
def read_fmatrix(FILE_MATRIX):
    global NUM_REPLICA
    files_path = [line.strip('\n') for line in open('{}'.format(FILE_MATRIX))]
    param_tmps = []
    for i in range(len(files_path)):
        param_tmp = []
        for j in range(NUM_REPLICA):
            param_tmp.append(np.loadtxt('{}.{}'.format(j,files_path[i])).transpose()) #choose this when the file column number is equal to the number of frames
            #param_tmp.append(np.loadtxt('{}.{}'.format(j,files_path[i]))) #choose this when the file column number is equal to the number of parameters
            #param_tmp.append(np.load('{}.{}'.format(j,files_path[i])).transpose()) #choose this when the npy file column number is equal to the number of frames
            #param_tmp.append(np.load('{}.{}'.format(j,files_path[i]))) #choose this when the npy file column number is equal to the number of parameters
        param_tmps.append(change_first_two_index(param_tmp))    
    return param_tmps,files_path

# Dimension: param_trj[which param][which replica][which snap]
param_trj,param_list = read_flist()

#Get num of snaps
NUM_SNAPS = len(param_trj[0][0])

#Initialze the number of snaps for each temperature in a given replica simulation
NT_KL = np.ones((NUM_REPLICA,NUM_REPLICA))*NUM_SNAPS/NUM_REPLICA #the list of temperature fraction at different replica

#Find the max and min energy
pot_min = np.min(param_trj[0])
pot_max = np.max(param_trj[0])
pot_max = pot_max + abs(pot_max*0.00001)
dU = (pot_max-pot_min)/NUM_BINS
def compute_gmk():
    #stop summing ater the first negative
    global NUM_BINS, NUM_REPLICA, NUM_SNAPS, param_trj, dU, MAX_STEPS, G_MK, NE_ML, HE_M, H_M, pot_min
    idx = NUM_SNAPS
    G_MK = np.zeros((NUM_BINS,NUM_REPLICA))
    NE_ML = np.zeros((NUM_BINS,NUM_REPLICA))
    HE_M = np.zeros(NUM_BINS)
    H_M = np.zeros(NUM_BINS)
    reach_first_negative = np.zeros(NUM_BINS)
    for k in range(NUM_REPLICA):
        U = ((param_trj[0][k] - pot_min)/dU).astype(int) #% NUM_BINS
        #print(U)
        HM = (1e-100)*np.ones(NUM_BINS)
        for n in range(idx):
            HM[U[n]] = HM[U[n]] + 1
            H_M[U[n]] = H_M[U[n]] + 1
        #print (HM)
        f = HM/idx
        f = f*(1-f)
        #print (f)
        #compute correlation
        max_step = np.copy(MAX_STEPS)
        tmp_max = int(idx**0.5)
        #print(tmp_max)
        if max_step > tmp_max:
            max_step = tmp_max
        #print(max_step)
        for istep in range(max_step):
            step = int(1+(istep*(istep+1))/2)
            weight = istep + 1
            #initialize
            Sum = np.zeros(NUM_BINS)
            for i in range(idx-step):
                if U[i] == U[i+step]:
                    Sum [U[i]] = Sum [U[i]] + 1
            for m in range(NUM_BINS):
                corr = Sum[m]/(idx-step) - (HM[m]/idx)**2
                if reach_first_negative[m] == 0 and corr < 0:
                    reach_first_negative[m] = 1
                if reach_first_negative[m] == 0:
                    #print(k, step, m, corr, f[m],"\n")
                    if f[m] != 0:
                        corr = corr/f[m]
                    else:
                        corr = 0
                    G_MK[m][k] = G_MK[m][k] + (1.0-step/idx)*corr*weight
                #print(k,step,m,corr,G_MK[m][k],"\n")
        for m in range(NUM_BINS):
            G_MK[m][k] = 1 + 2*G_MK[m][k]
            HE_M[m] = HE_M[m] + HM[m]/G_MK[m][k]
            for l in range(NUM_REPLICA):
                NE_ML[m][l] = NE_ML[m][l]  + NT_KL[l][k]/G_MK[m][k]

#compute the statistical inefficiency g_mk, effective Histrogram count and snaps for each temp
compute_gmk()
#iterative calculation of Density of State (DOS)
BETA = 500/TEMPS
F_L = dU*BETA*NUM_BINS/3

tol = 1.0
#iteration steps
iter_N = 1000
LOG_OMIGA_M = np.zeros(NUM_BINS)
for i in range(iter_N):
    #compute LOG_OMIGA_M
    for m in range(NUM_BINS):
        exp_base = F_L[0] - BETA[0]*dU*(m+0.5)
        tmp = np.sum(NE_ML[m]*dU*np.exp(F_L - BETA*dU*(m+0.5)-exp_base))
        if HE_M[m] > 0 and tmp > 0:
            LOG_OMIGA_M[m] = math.log(HE_M[m]) - math.log(tmp) - exp_base
        else:
            print("%s, %s: empty bin\n" %(m, HE_M[m]))
            print("or divid by zero tmp")
            for l in range(NUM_REPLICA):
                print("%s %s: %s\n" %(m, l, NE_ML[m][l]), F_L[l]-BETA[l]*dU*(m+0.5)-exp_base)
            exit(0)
    #compute F_L
    has_zero = 0
    tmp_fl = np.zeros(NUM_REPLICA)
    for l in range(NUM_REPLICA):
        tmp = 0
        exp_base = LOG_OMIGA_M[0] - BETA[l]*dU*0.5
        for m in range(NUM_BINS):
            tmp = tmp + dU*math.exp(LOG_OMIGA_M[m]-BETA[l]*dU*(m+0.5)-exp_base)
        if tmp > 0:
            tmp_fl[l] = -math.log(tmp)-exp_base
        else:
            tmp_fl[l] = np.inf
            has_zero = 1;
    if has_zero == 1:
        tmp_sum0=0
        tmp_sum1=0
        tmp_count=0
        for l in range(NUM_REPLICA):
            if tmp_fl[l] != np.inf:
                tmp_sum0 = tmp_sum0 + F_L[l]
                tmp_sum1 = tmp_sum1 + tmp_fl[l]
                tmp_count = tmp_count + 1
        tmp_sum0 = tmp_sum0/tmp_count
        tmp_sum1 = tmp_sum1/tmp_count
        tmp_diff = tmp_sum1 - tmp_sum0
        for l in range(NUM_REPLICA):
            if tmp_fl[l] == np.inf:
                tmp_fl[l] = F_L[l] + tmp_diff
    diff = tmp_fl[0] - F_L[0]
    tol = 0
    F_L[0] = tmp_fl[0]
    for l in range(1,NUM_REPLICA):
        tol = tol + (tmp_fl[l]-F_L[l]-diff)**2
        F_L[l] = tmp_fl[l]
    #print ("iteration: %s tol: %s\n" %(i,tol));
    if tol < TOLERANCE:
        break
#compute the statistical uncertainty of OMIGA
H_M = np.where(np.isin(H_M, 0), 1, H_M)
dTEMP = (MAX_TEMP-MIN_TEMP)/TEMP_BINS
f = open('{}'.format(CV_OUT),'w')
for t in range(TEMP_BINS):
    #compute CV
    #1.use DOS
    temp = MIN_TEMP + dTEMP*(0.5+t)
    beta = 500/temp
    Z = 0
    EZ = 0
    EZ2 = 0
    base_exp = LOG_OMIGA_M[0]-beta*dU*0.5
    for m in range(NUM_BINS):
        e = pot_min + dU*(m+0.5)
        Z = Z + dU*math.exp(LOG_OMIGA_M[m]-beta*dU*(m+0.5)-base_exp)
        EZ = EZ + dU*math.exp(LOG_OMIGA_M[m]-beta*dU*(m+0.5)-base_exp)*e
        EZ2 = EZ2 + dU*math.exp(LOG_OMIGA_M[m]-beta*dU*(m+0.5)-base_exp)*(e**2)
    EZ = EZ/Z
    EZ2 = EZ2/Z
    cv = (beta**2)*(EZ2-EZ**2)/500
    #compute uncertainty of Cv
    #Cv = kB*beta^2*(<E^2>-<E>^2)
    #2. Use the actual energy data to comptue Cv (TESTED: similar result as abov)
    if COMPUTE_DCV == 1:
        DELTA_W2=0
        DELTA_X2=0
        DELTA_Y2=0
        DELTA_WY=0
        DELTA_XY=0
        DELTA_WX=0
        W = 0
        X = 0
        Y = 0
        y_kn = np.zeros((NUM_REPLICA,NUM_SNAPS))
        x_kn = np.zeros((NUM_REPLICA,NUM_SNAPS))
        w_kn = np.zeros((NUM_REPLICA,NUM_SNAPS))
        for k in range(NUM_REPLICA):
            for n in range(NUM_SNAPS):
                m = ((param_trj[0][k][n]-pot_min)/dU).astype(int) #% NUM_BINS
                y_kn[k][n] = math.exp(LOG_OMIGA_M[m]-beta*dU*(m+0.5)-base_exp)/H_M[m]          
            y_kn[k] = y_kn[k]/NUM_SNAPS
            x_kn[k] = y_kn[k]*(param_trj[0][k]**2)
            w_kn[k] = y_kn[k]*param_trj[0][k]
            ave_w = np.sum(w_kn[k])/NUM_SNAPS
            ave_x = np.sum(x_kn[k])/NUM_SNAPS
            ave_y = np.sum(y_kn[k])/NUM_SNAPS
            ave_w2 = np.sum(w_kn[k]**2)/NUM_SNAPS
            ave_x2 = np.sum(x_kn[k]**2)/NUM_SNAPS
            ave_y2 = np.sum(y_kn[k]**2)/NUM_SNAPS
            ave_wy = np.sum(w_kn[k]*y_kn[k])/NUM_SNAPS
            ave_xy = np.sum(x_kn[k]*y_kn[k])/NUM_SNAPS
            ave_wx = np.sum(w_kn[k]*x_kn[k])/NUM_SNAPS
            W = np.sum(w_kn)
            X = np.sum(x_kn)
            Y = np.sum(y_kn)
            sigma_ww = ave_w2 - ave_w**2
            sigma_xx = ave_x2 - ave_x**2
            sigma_yy = ave_y2 - ave_y**2
            sigma_wy = ave_wy - ave_w*ave_y
            sigma_wx = ave_wx - ave_w*ave_x
            sigma_xy = ave_xy - ave_x*ave_y
            if sigma_ww==0 or sigma_xx==0 or sigma_yy==0 or sigma_wx==0 or sigma_wy==0 or sigma_xy==0:
                continue
            #compute the statistical inefficiency
            max_step = MAX_PARAM_STEPS
            tmp_max = int(math.sqrt(NUM_SNAPS))
            if max_step > tmp_max:
                max_step = tmp_max
            gww=0
            gxx=0
            gyy=0
            gwy=0
            gxy=0
            gwx=0
            first_ww_negative = 0
            first_xx_negative = 0
            first_yy_negative = 0
            first_xy_negative = 0
            first_wy_negative = 0
            first_wx_negative = 0
            for istep in range(max_step):
                step = int(1 + (istep*(istep+1))/2)
                weight = istep + 1
                sum_ww = 0
                sum_xx = 0
                sum_yy = 0
                sum_wx = 0
                sum_wy = 0
                sum_xy = 0
                for n in range(NUM_SNAPS-step):
                    sum_ww = sum_ww + (w_kn[k][n]-ave_w)*(w_kn[k][n+step]-ave_w)
                    sum_xx = sum_xx + (x_kn[k][n]-ave_x)*(x_kn[k][n+step]-ave_x)
                    sum_yy = sum_yy + (y_kn[k][n]-ave_y)*(y_kn[k][n+step]-ave_y)
                    sum_xy = sum_xy + ((x_kn[k][n]-ave_x)*(y_kn[k][n+step]-ave_y) + (y_kn[k][n]-ave_y)*(x_kn[k][n+step]-ave_x))
                    sum_wy = sum_wy + ((w_kn[k][n]-ave_w)*(y_kn[k][n+step]-ave_y) + (y_kn[k][n]-ave_y)*(w_kn[k][n+step]-ave_w))
                    sum_wx = sum_wx + ((w_kn[k][n]-ave_w)*(x_kn[k][n+step]-ave_x) + (x_kn[k][n]-ave_x)*(w_kn[k][n+step]-ave_w))
                sum_ww = sum_ww/(NUM_SNAPS-step)
                sum_xx = sum_xx/(NUM_SNAPS-step)
                sum_yy = sum_yy/(NUM_SNAPS-step)
                sum_wx = sum_wx/2/(NUM_SNAPS-step)
                sum_wy = sum_wy/2/(NUM_SNAPS-step)
                sum_xy = sum_xy/2/(NUM_SNAPS-step)
                corr_ww = sum_ww/sigma_ww
                corr_xx = sum_xx/sigma_xx
                corr_yy = sum_yy/sigma_yy
                corr_wx = sum_wx/sigma_wx
                corr_wy = sum_wy/sigma_wy
                corr_xy = sum_xy/sigma_xy
                if first_ww_negative == 0 and corr_ww<EPS:
                    first_ww_negative = 1
                if first_xx_negative == 0 and corr_xx<EPS:
                    first_xx_negative = 1
                if first_yy_negative == 0 and corr_yy<EPS:
                    first_yy_negative = 1
                if first_xy_negative == 0 and corr_xy<EPS:
                    first_xy_negative = 1
                if first_wx_negative == 0 and corr_wx<EPS:
                    first_wx_negative = 1
                if first_wy_negative == 0 and corr_wy<EPS:
                    first_wy_negative = 1
                if first_ww_negative == 0:
                    gww = gww + (1-step/NUM_SNAPS)*corr_ww*weight
                if first_xx_negative == 0:
                    gxx = gxx + (1-step/NUM_SNAPS)*corr_xx*weight
                if first_yy_negative == 0:
                    gyy = gyy + (1-step/NUM_SNAPS)*corr_yy*weight
                if first_wx_negative == 0:
                    gwx = gwx + (1-step/NUM_SNAPS)*corr_wx*weight
                if first_wy_negative == 0:
                    gwy = gwy + (1-step/NUM_SNAPS)*corr_wy*weight
                if first_xy_negative == 0:
                    gxy = gxy + (1-step/NUM_SNAPS)*corr_xy*weight
                if first_ww_negative == 1 and first_xx_negative == 1 and first_yy_negative == 1 and first_wx_negative == 1 and first_wy_negative == 1 and first_xy_negative == 1:
                    break
            gww = 1 + 2*gww
            gwx = 1 + 2*gwx
            gwy = 1 + 2*gwy
            gxx = 1 + 2*gxx
            gyy = 1 + 2*gyy
            gxy = 1 + 2*gxy
            DELTA_W2 = DELTA_W2 + gww*sigma_ww*NUM_SNAPS
            DELTA_X2 = DELTA_X2 + gxx*sigma_xx*NUM_SNAPS
            DELTA_Y2 = DELTA_Y2 + gyy*sigma_yy*NUM_SNAPS
            DELTA_WX = DELTA_WX + gwx*sigma_wx*NUM_SNAPS
            DELTA_WY = DELTA_WY + gwy*sigma_wy*NUM_SNAPS
            DELTA_XY = DELTA_XY + gxy*sigma_xy*NUM_SNAPS
        aveE2 = X/Y
        aveE = W/Y
        delE = (aveE**2)*(DELTA_W2/(W**2) + DELTA_Y2/(Y**2) - 2*DELTA_WY/(W*Y))
        delE2 = (aveE2**2)*(DELTA_X2/(X**2) + DELTA_Y2/(Y**2) - 2*DELTA_XY/(X*Y))
        dE_dE2 = aveE*aveE2*(DELTA_WX/(X*W) + DELTA_Y2/(Y**2) - DELTA_XY/(X*Y) - DELTA_WY/(W*Y))
        delta = delE2 + 4*(aveE**2)*delE - 4*aveE*dE_dE2
        the_cv = (beta**2)*(aveE2-aveE**2)/500
        dCV = (beta**2)*math.sqrt(abs(delta))/500
        print("%s   %s   %s   %s" %(temp,cv,the_cv,dCV),file = f)
        #print("%s   %s   %s   %s" %(temp,cv,the_cv,dCV))
    else:
        print("%s   %s"%(temp,cv), file = f)
        #print("%s   %s"%(temp,cv))
f.close()

#Compute the average parameters and standards
def ave_standa_param(param_trj):
    global NUM_BINS, TEMP_BINS, LOG_OMIGA_M, MIN_TEMP, dTEMP, COMPUTE_DPARAM, dU, H_M, pot_min
    NUM_PARAMS = len(param_trj)-1
    NUM_SNAPS = len(param_trj[0][0])
    if COMPUTE_DPARAM == 1:
        output = np.zeros((NUM_PARAMS,TEMP_BINS,3))
    else:
        output = np.zeros((NUM_PARAMS,TEMP_BINS,2))
    
    for t in range(TEMP_BINS):
        #computer param
        temp = MIN_TEMP + dTEMP*(0.5 + t)
        beta = 500/temp
        base_exp = LOG_OMIGA_M[0]-beta*dU*0.5
        Z = 0
        ave_param = np.zeros(NUM_PARAMS)
        fluc_param = np.zeros(NUM_PARAMS)
        DELTA_Y2 = 0
        DELTA_X2 = np.zeros(NUM_PARAMS)
        DELTA_XY = np.zeros(NUM_PARAMS)
        Y = 0
        X = np.zeros(NUM_PARAMS)
        for k in range(NUM_REPLICA):
            for n in range(NUM_SNAPS):
                m = ((param_trj[0][k][n]-pot_min)/dU).astype(int) #% NUM_BINS
                y_kn = np.exp(LOG_OMIGA_M[m]-beta*dU*(m+0.5)-base_exp)/H_M[m]
                Z = Z + y_kn
                for ip in range(NUM_PARAMS):
                    ave_param[ip] = ave_param[ip] + y_kn*param_trj[ip+1][k][n]
                    fluc_param[ip] = fluc_param[ip] + y_kn*param_trj[ip+1][k][n]**2
            if COMPUTE_DPARAM == 1:
                X_N = np.zeros((NUM_PARAMS,NUM_SNAPS))
                Y_N = np.zeros(NUM_PARAMS)
                ave_y = 0
                ave_y2 = 0
                ave_x = np.zeros(NUM_PARAMS)
                ave_x2 = np.zeros(NUM_PARAMS)
                ave_xy = np.zeros(NUM_PARAMS)
                for n in range(NUM_SNAPS):
                    m = (param_trj[0][k][n]).astype(int) #% NUM_BINS
                    Y_N[n] = np.exp(LOG_OMIGA_M[m]-beta*dU*(m+0.5)-base_exp)/H_M[m]/NUM_BINS
                    Y = Y + Y_N[n]
                    ave_y = ave_y + Y_N[n]
                    ave_y2 = ave_y2 + Y_N[n]**2
                    for p in range(NUM_PARAMS):
                        X_N[p][n] = Y_N[n]*param_trj[p+1][k][n]
                        ave_x[p] = ave_x[p] + X_N[p][n]
                        ave_x2[p] = ave_x2[p] + X_N[p][n]**2
                        ave_xy[p] = ave_xy[p] + X_N[p][n]*Y_N[n]
                        X[p] = X[p] + X_N[p][n]
                ave_y = ave_y/NUM_SNAPS
                ave_y2 = ave_y2/NUM_SNAPS
                sigma_yy = ave_y2 - ave_y**2
                if sigma_yy == 0:
                    continue
                sigma_xx = np.zeros(NUM_PARAMS)
                sigma_xy = np.zeros(NUM_PARAMS)
                ave_x = ave_x/NUM_SNAPS
                ave_x2 = ave_x2/NUM_SNAPS
                ave_xy = ave_xy/NUM_SNAPS
                sigma_xx = ave_x2 - ave_x**2
                sigma_xy = ave_xy - ave_x*ave_y
                if 0 in sigma_xx or 0 in sigma_xy:
                    continue
                max_step = MAX_PARAM_STEPS
                tmp_max = int(math.sqrt(NUM_SNAPS))
                if max_step > tmp_max:
                    max_step = tmp_max
                gxx = np.zeros(NUM_PARAMS)
                gxy = np.zeros(NUM_PARAMS)
                gyy = 0
                first_yy_negative = 0
                first_xx_negative = np.zeros(NUM_PARAMS)
                first_xy_negative = np.zeros(NUM_PARAMS)
                for istep in range(max_step):
                    step = 1+(istep*(istep+1))/2
                    weight = istep + 1
                    sum_yy = 0
                    sum_xx = np.zeros(NUM_PARAMS)
                    sum_xy = np.zeros(NUM_PARAMS)
                    for n in range(NUM_SNAPS):
                        sum_yy = sum_yy + (Y_N[n]-ave_y)*(Y_N[n+step]-ave_y)
                        for p in range(NUM_PARAMS):
                            sum_xx[p] = sum_xx[p] + (X_N[p][n]-ave_x[p])*(X_N[p][n+step]-ave_x[p])
                            sum_xy[p] = sum_xy[p] + (X_N[p][n]-ave_x[p])*(Y_N[n+step]-ave_y) + (X_N[p][n+step]-ave_x[p])*(Y_N[n]-ave_y)
                    sum_yy = sum_yy/(NUM_SNAPS - step)
                    corr_yy = sum_yy/sigma_yy
                    if first_yy_negative == 0 and corr_yy < EPS:
                        first_yy_negative = 1
                    corr_xx = np.zeros(NUM_PARAMS)
                    corr_xy = np.zeros(NUM_PARAMS)
                    sum_xx = sum_xx/(NUM_SNAPS - step)
                    sum_xy = sum_xx/(NUM_SNAPS - step)
                    corr_xx = sum_xx/sigma_xx
                    corr_xy = sum_xx/sigma_xy
                    for p in range(NUM_PARAMS):
                        if first_xy_negative[p] == 0 and first_xy_negative[p] < EPS:
                            first_xy_negative[p] = 1
                        if first_xx_negative[p] == 0 and first_xx_negative[p] < EPS:
                            first_xx_negative[p] = 1
                    if first_yy_negative == 0:
                        gyy = gyy + (1-step/NUM_SNAPS)*corr_yy*weight
                        first_negative = first_yy_negative
                    for p in range(NUM_PARAMS):
                        if first_xx_negative[p] == 0:
                            gxx[p] = gxx[p] + (1-step/NUM_SNAPS)*corr_xx[p]*weight
                        if first_xy_negative == 0:
                            gxy[p] = gxy[p] + (1-step/NUM_SNAPS)*corr_xy[p]*weight
                        first_negative = first_negative*first_xx_negative[p]
                        first_negative = first_negative*first_xy_negative[p]
                        if first_negative ==1:
                            break
                gyy = 1 + 2*gyy
                DELTA_Y2 = DELTA_Y2 + gyy*sigma_yy*NUM_SNAPS
                gxx = 1 + 2*gxx
                gxy = 1 + 2*gxy
                DELTA_X2 = DELTA_X2 + gxx*sigma_xx*NUM_SNAPS
                DELTA_XY = DELTA_XY + gxy*sigma_xy*NUM_SNAPS
        for i in range(NUM_PARAMS):
            ave = ave_param[i]/Z
            fluc = fluc_param[i]/Z
            output[i][t][0] = ave
            output[i][t][1] = fluc - ave**2
            if COMPUTE_DPARAM == 1:
                delta = (ave**2)*(DELTA_X2[p]/(X[p]*X[p])+DELTA_Y2/(Y*Y)-2*DELTA_XY[p]/(Y*X[p]))
                output[i][t][2] = math.sqrt(delta)
    return output
# Computa average of parameters at a certain temperature
def ave_param_at_temp(param_trj,temp):
    global NUM_BINS, TEMP_BINS, LOG_OMIGA_M, MIN_TEMP, dTEMP, COMPUTE_DPARAM, dU, H_M, pot_min
    NUM_PARAMS = len(param_trj)-1
    NUM_SNAPS = len(param_trj[0][0])
    output = np.zeros((NUM_PARAMS,2))
        #computer param
    beta = 500/temp
    base_exp = LOG_OMIGA_M[0]-beta*dU*0.5
    Z = 0
    ave_param = np.zeros(NUM_PARAMS)
    fluc_param = np.zeros(NUM_PARAMS)
    for k in range(NUM_REPLICA):
        for n in range(NUM_SNAPS):
            m = ((param_trj[0][k][n]-pot_min)/dU).astype(int) % NUM_BINS
            y_kn = np.exp(LOG_OMIGA_M[m]-beta*dU*(m+0.5)-base_exp)/H_M[m]
            Z = Z + y_kn
            for ip in range(NUM_PARAMS):
                ave_param[ip] = ave_param[ip] + y_kn*param_trj[ip+1][k][n]
                fluc_param[ip] = fluc_param[ip] + y_kn*param_trj[ip+1][k][n]**2
    for i in range(NUM_PARAMS):
        ave = ave_param[i]/Z
        fluc = fluc_param[i]/Z
        output[i][0] = ave
        output[i][1] = fluc - ave**2
    return output

#Output to a file
def coutput(output,param_list):
    global COMPUTE_DPARAM, MIN_TEMP, dTEMP
    NUM_PARAMS = len(param_list) - 1
    for i in range(NUM_PARAMS):
        f = open('ave_{}.dat'.format(param_list[i+1]),'w')
        for n in range(TEMP_BINS):
            temp = MIN_TEMP + dTEMP*(0.5+n)
            if COMPUTE_DPARAM == 1:
                print ("%s   %s   %s   %s" %(temp,output[i][n][0],output[i][n][1],output[i][n][2]),file=f)
            else:
                print("%s   %s   %s" %(temp,output[i][n][0],output[i][n][1]),file=f)
        f.close()

#Check dssp files and output calculations
if DSSP == 'YES':
    data_t = []
    sec_st_list = ['helix','betasheet','turn','unstr']
    if RESISTRU == 'NO':
        for j in range(NUM_REPLICA):
            data_t.append(dssp_convert(j))
    else:
        data_RB = []
        data_RH = []
        data_RT = []
        data_RU = []
        temp_rw_stru = [int(x) for x in df.iloc[np.flatnonzero(df[0] == 'RESISTRU')[0]].iat[2].split(',')]
        for j in range(NUM_REPLICA):
            chain4,RH,RB,RT,RU = dssp_convert(j)
            data_t.append(chain4)
            data_RH.append(RH)
            data_RB.append(RB)
            data_RT.append(RT)
            data_RU.append(RU)
        sec_RW_list = [x for x in range(1,len(data_RH[0])+1)]
        data_RH = change_first_two_index(data_RH);data_RH.insert(0,param_trj[0])
        data_RB = change_first_two_index(data_RB);data_RB.insert(0,param_trj[0])
        data_RT = change_first_two_index(data_RT);data_RT.insert(0,param_trj[0])
        data_RU = change_first_two_index(data_RU);data_RU.insert(0,param_trj[0])
        del RH,RB,RT,RU
    data_t =  change_first_two_index(data_t)
    param_list = param_list + sec_st_list
    param_trj = param_trj + data_t
    del sec_st_list,data_t

#Compute 1d pmf
def compute_pmf1d(name, bins, temp):
    global NUM_BINS, LOG_OMIGA_M, du, pot_min, param_trj, param_list, NUM_REPLICA, NUM_SNAPS, H_M
    beta = 500/temp
    base_exp = LOG_OMIGA_M[0]-beta*dU*0.5
    ip = param_list.index(name)
    param_min = np.min(param_trj[ip])
    param_max = np.max(param_trj[ip])
    param_max = param_max + abs(0.00001*param_max)
    dparam = (param_max - param_min)/bins
    Z = 0
    pmfout = np.zeros((bins,2))
    N = []
    for k in range(NUM_REPLICA):
        for idx in range(NUM_SNAPS):
            m = ((param_trj[0][k][idx]-pot_min)/dU).astype(int) #% NUM_BINS
            n = ((param_trj[ip][k][idx]-param_min)/dparam).astype(int) #% bins
            y_kn = math.exp(LOG_OMIGA_M[m]-beta*dU*(m+0.5)-base_exp)/H_M[m]
            Z = Z + y_kn
            pmfout[n][1] = pmfout[n][1] + y_kn
            N.append(n)
    N = list(set(N))
    pmfout_tmp = np.array([pmfout[x][1] for x in N])
    pmfout_tmp = -np.log(pmfout_tmp/Z)/beta
    for k in range(bins):
        pmfout[k][0] = param_min + (k+0.5)*dparam
        if k in N:
            pmfout[k][1] = -math.log(pmfout[k][1]/Z)/beta - np.min(pmfout_tmp)
        else:
            pmfout[k][1] = np.max(pmfout_tmp) - np.min(pmfout_tmp) + 10
    return pmfout

#Compute 2d pmf
def compute_pmf2d(names,bins,temps):
    global NUM_BINS, LOG_OMIGA_M, du, pot_min, param_trj, param_list, NUM_REPLICA, NUM_SNAPS, H_M, pot_max
    beta = 500/temps
    base_exp = LOG_OMIGA_M[0]-beta*dU*0.5
    Z = 0
    ip = param_list.index(names[0])
    jp = param_list.index(names[1])
    param1_min = pot_min
    param1_max = pot_max
    param2_min = pot_min
    param2_max = pot_max
    dparam1 = dU
    dparam2 = dU
    if ip > 0:
        param1_min = np.min(param_trj[ip])
        param1_max = np.max(param_trj[ip])
        param1_max = param1_max + abs(0.00001*param1_max)
        dparam1 = (param1_max - param1_min)/bins[0]
    if jp > 0:
        param2_min = np.min(param_trj[jp])
        param2_max = np.max(param_trj[jp])
        param2_max = param2_max + abs(0.00001*param2_max)
        dparam2 = (param2_max - param2_min)/bins[1]   
    prob_param = np.zeros((bins[0],bins[1]))
    for k in range(NUM_REPLICA):
        for idx in range(NUM_SNAPS):
            m = ((param_trj[0][k][idx]-pot_min)/dU).astype(int) #% NUM_BINS
            y_kn = math.exp(LOG_OMIGA_M[m]-beta*dU*(m+0.5)-base_exp)/H_M[m]
            Z = Z + y_kn
            n1 = m
            n2 = m
            if ip > 0:
                n1 = ((param_trj[ip][k][idx]-param1_min)/dparam1).astype(int) #% bins[0]
            if jp > 0:
                n2 = ((param_trj[jp][k][idx]-param2_min)/dparam2).astype(int) #% bins[1]          
            prob_param[n1][n2] = prob_param[n1][n2] + y_kn
    prob_max = np.min(prob_param[np.nonzero(prob_param)])
    prob_param = np.where(np.isin(prob_param, 0), prob_max, prob_param)
    prob_param = -np.log(prob_param/Z)/beta
    pmfout = np.zeros((bins[0]*bins[1],3))
    for k in range(bins[0]):
        for l in range(bins[1]):
            pmfout[k*bins[1]+l][0] = param1_min + (k+0.5)*dparam1
            pmfout[k*bins[1]+l][1] = param2_min + (l+0.5)*dparam2
            pmfout[k*bins[1]+l][2] = prob_param[k][l]
    return pmfout

#Output 2ndst calculation
try:
    data_RH
    data_RB
    data_RT
    data_RU
except:
    print("No 2ndst calculation for each residue")
else:
    for i in temp_rw_stru:
        np.savetxt('RW_HE_{}.txt'.format(i),ave_param_at_temp(data_RH,i))
        np.savetxt('RW_BT_{}.txt'.format(i),ave_param_at_temp(data_RB,i))
        np.savetxt('RW_TN_{}.txt'.format(i),ave_param_at_temp(data_RT,i))
        np.savetxt('RW_UT_{}.txt'.format(i),ave_param_at_temp(data_RU,i))

#Check any param matrix and out put average calculation
try:
    open('{}'.format(FILE_MATRIX))
except:
    print('No matrix input')
else:
    data_params,param_names = read_fmatrix(FILE_MATRIX)
    for i in range(len(param_names)):
        data_params[i].insert(0,param_trj[0])
    temp_rw_param = [int(x) for x in df.iloc[np.flatnonzero(df[0] == 'FILE_MATRIX')[0]].iat[2].split(',')]
    #output ave_param
    for i in temp_rw_param:
        for j in range(len(param_names)):
            np.savetxt('RW_{}_{}.txt'.format(param_names[j],i),ave_param_at_temp(data_params[j],i))

#Output regular ave_params
if AVE_PARAM == 'YES':
    coutput(ave_standa_param(param_trj),param_list)

#Output PMF 1D/2D computation
try:
    PMF1D
except:
    print("No PMF1D Computation")
else:
    for i in range(len(pmf1d_list)):
        np.savetxt('pmf1d_{}_{}.dat'.format(pmf1d_list[i],pmf1d_temps[i]),compute_pmf1d(pmf1d_list[i],pmf1d_bins[i],pmf1d_temps[i]))
try:
    PMF2D
except:
    print("No PMF2d Computation")
else:
    for i in range(len(pmf2d_list)):
        np.savetxt('pmf2d_{}_{}_{}.dat'.format(pmf2d_list[i][0],pmf2d_list[i][1],pmf2d_temps[i]),
                   compute_pmf2d(pmf2d_list[i],pmf2d_bins[i],pmf2d_temps[i]))
