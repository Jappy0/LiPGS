#!/usr/bin/env python

'''
This script merges all scripts related to PRS-CS into one file.
'''

"""
PRS-CS: a polygenic prediction method that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
using GWAS summary statistics and an external LD reference panel.

Reference: T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors.
           Nature Communications, 10:1776, 2019.


Usage:
python PRScs.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --out_dir=OUTPUT_DIR
                [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR
                 --chrom=CHROM --write_psi=WRITE_PSI --write_pst=WRITE_POSTERIOR_SAMPLES --seed=SEED]

"""


import os
import sys
import getopt
import numpy as np
from scipy.stats import norm
from scipy import linalg
import h5py
import math
from numpy import random
import numpy as np
from scipy import linalg 
from numpy import random

def psi(x, alpha, lam):
    f = -alpha*(math.cosh(x)-1.0)-lam*(math.exp(x)-x-1.0)
    return f

def dpsi(x, alpha, lam):
    f = -alpha*math.sinh(x)-lam*(math.exp(x)-1.0)
    return f

def g(x, sd, td, f1, f2):
    if (x >= -sd) and (x <= td):
        f = 1.0
    elif x > td:
        f = f1
    elif x < -sd:
        f = f2

    return f

def gigrnd(p, a, b):
    # setup -- sample from the two-parameter version gig(lam,omega)
    p = float(p); a = float(a); b = float(b)
    lam = p
    omega = math.sqrt(a*b)

    if lam < 0:
        lam = -lam
        swap = True
    else:
        swap = False

    alpha = math.sqrt(math.pow(omega,2)+math.pow(lam,2))-lam

    # find t
    x = -psi(1.0, alpha, lam)
    if (x >= 0.5) and (x <= 2.0):
        t = 1.0
    elif x > 2.0:
        if (alpha == 0) and (lam == 0):
            t = 1.0
        else:
            t = math.sqrt(2.0/(alpha+lam))
    elif x < 0.5:
        if (alpha == 0) and (lam == 0):
            t = 1.0
        else:
            t = math.log(4.0/(alpha+2.0*lam))

    # find s
    x = -psi(-1.0, alpha, lam)
    if (x >= 0.5) and (x <= 2.0):
        s = 1.0
    elif x > 2.0:
        if (alpha == 0) and (lam == 0):
            s = 1.0
        else:
            s = math.sqrt(4.0/(alpha*math.cosh(1)+lam))
    elif x < 0.5:
        if (alpha == 0) and (lam == 0):
            s = 1.0
        elif alpha == 0:
            s = 1.0/lam
        elif lam == 0:
            s = math.log(1.0+1.0/alpha+math.sqrt(1.0/math.pow(alpha,2)+2.0/alpha))
        else:
            s = min(1.0/lam, math.log(1.0+1.0/alpha+math.sqrt(1.0/math.pow(alpha,2)+2.0/alpha)))

    # find auxiliary parameters
    eta = -psi(t, alpha, lam)
    zeta = -dpsi(t, alpha, lam)
    theta = -psi(-s, alpha, lam)
    xi = dpsi(-s, alpha, lam)

    p = 1.0/xi
    r = 1.0/zeta

    td = t-r*eta
    sd = s-p*theta
    q = td+sd

    # random variate generation
    while True:
        U = random.random()
        V = random.random()
        W = random.random()
        if U < q/(p+q+r):
            rnd = -sd+q*V
        elif U < (q+r)/(p+q+r):
            rnd = td-r*math.log(V)
        else:
            rnd = -sd+p*math.log(V)

        f1 = math.exp(-eta-zeta*(rnd-t))
        f2 = math.exp(-theta+xi*(rnd+s))
        if W*g(rnd, sd, td, f1, f2) <= math.exp(psi(rnd, alpha, lam)):
            break

    # transform back to the three-parameter version gig(p,a,b)
    rnd = math.exp(rnd)*(lam/omega+math.sqrt(1.0+math.pow(lam,2)/math.pow(omega,2)))
    if swap:
        rnd = 1.0/rnd

    rnd = rnd/math.sqrt(a/b)
    return rnd

def mcmc(a, b, phi, sst_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, out_dir, beta_std, write_psi, write_pst, seed):
    print('... MCMC ...')

    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    beta_mrg = np.array(sst_dict['BETA'], ndmin=2).T
    maf = np.array(sst_dict['MAF'], ndmin=2).T
    n_pst = int((n_iter-n_burnin)/thin)
    p = len(sst_dict['SNP'])
    n_blk = len(ld_blk)

    # initialization
    beta = np.zeros((p,1))
    psi = np.ones((p,1))
    sigma = 1.0
    
    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    if write_pst == 'TRUE':
        beta_pst = np.zeros((p,n_pst))

    beta_est = np.zeros((p,1))
    psi_est = np.zeros((p,1))
    sigma_est = 0.0
    phi_est = 0.0

    # MCMC
    pp = 0
    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('--- iter-' + str(itr) + ' ---')

        mm = 0; quad = 0.0
        for kk in range(n_blk):
            if blk_size[kk] == 0:
                continue
            else:
                idx_blk = range(mm,mm+blk_size[kk])
                dinvt = ld_blk[kk]+np.diag(1.0/psi[idx_blk].T[0])
                dinvt_chol = linalg.cholesky(dinvt)
                beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[idx_blk], trans='T') + np.sqrt(sigma/n)*random.randn(len(idx_blk),1)
                beta[idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                quad += np.dot(np.dot(beta[idx_blk].T, dinvt), beta[idx_blk])
                mm += blk_size[kk]

        err = max(n/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n/2.0*sum(beta**2/psi))
        sigma = 1.0/random.gamma((n+p)/2.0, 1.0/err)

        delta = random.gamma(a+b, 1.0/(psi+phi))

        for jj in range(p):
            psi[jj] = gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)
        psi[psi>1] = 1.0

        if phi_updt == True:
            w = random.gamma(1.0, 1.0/(phi+1.0))
            phi = random.gamma(p*b+0.5, 1.0/(sum(delta)+w))

        # posterior
        if (itr>n_burnin) and (itr % thin == 0):
            beta_est = beta_est + beta/n_pst
            psi_est = psi_est + psi/n_pst
            sigma_est = sigma_est + sigma/n_pst
            phi_est = phi_est + phi/n_pst

            if write_pst == 'TRUE':
                beta_pst[:,[pp]] = beta
                pp += 1

    # convert standardized beta to per-allele beta
    if beta_std == 'FALSE':
        beta_est /= np.sqrt(2.0*maf*(1.0-maf))

        if write_pst == 'TRUE':
            beta_pst /= np.sqrt(2.0*maf*(1.0-maf))


    # write posterior effect sizes
    if phi_updt == True:
        eff_file = out_dir + '_pst_eff_a%d_b%.1f_phiauto_chr%d.txt' % (a, b, chrom)
    else:
        eff_file = out_dir + '_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt' % (a, b, phi, chrom)

    with open(eff_file, 'w') as ff:
        if write_pst == 'TRUE':
            for snp, bp, a1, a2, beta in zip(sst_dict['SNP'], sst_dict['BP'], sst_dict['A1'], sst_dict['A2'], beta_pst):
                ff.write(('%d\t%s\t%d\t%s\t%s' + '\t%.6e'*n_pst + '\n') % (chrom, snp, bp, a1, a2, *beta))
        else:
            for snp, bp, a1, a2, beta in zip(sst_dict['SNP'], sst_dict['BP'], sst_dict['A1'], sst_dict['A2'], beta_est):
                ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    # write posterior estimates of psi
    if write_psi == 'TRUE':
        if phi_updt == True:
            psi_file = out_dir + '_pst_psi_a%d_b%.1f_phiauto_chr%d.txt' % (a, b, chrom)
        else:
            psi_file = out_dir + '_pst_psi_a%d_b%.1f_phi%1.0e_chr%d.txt' % (a, b, phi, chrom)

        with open(psi_file, 'w') as ff:
            for snp, psi in zip(sst_dict['SNP'], psi_est):
                ff.write('%s\t%.6e\n' % (snp, psi))

    # print estimated phi
    if phi_updt == True:
        print('... Estimated global shrinkage parameter: %1.2e ...' % phi_est )

    print('... Done ...')

def parse_ref(ref_file, chrom):
    print('... parse reference file: %s ...' % ref_file)

    ref_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'MAF':[]}
    with open(ref_file) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom:
                ref_dict['CHR'].append(chrom)
                ref_dict['SNP'].append(ll[1])
                ref_dict['BP'].append(int(ll[2]))
                ref_dict['A1'].append(ll[3])
                ref_dict['A2'].append(ll[4])
                ref_dict['MAF'].append(float(ll[5]))

    print('... %d SNPs on chromosome %d read from %s ...' % (len(ref_dict['SNP']), chrom, ref_file))
    return ref_dict

def parse_bim(bim_file, chrom):
    print('... parse bim file: %s ...' % (bim_file + '.bim'))

    vld_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(bim_file + '.bim') as ff:
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom:
                vld_dict['SNP'].append(ll[1])
                vld_dict['A1'].append(ll[4])
                vld_dict['A2'].append(ll[5])

    print('... %d SNPs on chromosome %d read from %s ...' % (len(vld_dict['SNP']), chrom, bim_file + '.bim'))
    return vld_dict

def parse_sumstats(ref_dict, vld_dict, sst_file, n_subj):
    print('... parse sumstats file: %s ...' % sst_file)

    ATGC = ['A', 'T', 'G', 'C']
    sst_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(sst_file) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            if ll[1] in ATGC and ll[2] in ATGC:
                sst_dict['SNP'].append(ll[0])
                sst_dict['A1'].append(ll[1])
                sst_dict['A2'].append(ll[2])

    print('... %d SNPs read from %s ...' % (len(sst_dict['SNP']), sst_file))


    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    vld_snp = set(zip(vld_dict['SNP'], vld_dict['A1'], vld_dict['A2']))

    ref_snp = set(zip(ref_dict['SNP'], ref_dict['A1'], ref_dict['A2'])) | set(zip(ref_dict['SNP'], ref_dict['A2'], ref_dict['A1'])) | \
              set(zip(ref_dict['SNP'], [mapping[aa] for aa in ref_dict['A1']], [mapping[aa] for aa in ref_dict['A2']])) | \
              set(zip(ref_dict['SNP'], [mapping[aa] for aa in ref_dict['A2']], [mapping[aa] for aa in ref_dict['A1']]))
    
    sst_snp = set(zip(sst_dict['SNP'], sst_dict['A1'], sst_dict['A2'])) | set(zip(sst_dict['SNP'], sst_dict['A2'], sst_dict['A1'])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A1']], [mapping[aa] for aa in sst_dict['A2']])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A2']], [mapping[aa] for aa in sst_dict['A1']]))

    comm_snp = vld_snp & ref_snp & sst_snp

    print('... %d common SNPs in the reference, sumstats, and validation set ...' % len(comm_snp))


    n_sqrt = np.sqrt(n_subj)
    sst_eff = {}
    with open(sst_file) as ff:
        header = (next(ff).strip()).split()
        header = [col.upper() for col in header]
        for line in ff:
            ll = (line.strip()).split()
            snp = ll[0]; a1 = ll[1]; a2 = ll[2]
            if a1 not in ATGC or a2 not in ATGC:
                continue
            if (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                if 'BETA' in header:
                    beta = float(ll[3])
                elif 'OR' in header:
                    beta = np.log(float(ll[3]))

                if 'SE' in header:
                    se = float(ll[4])
                    beta_std = beta/se/n_sqrt
                elif 'P' in header:
                    p = max(float(ll[4]), 1e-323)
                    beta_std = np.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt

                sst_eff.update({snp: beta_std})

            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                if 'BETA' in header:
                    beta = float(ll[3])
                elif 'OR' in header:
                    beta = np.log(float(ll[3]))

                if 'SE' in header:
                    se = float(ll[4])
                    beta_std = -1*beta/se/n_sqrt
                elif 'P' in header:
                    p = max(float(ll[4]), 1e-323)
                    beta_std = -1*np.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt

                sst_eff.update({snp: beta_std})


    sst_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'MAF':[], 'BETA':[], 'FLP':[]}
    for (ii, snp) in enumerate(ref_dict['SNP']):
        if snp in sst_eff:
            sst_dict['SNP'].append(snp)
            sst_dict['CHR'].append(ref_dict['CHR'][ii])
            sst_dict['BP'].append(ref_dict['BP'][ii])
            sst_dict['BETA'].append(sst_eff[snp])

            a1 = ref_dict['A1'][ii]; a2 = ref_dict['A2'][ii]
            if (snp, a1, a2) in comm_snp:
                sst_dict['A1'].append(a1)
                sst_dict['A2'].append(a2)
                sst_dict['MAF'].append(ref_dict['MAF'][ii])
                sst_dict['FLP'].append(1)
            elif (snp, a2, a1) in comm_snp:
                sst_dict['A1'].append(a2)
                sst_dict['A2'].append(a1)
                sst_dict['MAF'].append(1-ref_dict['MAF'][ii])
                sst_dict['FLP'].append(-1)
            elif (snp, mapping[a1], mapping[a2]) in comm_snp:
                sst_dict['A1'].append(mapping[a1])
                sst_dict['A2'].append(mapping[a2])
                sst_dict['MAF'].append(ref_dict['MAF'][ii])
                sst_dict['FLP'].append(1)
            elif (snp, mapping[a2], mapping[a1]) in comm_snp:
                sst_dict['A1'].append(mapping[a2])
                sst_dict['A2'].append(mapping[a1])
                sst_dict['MAF'].append(1-ref_dict['MAF'][ii])
                sst_dict['FLP'].append(-1)

    return sst_dict

def parse_ldblk(ldblk_dir, sst_dict, chrom):
    print('... parse reference LD on chromosome %d ...' % chrom)

    if '1kg' in os.path.basename(ldblk_dir):
        chr_name = ldblk_dir + '/ldblk_1kg_chr' + str(chrom) + '.hdf5'
    elif 'ukbb' in os.path.basename(ldblk_dir):
        chr_name = ldblk_dir + '/ldblk_ukbb_chr' + str(chrom) + '.hdf5'

    hdf_chr = h5py.File(chr_name, 'r')
    n_blk = len(hdf_chr)
    ld_blk = [np.array(hdf_chr['blk_'+str(blk)]['ldblk']) for blk in range(1,n_blk+1)]

    snp_blk = []
    for blk in range(1,n_blk+1):
        snp_blk.append([bb.decode("UTF-8") for bb in list(hdf_chr['blk_'+str(blk)]['snplist'])])

    blk_size = []
    mm = 0
    for blk in range(n_blk):
        idx = [ii for (ii, snp) in enumerate(snp_blk[blk]) if snp in sst_dict['SNP']]
        blk_size.append(len(idx))
        if idx != []:
            idx_blk = range(mm,mm+len(idx))
            flip = [sst_dict['FLP'][jj] for jj in idx_blk]
            ld_blk[blk] = ld_blk[blk][np.ix_(idx,idx)]*np.outer(flip,flip)

            _, s, v = linalg.svd(ld_blk[blk])
            h = np.dot(v.T, np.dot(np.diag(s), v))
            ld_blk[blk] = (ld_blk[blk]+h)/2            

            mm += len(idx)
        else:
            ld_blk[blk] = np.array([])

    return ld_blk, blk_size

def parse_param():
    long_opts_list = ['ref_dir=', 'bim_prefix=', 'sst_file=', 'a=', 'b=', 'phi=', 'n_gwas=',
                      'n_iter=', 'n_burnin=', 'thin=', 'out_dir=', 'chrom=', 'beta_std=', 'write_psi=', 'write_pst=', 'seed=', 'help']

    param_dict = {'ref_dir': None, 'bim_prefix': None, 'sst_file': None, 'a': 1, 'b': 0.5, 'phi': None, 'n_gwas': None,
                  'n_iter': 1000, 'n_burnin': 500, 'thin': 5, 'out_dir': None, 'chrom': range(1,23),
                  'beta_std': 'FALSE', 'write_psi': 'FALSE', 'write_pst': 'FALSE', 'seed': None}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)          
        except:
            print('Option not recognized.')
            print('Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--ref_dir": param_dict['ref_dir'] = arg
            elif opt == "--bim_prefix": param_dict['bim_prefix'] = arg
            elif opt == "--sst_file": param_dict['sst_file'] = arg
            elif opt == "--a": param_dict['a'] = float(arg)
            elif opt == "--b": param_dict['b'] = float(arg)
            elif opt == "--phi": param_dict['phi'] = float(arg)
            elif opt == "--n_gwas": param_dict['n_gwas'] = int(arg)
            elif opt == "--n_iter": param_dict['n_iter'] = int(arg)
            elif opt == "--n_burnin": param_dict['n_burnin'] = int(arg)
            elif opt == "--thin": param_dict['thin'] = int(arg)
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--chrom": param_dict['chrom'] = arg.split(',')
            elif opt == "--beta_std": param_dict['beta_std'] = arg.upper()
            elif opt == "--write_psi": param_dict['write_psi'] = arg.upper()
            elif opt == "--write_pst": param_dict['write_pst'] = arg.upper()
            elif opt == "--seed": param_dict['seed'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)

    if param_dict['ref_dir'] == None:
        print('* Please specify the directory to the reference panel using --ref_dir\n')
        sys.exit(2)
    elif param_dict['bim_prefix'] == None:
        print('* Please specify the directory and prefix of the bim file for the target dataset using --bim_prefix\n')
        sys.exit(2)
    elif param_dict['sst_file'] == None:
        print('* Please specify the summary statistics file using --sst_file\n')
        sys.exit(2)
    elif param_dict['n_gwas'] == None:
        print('* Please specify the sample size of the GWAS using --n_gwas\n')
        sys.exit(2)
    elif param_dict['out_dir'] == None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict

def main():
    param_dict = parse_param()

    for chrom in param_dict['chrom']:
        print('##### process chromosome %d #####' % int(chrom))

        if '1kg' in os.path.basename(param_dict['ref_dir']):
            ref_dict = parse_ref(param_dict['ref_dir'] + '/snpinfo_1kg_hm3', int(chrom))
        elif 'ukbb' in os.path.basename(param_dict['ref_dir']):
            ref_dict = parse_ref(param_dict['ref_dir'] + '/snpinfo_ukbb_hm3', int(chrom))

        vld_dict = parse_bim(param_dict['bim_prefix'], int(chrom))

        sst_dict = parse_sumstats(ref_dict, vld_dict, param_dict['sst_file'], param_dict['n_gwas'])

        ld_blk, blk_size = parse_ldblk(param_dict['ref_dir'], sst_dict, int(chrom))

        mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], sst_dict, param_dict['n_gwas'], ld_blk, blk_size,
            param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], int(chrom), param_dict['out_dir'], param_dict['beta_std'],
	    param_dict['write_psi'], param_dict['write_pst'], param_dict['seed'])

        print('\n')

if __name__ == '__main__':
    main()