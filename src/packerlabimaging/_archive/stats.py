# class for collecting statistical tests for analysing cellsdata
import numpy as np
import scipy.stats as stats
import statsmodels.api
import statsmodels as sm

from packerlabimaging import AllOpticalTrial

"""STATISTICS FOR TESTING FOR SIGNIFICANCE OF PHOTOSTIM RESPONSES"""


def runWilcoxonsTest(array1, array2):
    """
    Run wilcoxons statistical comparison test on two input arrays. Array 1 and array 2 must be of same length and ordered in paired format.
    """
    # check if the two distributions of flu values (pre/post) are different
    assert array1.shape == array2.shape, 'shapes for .__pre_array and .__post_array need to be the same for wilcoxon test'
    wilcoxons = np.empty(len(array1))  # [cell (p-value)]

    for cell in range(len(array1)):
        wilcoxons[cell] = stats.wilcoxon(array2[cell], array1[cell])[1]

    return wilcoxons


def sigTestAvgResponse(trial: AllOpticalTrial, p_vals: list, alpha=0.1):
    """
    Uses the p values and a threshold for the Benjamini-Hochberg correction to return which
    cells are still significant after correcting for multiple significance testing

    :param trial:
    :param p_vals: list of p-values, each corresponding to a given cell's comparison
    :param alpha:
    :return:
    """
    print('\n----------------------------------------------------------------')
    print('running statistical significance testing for nontargets response arrays ')
    print('----------------------------------------------------------------')

    sig_units = np.full_like(p_vals, False, dtype=bool)

    try:
        sig_units, _, _, _ = sm.stats.multitest.multipletests(p_vals, alpha=alpha, method='fdr_bh',
                                                              is_sorted=False, returnsorted=False)
    except ZeroDivisionError:
        print('no cells responding')

    # p values without bonferroni correction
    no_bonf_corr = [i for i, p in enumerate(p_vals) if p < 0.05]
    nomulti_sig_units = np.zeros(len(p_vals), dtype='bool')
    nomulti_sig_units[no_bonf_corr] = True

    # p values after bonferroni correction
    bonf_corr = [i for i, p in enumerate(p_vals) if p < 0.05 / len(p_vals)]
    sig_units = np.zeros(trial.Suite2p.n_units, dtype='bool')
    sig_units[bonf_corr] = True

    return sig_units
