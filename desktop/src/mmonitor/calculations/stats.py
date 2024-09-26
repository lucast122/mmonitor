from typing import Iterable, Tuple, Dict

from numpy import arctanh
from pandas import Series
from scipy.stats import pearsonr, spearmanr, kendalltau, bootstrap

# available SciPy stats functions
_SCIPY_FUNS = {
    'Pearson': pearsonr,
    'Spearman': spearmanr,
    'Kendall': kendalltau
}


def scipy_correlation(xs: Series, ys: Series, tests: Iterable[str], scipy_fun: str) -> Tuple[float, Dict[str, float]]:
    """
    Calculate a statistical correlation between two series and its probability available in SciPy.
    """

    # calculate correlation and T-Test probability
    fun = _SCIPY_FUNS[scipy_fun]
    score, t_test_score = fun(xs, ys)

    # calculate other selected probabilities
    test_scores = dict()
    # T-Test score is already calculated by calling the SciPy function
    if 'T-Test' in tests:
        test_scores['T-Test'] = t_test_score

    # Bootstrapping is very slow!
    if 'Bootstrap' in tests:
        def statistic(l1, l2):
            return fun(l1, l2)[0]

        res = bootstrap((xs, ys), statistic, vectorized=False, paired=True)
        bs_score = (res.confidence_interval.low + res.confidence_interval.high) / 2
        test_scores['Bootstrap'] = bs_score

    # Fisher Z Transform is equivalent to arctanh
    if 'Fisher Z Transform' in tests:
        test_scores['Fisher Z Transform'] = arctanh(score)

    return score, test_scores
