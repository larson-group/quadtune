import numpy as np
from scipy.stats import norm


def bootstrapCalculations(numSamples,
                          metricsWeights,
                          metricsNames,
                          paramsNames,
                          numMetrics,
                          numMetricsToTune,
                          normMetricValsCol,
                          magParamValsRow,
                          defaultParamValsOrigRow,
                          normlzdSensMatrixPoly,
                          normlzdDefaultBiasesCol,
                          normlzdCurvMatrix,
                          reglrCoef,
                          defaultBiasesCol):
    """
    Performs bootstrap resampling on the metric set and solves for parameter values for each bootstrap sample.
    Also computes confidence intervals, residuals, and trade-off diagnostics.

    Parameters:
        numSamples (int): Number of bootstrap samples to generate.
        metricsWeights (np.ndarray): Metric weights, with shape (n_metrics, 1).
        metricsNames (np.ndarray): Names of the metrics, with shape (n_metrics,).
        paramsNames (np.ndarray): Names of the parameters, with shape (n_params,).
        numMetrics (int): Total number of metrics.
        numMetricsToTune (int): Number of metrics used for tuning (sampled during bootstrapping).
        normMetricValsCol (np.ndarray): Normalized observed metric values, with shape (n_metrics, 1).
        magParamValsRow (np.ndarray): Magnitudes of parameter values, with shape (1, n_params).
        defaultParamValsOrigRow (np.ndarray): Default/original parameter values, with shape (1, n_params).
        normlzdSensMatrixPoly (np.ndarray): Normalized sensitivity matrix, with shape (n_metrics, n_params).
        normlzdDefaultBiasesCol (np.ndarray): Normalized model biases for default parameters, with shape (n_metrics, 1).
        normlzdCurvMatrix (np.ndarray): Normalized curvature matrix, with shape (n_metrics, n_params).
        reglrCoef (float): Regularization coefficient.
        defaultBiasesCol (np.ndarray): Raw model biases for default parameters, with shape (n_metrics, 1).

    Returns:
        paramsBoot (np.ndarray): Bootstrap samples of parameters, with shape (n_samples, n_params).
        paramsTuned (np.ndarray): Tuned parameters from the full dataset, with shape (n_params,).
        residualsDefaultCol (np.ndarray): Residuals for default parameters, with shape (n_metrics, 1).
        residualsTunedCol (np.ndarray): Residuals for tuned parameters, with shape (n_metrics, 1).
        residualsBootstrapMatrix (np.ndarray): Residuals from bootstrap samples, with shape (n_samples, n_metrics).
        paramBoundsBoot (np.ndarray): Percentile and BCa confidence interval bounds,
                                      with shape (2, 3, n_params) where axis 0 indexes [Percentile, BCa],
                                      axis 1 indexes [lower, tuned, upper], and axis 2 is over parameters.
        normResidualPairsMatrix (np.ndarray): Matrix of pairwise minimum norms of residual vectors,
                                              with shape (n_metrics, n_metrics).
        tradeoffBinaryMatrix (np.ndarray): Binary matrix indicating metric trade-offs,
                                           with shape (n_metrics, n_metrics).
    """
    from quadtune_driver import solveUsingNonlin, fwdFnc
    # In order to do sampling with replacement of the metrics we sample the corresponding indices with replacement
    metricsSampleIdxMatrix = np.random.randint(0, numMetricsToTune, size=(numSamples, numMetricsToTune))

    paramsBoot = np.zeros((numSamples, len(paramsNames), 1))
    defaultBiasesApproxNonlinMatrix = np.full((numSamples, numMetrics, 1), np.nan)

    # Loop through every bootstrap sample and calculate parameter values
    for sampleIdx, metricsSampleIdxRow in enumerate(metricsSampleIdxMatrix):
        if sampleIdx % 10 == 0:
            print(f"Progress: {sampleIdx}/{numSamples}")
        defaultBiasesApproxNonlin, dnormlzdParamsSolnNonlin, paramsBoot[sampleIdx], \
            dnormlzdParamsSolnLin, paramsSolnLin, defaultBiasesApproxNonlin2x, \
            defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv = (
            solveUsingNonlin(metricsNames[metricsSampleIdxRow],
                             metricsWeights[metricsSampleIdxRow],
                             normMetricValsCol[metricsSampleIdxRow],
                             magParamValsRow,
                             defaultParamValsOrigRow,
                             normlzdSensMatrixPoly[metricsSampleIdxRow, :],
                             normlzdDefaultBiasesCol[metricsSampleIdxRow],
                             normlzdCurvMatrix[metricsSampleIdxRow, :],
                             reglrCoef,
                             beVerbose=False))

        defaultBiasesApproxNonlinMatrix[sampleIdx, :] = fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrixPoly,
                                                               normlzdCurvMatrix, numMetrics) * np.abs(
            normMetricValsCol)
    print(f"Progress: {numSamples}/{numSamples}")
    paramsBoot = paramsBoot[:, :, 0]

    # get solution of the full dataset
    biasesTuned, _, paramsTuned, *_ = solveUsingNonlin(metricsNames,
                                                       metricsWeights,
                                                       normMetricValsCol,
                                                       magParamValsRow,
                                                       defaultParamValsOrigRow,
                                                       normlzdSensMatrixPoly,
                                                       normlzdDefaultBiasesCol,
                                                       normlzdCurvMatrix,
                                                       reglrCoef,
                                                       beVerbose=False)
    paramsTuned = paramsTuned[:, 0]
    # compute residuals
    residualsDefaultCol = - defaultBiasesCol
    residualsTunedCol = -biasesTuned - defaultBiasesCol
    residualsBootstrapMatrix = -defaultBiasesApproxNonlinMatrix[:, :, 0] - defaultBiasesCol.T

    # lower and upper bounds for error pars plot
    ciLowerPercentile, ciUpperPercentile = percentileIntervals(paramsBoot)
    jackknife_params = computeJackknifeParams(metricsNames, paramsNames, metricsWeights, normMetricValsCol,
                                              magParamValsRow, defaultParamValsOrigRow, normlzdSensMatrixPoly,
                                              normlzdDefaultBiasesCol, normlzdCurvMatrix, reglrCoef)
    ciLowerBca, ciUpperBca = bcaIntervals(paramsBoot, paramsTuned, jackknife_params)
    paramBoundsBoot = np.array(
        [[ciLowerPercentile, paramsTuned, ciUpperPercentile], [ciLowerBca, paramsTuned, ciUpperBca]])

    normResidualPairsMatrix = computeNormResidualPairsMatrix(residualsBootstrapMatrix)
    tradeoffBinaryMatrix = computeTradeoffBinaryMatrix(residualsBootstrapMatrix, normResidualPairsMatrix)

    return paramsBoot, paramsTuned, residualsDefaultCol, residualsTunedCol, residualsBootstrapMatrix, paramBoundsBoot, normResidualPairsMatrix, tradeoffBinaryMatrix


def computeJackknifeParams(metricsNames, paramsNames, metricsWeights, normMetricValsCol, magParamValsRow,
                           defaultParamValsOrigRow, normlzdSensMatrixPoly, normlzdDefaultBiasesCol, normlzdCurvMatrix,
                           reglrCoef):
    """
    Computes jackknife estimates of tuned parameters by excluding each metric once and resolving the tuning.

    Parameters:
        metricsNames (np.ndarray): Names of metrics, with shape (n_metrics,).
        paramsNames (np.ndarray): Names of parameters, with shape (n_params,).
        metricsWeights (np.ndarray): Weights for each metric, with shape (n_metrics, 1).
        normMetricValsCol (np.ndarray): Normalized observed metric values, with shape (n_metrics, 1).
        magParamValsRow (np.ndarray): Magnitudes of parameter values, with shape (1, n_params).
        defaultParamValsOrigRow (np.ndarray): Default parameter values, with shape (1, n_params).
        normlzdSensMatrixPoly (np.ndarray): Normalized sensitivity matrix, with shape (n_metrics, n_params).
        normlzdDefaultBiasesCol (np.ndarray): Normalized default model biases, with shape (n_metrics, 1).
        normlzdCurvMatrix (np.ndarray): Normalized curvature matrix, with shape (n_metrics, n_params).
        reglrCoef (float): Regularization coefficient.

    Returns:
        jacknife_params (np.ndarray): Jackknife parameter estimates, with shape (n_metrics, n_params).
    """
    from quadtune_driver import solveUsingNonlin
    print("Computing jackknife estimates . . .")
    jacknifeParams = np.zeros((len(metricsNames), len(paramsNames), 1))
    for i in range(len(metricsNames)):
        _, _, jacknifeParams[i], *_ = solveUsingNonlin(np.delete(metricsNames, i),
                                                        np.delete(metricsWeights, i, axis=0),
                                                        np.delete(normMetricValsCol, i, axis=0),
                                                        magParamValsRow,
                                                        defaultParamValsOrigRow,
                                                        np.delete(normlzdSensMatrixPoly, i, axis=0),
                                                        np.delete(normlzdDefaultBiasesCol, i, axis=0),
                                                        np.delete(normlzdCurvMatrix, i, axis=0),
                                                        reglrCoef,
                                                        beVerbose=False)
    return jacknifeParams[:, :, 0]


def bcaIntervals(paramsBoot, paramsTuned, jackknifeParams, alpha=0.05):
    """
    Computes BCa (bias-corrected and accelerated) bootstrap confidence intervals for parameter estimates.

    Parameters:
        paramsBoot (np.ndarray): Bootstrap samples of parameters, with shape (n_samples, n_params).
        paramsTuned (np.ndarray): Tuned parameter values used as reference points, with shape (n_params,).
        jackknifeParams (np.ndarray): Jackknife estimates of parameters, with shape (n_metrics, n_params).
        alpha (float): Significance level for the confidence interval (default is 0.05 for a 95% CI).

    Returns:
        ciLower (np.ndarray): Lower bounds of the BCa confidence intervals, with shape (n_params,).
        ciUpper (np.ndarray): Upper bounds of the BCa confidence intervals, with shape (n_params,).
    """
    B, p = paramsBoot.shape

    ciLower = np.zeros(p)
    ciUpper = np.zeros(p)

    for j in range(p):
        thetaBoot = paramsBoot[:, j]
        thetaHat = paramsTuned[j]
        thetaJack = jackknifeParams[:, j]

        # 1. Bias correction z0
        proportion = np.mean(thetaBoot < thetaHat)
        z0 = norm.ppf(proportion)

        # 2. Acceleration a (jackknife skewness)
        theta_bar = np.mean(thetaJack)
        num = np.sum((theta_bar - thetaJack) ** 3)
        denom = 6 * (np.sum((theta_bar - thetaJack) ** 2) ** 1.5)
        a = num / denom if denom != 0 else 0.0  # avoid division by zero

        # 3. Adjusted alpha levels
        zAlphaLow = norm.ppf(alpha / 2)
        zAlphaHigh = norm.ppf(1 - alpha / 2)

        pctLow = norm.cdf(z0 + (z0 + zAlphaLow) / (1 - a * (z0 + zAlphaLow)))
        pctHigh = norm.cdf(z0 + (z0 + zAlphaHigh) / (1 - a * (z0 + zAlphaHigh)))

        ciLower[j] = np.percentile(thetaBoot, 100 * pctLow)
        ciUpper[j] = np.percentile(thetaBoot, 100 * pctHigh)

    return ciLower, ciUpper


def percentileIntervals(paramsBoot, alpha=0.05):
    """
    Computes percentile confidence intervals for bootstrap parameter samples.

    Parameters:
        paramsBoot (np.ndarray): Bootstrap samples of parameters, with shape (n_samples, n_params).
        alpha (float): Significance level for the confidence interval (default is 0.05 for a 95% CI).

    Returns:
        ciLower (np.ndarray): Lower bounds of the confidence intervals, with shape (n_params,).
        ciUpper (np.ndarray): Upper bounds of the confidence intervals, with shape (n_params,).
    """
    ciLower = np.percentile(paramsBoot, 100 * (alpha / 2), axis=0)
    ciUpper = np.percentile(paramsBoot, 100 * (1 - alpha / 2), axis=0)
    return ciLower, ciUpper


def computeNormResidualPairsMatrix(residualsBootstrapMatrix):
    """
    Computes a lower triangular matrix of minimum pairwise norms between residual vectors.

    Parameters:
        residualsBootstrapMatrix (np.ndarray): Bootstrap residuals, with shape (n_samples, n_metrics).

    Returns:
        result (np.ndarray): Lower triangular matrix of pairwise residual norms,
                             with shape (n_metrics, n_metrics), filled with NaNs above the diagonal.
    """
    n_cols = residualsBootstrapMatrix.shape[1]
    result = np.full((n_cols, n_cols), np.nan)

    for i in range(n_cols):
        for j in range(0, i + 1):
            vecs = np.stack((residualsBootstrapMatrix[:, i], residualsBootstrapMatrix[:, j]),
                            axis=1)  # shape: (n_rows, 2)
            norms = np.linalg.norm(vecs, axis=1)
            result[i, j] = np.min(norms)
    return result


def computeTradeoffBinaryMatrix(residualsBootstrapMatrix, normResidualPairsMatrix):
    """
    Computes a lower triangular binary matrix indicating metric pairs that exhibit potential trade-offs.

    Parameters:
        residualsBootstrapMatrix (np.ndarray): Bootstrap residuals, with shape (n_samples, n_metrics).
        normResidualPairsMatrix (np.ndarray): Pairwise min norm matrix of residuals, with shape (n_metrics, n_metrics).

    Returns:
        tradeoffBinaryMatrix (np.ndarray): Binary matrix indicating trade-off presence between metric pairs,
                                           with shape (n_metrics, n_metrics).
    """
    n = normResidualPairsMatrix.shape[0]
    threshold = 1 / np.sqrt(2)
    correlation_threshold = 0.8
    valid_pairs = []

    # Binary matrix to mark valid pairs
    tradeoffBinaryMatrix = np.zeros((n, n), dtype=int)

    # Check all lower diagonal pairs
    for i in range(n):
        for j in range(i):
            d_i = normResidualPairsMatrix[i, i]
            d_j = normResidualPairsMatrix[j, j]
            off_diag = normResidualPairsMatrix[i, j]

            if abs(d_i) < threshold and abs(d_j) < threshold and abs(off_diag) > threshold:
                corr = np.corrcoef(residualsBootstrapMatrix[:, i], residualsBootstrapMatrix[:, j])[0, 1]
                if abs(corr) > correlation_threshold:
                    tradeoffBinaryMatrix[i, j] = 1
                    valid_pairs.append((i, j))

    return tradeoffBinaryMatrix
