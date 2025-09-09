# -*- coding: utf-8 -*-

# Run this app with `python3 quadtune_driver.py` and
# view the plots at http://127.0.0.1:8050/ in your web browser.
# (To open a web browser on a larson-group computer,
# login to malan with `ssh -X` and then type `firefox &`.)

import numpy as np
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.interpolate import UnivariateSpline
from sklearn import linear_model

import matplotlib.pyplot as plt

def main():
    """
    Main driver for QuadTune.
    It calls routines to feeds in input configuration data,
    find optimal parameter values, and
    create diagnostic plots.
    """

    #from config_default import setUpConfig
    #from config_example import setUpConfig
    from config_example_eam import setUpConfig
    from set_up_inputs \
        import setUpColAndRowVectors, \
               setUpDefaultMetricValsCol, \
               createInteractIdxs, \
               calcNormlzdInteractBiasesCols, \
               readDnormlzdParamsInteract, \
               calcInteractDerivs, \
               checkInteractParamVals, \
               printInteractDiagnostics, \
               checkInteractDerivs

    from create_nonbootstrap_figs import createFigs
    from create_bootstrap_figs import bootstrapPlots
    from do_bootstrap_calcs import bootstrapCalculations


    print("Set up inputs . . .")

    # The user should input all tuning configuration info into file set_up_inputs.py
    (numMetricsNoCustom, numMetricsToTune,
     metricsNames, metricsNamesNoprefix,
     varPrefixes, mapVarIdx, boxSize,
     highlightedMetricsToPlot, createPlotType,
     metricsWeights, metricsNorms,
     obsMetricValsDict,
     obsOffsetCol, obsGlobalAvgCol, doObsOffset,
     paramsNames, paramsScales,
     transformedParamsNames,
     prescribedParamsNames, prescribedParamsScales,
     prescribedTransformedParamsNames,
     prescribedParamValsRow,
     prescribedSensNcFilenames, prescribedSensNcFilenamesExt,
     sensNcFilenames, sensNcFilenamesExt,
     sensSST4KNcFilenames, sensSST4KNcFilenamesExt,
     defaultNcFilename, globTunedNcFilename,
     defaultSST4KNcFilename,
     interactParamsNamesAndFilenames,
     doPiecewise,
     reglrCoef, doBootstrapSampling, numBootstrapSamples) \
    = \
        setUpConfig(beVerbose=False)

    # Number of regional metrics, including all of varPrefixes including the metrics we're not tuning, plus custom regions.
    numMetrics = len(metricsNames)

    # We apply a tiny weight to the final metrics.
    #    Those metrics will appear in the diagnostics
    #    but their errors will not be accounted for in tuning.
    metricsWeights[numMetricsToTune:] = 1e-12

    print("Set up preliminaries . . .")

    obsMetricValsCol, normMetricValsCol, \
    defaultBiasesCol, \
    defaultParamValsOrigRow, \
    sensParamValsRow, sensParamValsRowExt, \
    dnormlzdSensParams, \
    magParamValsRow, \
    dnormlzdPrescribedParams, \
    magPrescribedParamValsRow, \
    = setUpColAndRowVectors(metricsNames, metricsNorms,
                            obsMetricValsDict,
                            obsOffsetCol, obsGlobalAvgCol, doObsOffset,
                            paramsNames, transformedParamsNames, prescribedParamsNames, prescribedParamValsRow,
                            prescribedTransformedParamsNames,
                            sensNcFilenames, sensNcFilenamesExt,
                            defaultNcFilename
                            )

    # Construct numMetrics x numParams matrix of second derivatives, d2metrics/dparams2.
    #     The derivatives are normalized by observed metric values and max param values.
    # Also construct a linear sensitivity matrix, dmetrics/dparams.
    normlzdCurvMatrix, normlzdSensMatrixPoly, normlzdConstMatrix, \
    normlzdOrdDparamsMin, normlzdOrdDparamsMax, \
    normlzd_pMid, normlzdLeftSensMatrix, normlzdRightSensMatrix = \
        constructNormlzdSensCurvMatrices(metricsNames, paramsNames, transformedParamsNames,
                                   normMetricValsCol, magParamValsRow,
                                   sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)

    # In order to weight certain metrics, multiply each row of normlzdSensMatrixPoly
    # by metricsWeights
    normlzdWeightedSensMatrixPoly = np.diag(np.transpose(metricsWeights)[0]) @ normlzdSensMatrixPoly

    # interactIdxs = array of numInteractTerms (j,k) tuples of parameter indices of interaction terms
    interactIdxs = createInteractIdxs(interactParamsNamesAndFilenames, paramsNames)

    checkInteractParamVals(
        interactIdxs, interactParamsNamesAndFilenames,
        sensParamValsRow, sensParamValsRowExt,
        defaultParamValsOrigRow,
        paramsNames, transformedParamsNames,
        len(paramsNames))

    defaultMetricValsCol = obsMetricValsCol+defaultBiasesCol

    # normlzdInteractBiasesCols = numMetrics x numInteractTerms array
    normlzdInteractBiasesCols = \
              calcNormlzdInteractBiasesCols( defaultMetricValsCol,
                              normMetricValsCol,
                              metricsNames,
                              interactParamsNamesAndFilenames)

    # dnormlzdParamsInteract = array of numInteractTerms tuples of parameter *values*
    dnormlzdParamsInteract = \
        readDnormlzdParamsInteract(interactParamsNamesAndFilenames, interactIdxs,
                               defaultParamValsOrigRow, magParamValsRow,
                               paramsNames, transformedParamsNames, len(paramsNames))

    # normlzdInteractDerivs = numMetrics x numInteractTerms array
    normlzdInteractDerivs = calcInteractDerivs(interactIdxs,
                       dnormlzdParamsInteract,
                       normlzdInteractBiasesCols,
                       normlzdCurvMatrix, normlzdSensMatrixPoly,
                       numMetrics)

    checkInteractDerivs(normlzdInteractBiasesCols,
                        dnormlzdParamsInteract,
                        len(paramsNames),
                        normlzdSensMatrixPoly, normlzdCurvMatrix,
                        numMetrics,
                        normlzdInteractDerivs, interactIdxs)

    # For prescribed parameters, construct numMetrics x numParams matrix of second derivatives, d2metrics/dparams2.
    # The derivatives are normalized by observed metric values and max param values.
    normlzdPrescribedCurvMatrix, normlzdPrescribedSensMatrixPoly, normlzdPrescribedConstMatrix, \
    normlzdPrescribedOrdDparamsMin, normlzdPrescribedOrdDparamsMax, \
    normlzd_pMidPrescribed, normlzdPrescribedLeftSensMatrix, normlzdPrescribedRightSensMatrix = \
        constructNormlzdSensCurvMatrices(metricsNames, prescribedParamsNames, prescribedTransformedParamsNames,
                                   normMetricValsCol, magPrescribedParamValsRow,
                                   prescribedSensNcFilenames, prescribedSensNcFilenamesExt, defaultNcFilename)

    # This is the prescribed correction to the metrics that appears on the left-hand side of the Taylor equation.
    #   It is not a bias from the obs.  It is a correction to the simulated default metric values
    #   based on prescribed param values.

    ### THIS CALL DOESN'T ACCOUNT FOR INTERACTIONS!!!
    normlzdPrescribedBiasesCol = \
         fwdFnc(dnormlzdPrescribedParams, normlzdPrescribedSensMatrixPoly, normlzdPrescribedCurvMatrix, numMetrics,
                normlzdInteractDerivs= np.empty(0), interactIdxs = np.empty(0))

    prescribedBiasesCol = normlzdPrescribedBiasesCol * np.abs(normMetricValsCol)

    # defaultBiasesCol + prescribedBiasesCol = -fwdFnc_tuned_params  (see lossFnc).
    #     This lumps the prescribed-parameter adjustment into defaultBiasesCol.
    #        but it may be clearer to separate them out.
    # defaultBiasesCol = default simulation - observations
    defaultBiasesCol = defaultBiasesCol + prescribedBiasesCol

    normlzdDefaultBiasesCol = defaultBiasesCol / np.abs(normMetricValsCol)

    print("Optimizing parameter values . . . ")

    # sValsRatio = a threshold ratio of largest singular value
    #              to the smallest retained singular value.
    # If sValsRatio is large enough, then all singular vectors will be kept.
    # If sValsRatio is 1, then only the first singular vector will be kept.
    sValsRatio = 800.
    normlzdSensMatrixPolySvd = \
        approxMatrixWithSvd(normlzdSensMatrixPoly, sValsRatio, sValsNumToKeep=None, beVerbose=False)
    normlzdCurvMatrixSvd = \
        approxMatrixWithSvd(normlzdCurvMatrix, sValsRatio, sValsNumToKeep=None, beVerbose=False)

    #######################################################################################################
    #
    # Calculate an ensemble of parameter values by doing bootstrap sampling of the regional metrics.
    #
    #######################################################################################################

    if doBootstrapSampling:

        print("Starting boostrap sampling . . .")

        # SST4K:  call constructNormlzdSensCurvMatrices with SST4K sensFiles.

        # For SST4K runs,
        #     construct numMetrics x numParams matrix of second derivatives, d2metrics/dparams2.
        #     The derivatives are normalized by observed metric values and max param values.
        # Also construct a linear sensitivity matrix, dmetrics/dparams.
        normlzdCurvMatrixSST4K, normlzdSensMatrixPolySST4K, normlzdConstMatrixSST4K, \
        normlzdOrdDparamsMinSST4K, normlzdOrdDparamsMaxSST4K, \
        normlzd_pMidSST4K, normlzdLeftSensMatrixSST4K, normlzdRightSensMatrixSST4K = \
            constructNormlzdSensCurvMatrices(metricsNames, paramsNames, transformedParamsNames,
                                             normMetricValsCol, magParamValsRow,
                                             sensSST4KNcFilenames, sensSST4KNcFilenamesExt, defaultSST4KNcFilename)

        # SST4K: Here feed normlzdSensMatrixPoly and normlzdCurvMatrix from SST4K runs into bootstrapCalculations.

        ( paramsBoot, paramsTuned, residualsDefaultCol, residualsTunedCol,
          residualsBootstrapMatrix,
          defaultBiasesApproxNonlinMatrixSST4K,
          dDefaultBiasesApproxNonlinMatrixSST4K,
          paramBoundsBoot,
          normResidualPairsMatrix, tradeoffBinaryMatrix ) = \
        bootstrapCalculations(numBootstrapSamples,
                              metricsWeights,
                              metricsNames,
                              paramsNames,
                              numMetrics,  # numMetrics is redundant, given that we're feeding in metricsNames
                              numMetricsToTune,
                              normMetricValsCol,
                              magParamValsRow,
                              defaultParamValsOrigRow,
                              normlzdSensMatrixPoly,
                              normlzdSensMatrixPolySST4K,
                              normlzdDefaultBiasesCol,
                              normlzdCurvMatrix,
                              normlzdCurvMatrixSST4K,
                              doPiecewise, normlzd_pMid,
                              normlzdLeftSensMatrix, normlzdRightSensMatrix,
                              reglrCoef,
                              defaultBiasesCol)

        print(f"Sample avg of paramsBoot = {np.mean(paramsBoot, axis=0)}")

        bootstrapPlots(numMetricsToTune,  # Should we feed in numMetrics instead??
                       boxSize,
                       metricsNames,
                       residualsBootstrapMatrix,
                       defaultBiasesApproxNonlinMatrixSST4K,
                       dDefaultBiasesApproxNonlinMatrixSST4K,
                       residualsTunedCol,
                       residualsDefaultCol,
                       paramsNames,
                       paramsBoot,
                       paramsTuned,
                       defaultParamValsOrigRow,
                       paramBoundsBoot,
                       normResidualPairsMatrix,
                       tradeoffBinaryMatrix)
    else:
        paramBoundsBoot = None

    #end if doBootstrapSampling

    ########################################
    #
    # Resume non-bootstrap calculations
    #
    #########################################


    defaultBiasesApproxNonlin, \
    dnormlzdParamsSolnNonlin, paramsSolnNonlin, \
    dnormlzdParamsSolnLin, paramsSolnLin, \
    defaultBiasesApproxNonlin2x, \
    defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv = \
        solveUsingNonlin(metricsNames,
                         metricsWeights, normMetricValsCol, magParamValsRow,
                         defaultParamValsOrigRow,
                         normlzdSensMatrixPolySvd, normlzdDefaultBiasesCol,
                         normlzdCurvMatrixSvd,
                         doPiecewise, normlzd_pMid,
                         normlzdLeftSensMatrix, normlzdRightSensMatrix,
                         normlzdInteractDerivs, interactIdxs,
                         reglrCoef,
                         beVerbose=False)

    y_hat_i = defaultBiasesApproxNonlin + defaultBiasesCol + obsMetricValsCol

    #print("Tuned parameter perturbation values (dnormzldParamsSolnNonlin)")
    #for idx in range(0,len(paramsNames)): \
    #    print("{:33s} {:7.7g}".format(paramsNames[idx], dnormlzdParamsSolnNonlin[idx][0] ) )
    print("Tuned parameter values (paramsSolnNonlin)")
    for idx in range(0,len(paramsNames)): \
        print("{:33s} {:7.7g}".format(paramsNames[idx], paramsSolnNonlin[idx][0] ) )

    # Check whether the minimizer actually reduces chisqd
    # Initial value of chisqd, which assumes parameter perturbations are zero
    #normlzdWeightedDefaultBiasesCol = metricsWeights * normlzdDefaultBiasesCol
    chisqdZero = lossFnc(np.zeros_like(defaultParamValsOrigRow),
                         normlzdSensMatrixPoly, normlzdDefaultBiasesCol, metricsWeights,
                         normlzdCurvMatrix, reglrCoef, numMetrics,
                         normlzdInteractDerivs, interactIdxs)  # Should I feed in numMetricsToTune instead??
                                                                    #   But metricsWeights is already set to eps for un-tuned metrics.
    # Optimized value of chisqd, which uses optimal values of parameter perturbations
    chisqdMin = lossFnc(dnormlzdParamsSolnNonlin.T,
                        normlzdSensMatrixPoly, normlzdDefaultBiasesCol, metricsWeights,
                        normlzdCurvMatrix, reglrCoef, numMetrics,
                        normlzdInteractDerivs, interactIdxs)  # Should I feed in numMetricsToTune instead??

    print("chisqdMinRatio (all metrics, non-unity metricsWeights) =", chisqdMin/chisqdZero)

    chisqdUnweightedZero = lossFnc(np.zeros_like(defaultParamValsOrigRow),
                                   normlzdSensMatrixPoly, normlzdDefaultBiasesCol, np.ones_like(metricsWeights),
                                   normlzdCurvMatrix, reglrCoef, numMetrics,
                                   normlzdInteractDerivs, interactIdxs)  # Should I feed in numMetricsToTune instead??
    # Optimized value of chisqd, which uses optimal values of parameter perturbations
    chisqdUnweightedMin = lossFnc(dnormlzdParamsSolnNonlin.T,
                                  normlzdSensMatrixPoly, normlzdDefaultBiasesCol, np.ones_like(metricsWeights),
                                  normlzdCurvMatrix, reglrCoef, numMetrics,
                                  normlzdInteractDerivs, interactIdxs)  # Should I feed in numMetricsToTune instead??

    print("chisqdUnweightedMinRatio (all metrics, metricsWeights=1) =", chisqdUnweightedMin/chisqdUnweightedZero)

    # Set up a column vector of metric values from the global simulation based on optimized
    #     parameter values.
    globTunedMetricValsCol = setUpDefaultMetricValsCol(metricsNames, globTunedNcFilename)

    # Store biases in default simulation, ( global_model - obs )
    globTunedBiasesCol = np.subtract(globTunedMetricValsCol, obsMetricValsCol)
    #globTunedBiasesCol = globTunedBiasesCol + prescribedBiasesCol

    # Check whether the minimizer actually reduces chisqd
    # Initial value of chisqd, which assumes parameter perturbations are zero
    normlzdGlobTunedBiasesCol = globTunedBiasesCol/np.abs(normMetricValsCol)
    chisqdGlobTunedMin = lossFnc(np.zeros_like(defaultParamValsOrigRow),
                                 normlzdSensMatrixPoly, normlzdGlobTunedBiasesCol, metricsWeights,
                                 normlzdCurvMatrix, reglrCoef, numMetrics, # Should I feed in numMetricsToTune instead??
                                 normlzdInteractDerivs, interactIdxs)

    print("chisqdGlobTunedMinRatio =", chisqdGlobTunedMin/chisqdZero)

    chisqdUnweightedGlobTunedMin = lossFnc(np.zeros_like(defaultParamValsOrigRow),
                                           normlzdSensMatrixPoly, normlzdGlobTunedBiasesCol, np.ones_like(metricsWeights),
                                           normlzdCurvMatrix, reglrCoef, numMetrics,  # Should I feed in numMetricsToTune instead??
                                           normlzdInteractDerivs, interactIdxs)

    print("chisqdUnweightedGlobTunedMinRatio =", chisqdUnweightedGlobTunedMin/chisqdUnweightedZero)
    print("-----------------------------------------------------")

    if True:
        printInteractDiagnostics(interactIdxs,
                       normlzdInteractDerivs,
                       dnormlzdParamsSolnNonlin,
                       normlzdCurvMatrix, normlzdSensMatrixPoly,
                       paramsNames, numMetrics)

    ##############################################
    #
    #    Create plots
    #
    ##############################################

    # Find best-fit params by use of the Elastic Net algorithm
    defaultBiasesApproxElastic, defaultBiasesApproxElasticNonlin, \
    dnormlzdParamsSolnElastic, paramsSolnElastic = \
        findParamsUsingElastic(normlzdSensMatrixPoly, normlzdWeightedSensMatrixPoly,
                     defaultBiasesCol, normMetricValsCol, metricsWeights,
                     magParamValsRow, defaultParamValsOrigRow,
                     normlzdCurvMatrix,
                     beVerbose=False)
    #defaultBiasesApproxElasticCheck = ( normlzdWeightedSensMatrixPoly @ dnormlzdParamsSolnElastic ) \
    #                        * np.reciprocal(metricsWeights) * np.abs(normMetricValsCol)
    #print("defaultBiasesApproxElastic = ", defaultBiasesApproxElastic)
    #print("defaultBiasesApproxElasticCheck = ", defaultBiasesApproxElasticCheck)

    normlzdLinplusSensMatrixPoly = normlzdSemiLinMatrixFnc(
                                        dnormlzdParamsSolnNonlin, normlzdSensMatrixPoly,
                                        normlzdCurvMatrix, numMetrics)
    #normlzdWeightedLinplusSensMatrixPoly = np.diag(np.transpose(metricsWeights)[0]) \
    #                                          @ normlzdLinplusSensMatrixPoly

    createFigs(numMetricsNoCustom, metricsNames, metricsNamesNoprefix,
               numMetricsToTune,
               varPrefixes, mapVarIdx, boxSize,
               highlightedMetricsToPlot,
               paramsNames, transformedParamsNames, paramsScales,
               metricsWeights, obsMetricValsCol, normMetricValsCol, magParamValsRow,
               defaultBiasesCol, defaultBiasesApproxNonlin, defaultBiasesApproxElastic,
               defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv,
               normlzdDefaultBiasesCol,
               normlzdCurvMatrix, normlzdSensMatrixPoly, normlzdConstMatrix,
               normlzdInteractDerivs, interactIdxs,
               normlzdOrdDparamsMin, normlzdOrdDparamsMax,
               normlzdWeightedSensMatrixPoly,
               dnormlzdParamsSolnNonlin,
               defaultParamValsOrigRow,
               normlzdGlobTunedBiasesCol, normlzdLinplusSensMatrixPoly,
               paramsSolnLin, dnormlzdParamsSolnLin,
               paramsSolnNonlin,
               paramsSolnElastic, dnormlzdParamsSolnElastic,
               sensNcFilenames, sensNcFilenamesExt, defaultNcFilename,
               createPlotType,
               beVerbose=False, useLongTitle=False, paramBoundsBoot=paramBoundsBoot)

    return

def normlzdSemiLinMatrixFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics):
    """Calculate semi-linear matrix, sensMatrix + curvMatrix*dp, for use in forward solution"""

    normlzdSemiLinMatrix = \
        normlzdSensMatrix \
        + 0.5 * normlzdCurvMatrix * ( np.ones((numMetrics,1)) @ dnormlzdParams.T ).reshape(numMetrics,len(dnormlzdParams))

    return normlzdSemiLinMatrix

def normlzdPiecewiseLinMatrixFnc(dnormlzdParams, normlzd_pMid, normlzdDefaultParamValsRow,
                                 normlzdLeftSensMatrix, normlzdRightSensMatrix):
    """Calculate piecewise-linear matrix for use in forward solution"""

    normlzdPiecewiseLinMatrix = np.zeros_like(normlzdLeftSensMatrix)

    for col in np.arange(len(normlzdDefaultParamValsRow)):
        if ( dnormlzdParams >= ( normlzd_pMid - normlzdDefaultParamValsRow ) ):
            normlzdPiecewiseLinMatrix[:,col] = normlzdRightSensMatrix[:,col]
        else:
            normlzdPiecewiseLinMatrix[:, col] = normlzdLeftSensMatrix[:, col]

    return normlzdPiecewiseLinMatrix

def fwdFncNoInteract(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics):
    """Calculate forward nonlinear solution, normalized but not weighted"""

    normlzdDefaultBiasesApproxNonlin = \
        normlzdSemiLinMatrixFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics) @ dnormlzdParams

    return normlzdDefaultBiasesApproxNonlin

def fwdFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics,
           normlzdInteractDerivs = np.empty(0), interactIdxs = np.empty(0)):
    """Calculate forward nonlinear solution, normalized but not weighted"""

    normlzdDefaultBiasesApproxNonlin = \
        fwdFncNoInteract(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics)

    # dnormlzd_dpj_dpk = ( dp_j * dp_k ) for each interaction term
    dnormlzd_dpj_dpk = calc_dnormlzd_dpj_dpk(dnormlzdParams, interactIdxs)
    interactTerms = normlzdInteractDerivs @ dnormlzd_dpj_dpk

    normlzdDefaultBiasesApproxNonlin += interactTerms

    #normlzdDefaultBiasesApproxNonlin = \
    #    normlzdSemiLinMatrixFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics) @ dnormlzdParams

            #normlzdSensMatrix @ dnormlzdParams \
            #+ 0.5 * normlzdCurvMatrix @ (dnormlzdParams * dnormlzdParams) 

    return normlzdDefaultBiasesApproxNonlin

def calc_dnormlzd_dpj_dpk(dnormlzdParams, interactIdxs):

    '''
    Input:
    dnormlzdParams = A numParams numpy row array of dp
    Output:
    dnormlzd_dpj_dpk = A numInteract numpy col array = [ dp_k * dp_j, ... , ]
    '''
    dnormlzd_dpj_dpk = np.zeros((len(interactIdxs),1))
    for idx, jkTuple in np.ndenumerate(interactIdxs):
        dnormlzd_dpj_dpk[idx,0] = dnormlzdParams[jkTuple[0]][0] * dnormlzdParams[jkTuple[1]][0]

    return dnormlzd_dpj_dpk

def lossFncMetrics(dnormlzdParams, normlzdSensMatrix,
                   normlzdDefaultBiasesCol, metricsWeights,
                   normlzdCurvMatrix, numMetrics,
                   normlzdInteractDerivs = np.empty(0), interactIdxs = np.empty(0)):
    """Each regional component of loss function (including squares)"""

    weightedBiasDiffSqdCol = \
        np.square( (-normlzdDefaultBiasesCol
                    - fwdFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics,
                             normlzdInteractDerivs, interactIdxs)
         ) * metricsWeights )

    return weightedBiasDiffSqdCol

def lossFnc(dnormlzdParams, normlzdSensMatrix, normlzdDefaultBiasesCol, metricsWeights,
            normlzdCurvMatrix, reglrCoef, numMetrics,
            normlzdInteractDerivs = np.empty(0), interactIdxs = np.empty(0)):
    """Define objective function (a.k.a. loss function) that is to be minimized."""

    dnormlzdParams = np.atleast_2d(dnormlzdParams).T # convert from 1d row array to 2d column array
    weightedBiasDiffSqdCol = \
        lossFncMetrics(dnormlzdParams, normlzdSensMatrix,
                       normlzdDefaultBiasesCol, metricsWeights,
                       normlzdCurvMatrix, numMetrics,
                       normlzdInteractDerivs, interactIdxs)
    #weightedBiasDiffSqdCol = \
    #    np.square( (-normlzdDefaultBiasesCol \
    #     - fwdFnc(dnormlzdParams, normlzdSensMatrix, normlzdCurvMatrix, numMetrics) \
    #     ) * metricsWeights )
    # This is the chisqd fnc listed in Eqn. (15.2.2) of Numerical Recipes, 1992.
    # It is like MSE (not RMSE), except that it sums the squares rather than averaging them.
    chisqd = np.sum(weightedBiasDiffSqdCol) \
             + reglrCoef * np.linalg.norm(dnormlzdParams, ord=1)
    #chisqd = np.sqrt(np.sum(weightedBiasDiffSqdCol)) \
    #         + reglrCoef * np.linalg.norm(dnormlzdParams, ord=1)
    #chisqd = np.linalg.norm( weightedBiasDiffCol, ord=2 )**1  \
    #            + reglrCoef * np.linalg.norm( dnormlzdParams, ord=1 )

    return chisqd

def solveUsingNonlin(metricsNames,
                     metricsWeights, normMetricValsCol, magParamValsRow,
                     defaultParamValsOrigRow,
                     normlzdSensMatrix, normlzdDefaultBiasesCol,
                     normlzdCurvMatrix,
                     doPiecewise, normlzd_pMid,
                     normlzdLeftSensMatrix, normlzdRightSensMatrix,
                     normlzdInteractDerivs = np.empty(0), interactIdxs = np.empty(0),
                     reglrCoef = 0.0,
                     beVerbose = False):
    """Find optimal parameter values by minimizing quartic loss function"""

    numMetrics = len(metricsNames)


    # Don't let parameter values go negative
    lowerBoundsCol =  -defaultParamValsOrigRow[0]/magParamValsRow[0]

    #x0TwoYr = np.array([-0.1400083, -0.404022, 0.2203307, -0.9838958, 0.391993, -0.05910007, 1.198831])
    #x0TwoYr = np.array([0.5805136, -0.1447917, -0.2722521, -0.8183079, 0.3150205, -0.4794127, 0.1104284])
    x0TwoYr = np.array([0.5805136, -0.2722521, -0.8183079, 0.3150205, -0.4794127])
    # Perform nonlinear optimization
    #normlzdDefaultBiasesCol = defaultBiasesCol/np.abs(normMetricValsCol)
    #dnormlzdParamsSolnNonlin = minimize(lossFnc,x0=np.ones_like(np.transpose(defaultParamValsOrigRow)), \
    dnormlzdParamsSolnNonlin = (minimize(lossFnc, x0=np.zeros_like(np.transpose(defaultParamValsOrigRow[0])),
                                         #dnormlzdParamsSolnNonlin = minimize(lossFnc,x0=x0TwoYr, \
                                         #dnormlzdParamsSolnNonlin = minimize(lossFnc,dnormlzdParamsSoln, \
                                         args=(normlzdSensMatrix, normlzdDefaultBiasesCol, metricsWeights,
                                               normlzdCurvMatrix, reglrCoef, numMetrics, normlzdInteractDerivs, interactIdxs),
                                         method='Powell', tol=1e-12
                                         ))
                                #,)
                                #        bounds=Bounds(lb=lowerBoundsCol))
    dnormlzdParamsSolnNonlin = np.atleast_2d(dnormlzdParamsSolnNonlin.x).T


    dparamsSolnNonlin = dnormlzdParamsSolnNonlin * np.transpose(magParamValsRow)
    paramsSolnNonlin = np.transpose(defaultParamValsOrigRow) + dparamsSolnNonlin
    if beVerbose:
        print("paramsSolnNonlin.T=", paramsSolnNonlin.T)
        print("normlzdSensMatrix@dnPS.x.T=", normlzdSensMatrix @ dnormlzdParamsSolnNonlin)
        print("normlzdDefaultBiasesCol.T=", normlzdDefaultBiasesCol.T)
        print("normlzdSensMatrix=", normlzdSensMatrix)

    normlzdWeightedDefaultBiasesApproxNonlin = \
             fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrix, normlzdCurvMatrix, numMetrics,
                    normlzdInteractDerivs, interactIdxs) \
             * metricsWeights

    ### WHAT IS THIS METRIC USED FOR???
    scale = 2
    normlzdWeightedDefaultBiasesApproxNonlin2x = \
             fwdFnc(scale * dnormlzdParamsSolnNonlin, normlzdSensMatrix, 1 * normlzdCurvMatrix, numMetrics,
                    normlzdInteractDerivs, interactIdxs) \
             * metricsWeights

    # Relationship between QuadTune variable names and math symbols:
    # defaultBiasesApproxNonlin = (       forward model soln       - default soln )
    #                           = ( f0 +      fwdFnc               - default soln )
    #                           = ( f0 + df/dp*dp + 0.5d2f/dp2*dp2 -       f0     )
    # residual = (   y_i -                y_hat_i                        )
    #          = (   y_i - ( f0    +   df/dp_i*dp + 0.5d2f/dp2_i*dp2 )   )
    #          =   ( y_i -   f0 )  - ( df/dp_i*dp + 0.5d2f/dp2_i*dp2 )
    #          = -defaultBiasesCol - (   defaultBiasesApproxNonlin   )
    #          = -defaultBiasesCol -              fwdFnc
    #          = normlzdResid * abs(normMetricValsCol)
    #  where f0 = defaultBiasesCol + obsMetricValsCol,
    #        y_i = obsMetricValsCol.
    #  globTunedBiases = forward global model soln - obs
    #                =                    -global_resid
    defaultBiasesApproxNonlin = normlzdWeightedDefaultBiasesApproxNonlin \
                                * np.reciprocal(metricsWeights) * np.abs(normMetricValsCol)

    defaultBiasesApproxNonlin2x = normlzdWeightedDefaultBiasesApproxNonlin2x \
                                * np.reciprocal(metricsWeights) * np.abs(normMetricValsCol)

    # To provide error bars, calculate solution with no nonlinear term and double the nonlinear term
    defaultBiasesApproxNonlinNoCurv = \
             fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrix, 0 * normlzdCurvMatrix, numMetrics,
                    normlzdInteractDerivs, interactIdxs) \
             * np.abs(normMetricValsCol)

    defaultBiasesApproxNonlin2xCurv = \
             fwdFnc(dnormlzdParamsSolnNonlin, normlzdSensMatrix, 2 * normlzdCurvMatrix, numMetrics,
                    normlzdInteractDerivs, interactIdxs) \
             * np.abs(normMetricValsCol)



    dnormlzdParamsSolnLin = minimize(lossFnc, x0=np.zeros_like(np.transpose(defaultParamValsOrigRow[0])),
                                     args=(normlzdSensMatrix, normlzdDefaultBiasesCol, metricsWeights,
                               0*normlzdCurvMatrix, reglrCoef, numMetrics),
                                     method='Powell')
    dnormlzdParamsSolnLin = np.atleast_2d(dnormlzdParamsSolnLin.x).T
    dparamsSolnLin = dnormlzdParamsSolnLin * np.transpose(magParamValsRow)
    paramsSolnLin = np.transpose(defaultParamValsOrigRow) + dparamsSolnLin


    return (defaultBiasesApproxNonlin,
            dnormlzdParamsSolnNonlin, paramsSolnNonlin,
            dnormlzdParamsSolnLin, paramsSolnLin,
            defaultBiasesApproxNonlin2x,
            defaultBiasesApproxNonlinNoCurv, defaultBiasesApproxNonlin2xCurv
           )


def constructNormlzdSensCurvMatrices(metricsNames, paramsNames, transformedParamsNames,
                               normMetricValsCol, magParamValsRow,
                               sens1NcFilenames, sens2NcFilenames, defaultNcFilename):
    """
    For nonlinear 2nd-order term of Taylor series: 0.5*dp^2*d2m/dp2+...,
    construct a numMetrics x numParams matrix of 2nd-order derivatives, d2m/dp2.
    Each row is a different metric.  Each column is a different parameter.
    The matrix is nondimensionalized by the observed values of metrics
    and maximum values of parameters.
    """


    from set_up_inputs import setupSensArrays
    from set_up_inputs import setupDefaultParamVectors, \
                              setUpDefaultMetricValsCol


    if ( len(paramsNames) != len(sens1NcFilenames)   ):
        print("Number of parameters does not equal number of netcdf files.")
        quit()

    # Number of tunable parameters
    numParams = len(paramsNames)

    # Number of metrics
    numMetrics = len(metricsNames)

    # For use in normalizing metrics matrices
    invrsObsMatrix = np.reciprocal(np.abs(normMetricValsCol)) @ np.ones((1,numParams))

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setUpDefaultMetricValsCol(metricsNames, defaultNcFilename)
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))
    normlzdDefaultMetricValsMatrix = defaultMetricValsMatrix * invrsObsMatrix

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultParamValsRow, defaultParamValsOrigRow = \
        setupDefaultParamVectors(paramsNames, transformedParamsNames,
                                 numParams,
                                 defaultNcFilename)
    normlzdDefaultParamValsRow = defaultParamValsRow * np.reciprocal(magParamValsRow)

    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sens1MetricValsMatrix, sens1ParamValsRow, sens1ParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sens1NcFilenames,
                        beVerbose=False)
    normlzdSens1ParamValsRow = sens1ParamValsRow * np.reciprocal(magParamValsRow)
    normlzdSens1MetricValsMatrix = sens1MetricValsMatrix * invrsObsMatrix

    # Set up sensitivity-simulation matrices from the extended sensitivity simulation
    sens2MetricValsMatrix, sens2ParamValsRow, sens2ParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sens2NcFilenames,
                        beVerbose=False)
    normlzdSens2ParamValsRow = sens2ParamValsRow * np.reciprocal(magParamValsRow)
    normlzdSens2MetricValsMatrix = sens2MetricValsMatrix * invrsObsMatrix

    # Initialize matrix to store second derivatives of metrics w.r.t. parameters
    normlzdCurvMatrix = np.zeros_like(sens1MetricValsMatrix)
    normlzdCurvMatrix2 = np.zeros_like(sens1MetricValsMatrix)  # 2nd way of calculating derivs
    normlzdCurvMatrixPoly = np.zeros_like(sens1MetricValsMatrix)  # 3rd way of calculating derivs
    normlzdSensMatrixPoly = np.zeros_like(sens1MetricValsMatrix)  # Approx of linear sensitivity
    normlzdConstMatrixPoly = np.zeros_like(sens1MetricValsMatrix)  # Approx of linear sensitivity
    normlzdOrdDparamsMin = np.zeros_like(sens1MetricValsMatrix)
    normlzdOrdDparamsMax = np.zeros_like(sens1MetricValsMatrix)

    #pdb.set_trace()

    # Compute quadratic coefficients using a polynomial fit to metric and parameters
    # normlzdOrdMetrics and normlzdOrdParams are length-3 python lists
    for arrayCol in np.arange(numParams):
        for arrayRow in np.arange(numMetrics):

            # Set up three (x,y) points whose 2nd-order derivative we wish to calculate.
            # For the spline code below, the x points need to be ordered from least to greatest.
            if normlzdSens1ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdSens1ParamValsRow[0,arrayCol],
                       normlzdDefaultParamValsRow[0,arrayCol],
                       normlzdSens2ParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdSens1MetricValsMatrix[arrayRow,arrayCol],
                        normlzdDefaultMetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens2MetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdSens2ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdSens2ParamValsRow[0,arrayCol],
                       normlzdDefaultParamValsRow[0,arrayCol],
                       normlzdSens1ParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdSens2MetricValsMatrix[arrayRow,arrayCol],
                        normlzdDefaultMetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens1MetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdDefaultParamValsRow[0,arrayCol],
                       normlzdSens1ParamValsRow[0,arrayCol],
                       normlzdSens2ParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdDefaultMetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens1MetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens2MetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdDefaultParamValsRow[0,arrayCol],
                       normlzdSens2ParamValsRow[0,arrayCol],
                       normlzdSens1ParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdDefaultMetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens2MetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens1MetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdSens1ParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdSens1ParamValsRow[0,arrayCol],
                       normlzdSens2ParamValsRow[0,arrayCol],
                       normlzdDefaultParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdSens1MetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens2MetricValsMatrix[arrayRow,arrayCol],
                        normlzdDefaultMetricValsMatrix[arrayRow,arrayCol] ]
            elif normlzdSens2ParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol]:
                normlzdOrdParams = [ normlzdSens2ParamValsRow[0,arrayCol],
                       normlzdSens1ParamValsRow[0,arrayCol],
                       normlzdDefaultParamValsRow[0,arrayCol] ]
                normlzdOrdDparams = normlzdOrdParams - normlzdDefaultParamValsRow[0,arrayCol]
                normlzdOrdMetrics = [ normlzdSens2MetricValsMatrix[arrayRow,arrayCol],
                        normlzdSens1MetricValsMatrix[arrayRow,arrayCol],
                        normlzdDefaultMetricValsMatrix[arrayRow,arrayCol] ]
            else:
                print("Error: Sensitivity parameter values are equal to each other or the default value in constructNormlzdSensCurvMatrices.")
                print( "normlzdSens1ParamValsRow=   ",
                      np.array2string(normlzdSens1ParamValsRow, formatter={'float_kind': lambda x: f"{x:12.8f}"}) )
                print( "normlzdSens2ParamValsRow=   ",
                      np.array2string(normlzdSens2ParamValsRow, formatter={'float_kind': lambda x: f"{x:12.8f}"}) )
                print( "normlzdDefaultParamValsRow= ",
                      np.array2string(normlzdDefaultParamValsRow, formatter={'float_kind': lambda x: f"{x:12.8f}"}) )
                quit()

            normlzdOrdDparamsMin[arrayRow,arrayCol] = np.min(normlzdOrdDparams)
            normlzdOrdDparamsMax[arrayRow,arrayCol]  = np.max(normlzdOrdDparams)

            # Calculate second-order spline based on three given (x,y) points.
            metricValsSpline = UnivariateSpline(normlzdOrdParams,normlzdOrdMetrics,s=0,k=2)
            # Based on spline, find 2nd derivative at arbitrary point (1).
            # I hope that the derivative has the same value at all points,
            #    since it is a parabola.
            normlzdCurvMatrix[arrayRow,arrayCol] = metricValsSpline.derivative(n=2)(1)

            # Check results using a 2nd calculation
            polyCoefs = np.polyfit(normlzdOrdParams, normlzdOrdMetrics, 2)
            normlzdCurvMatrix2[arrayRow,arrayCol] = 2*polyCoefs[0]
            #pdb.set_trace()

            # Check results using a 3rd calculation
            polyCoefs = np.polyfit(normlzdOrdDparams, normlzdOrdMetrics, 2)
            # The curvature matrix is d2m/dp2, not 0.5*d2m/dp2
            #     because of this 2:
            normlzdCurvMatrixPoly[arrayRow][arrayCol] = 2. * polyCoefs[0]
            normlzdSensMatrixPoly[arrayRow,arrayCol] = polyCoefs[1]
            normlzdConstMatrixPoly[arrayRow,arrayCol] = polyCoefs[2]

    # Read in information for piecewise linear forward function

    # Value of parameter between the high and low parameter values
    #     (usually, but not always, the default value)
    normlzd_pMid = np.zeros_like(defaultParamValsRow)
    # Sensitivity (slope) of metrics to the left of normlzd_pMid
    normlzdLeftSensMatrix = np.zeros_like(sens1MetricValsMatrix)
    # Sensitivity (slope) of metrics to the left of normlzd_pMid
    normlzdRightSensMatrix = np.zeros_like(sens1MetricValsMatrix)


    for arrayCol in np.arange(numParams):

        # normlzd_pMid, LeftSens, and RightSens depend on relative values of parameters
        if normlzdSens1ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol]:

            normlzd_pMid[0,arrayCol] = normlzdDefaultParamValsRow[0,arrayCol]
            normlzdLeftSensMatrix[arrayCol,:] = \
                ( normlzdDefaultMetricValsMatrix[arrayCol,:] - normlzdSens1MetricValsMatrix[arrayCol,:] ) \
                / ( normlzdDefaultParamValsRow[0,arrayCol] - normlzdSens1ParamValsRow[0,arrayCol] )
            normlzdRightSensMatrix[arrayCol,:] = \
                ( normlzdSens2MetricValsMatrix[arrayCol,:] - normlzdDefaultMetricValsMatrix[arrayCol,:] ) \
                / ( normlzdSens2ParamValsRow[0,arrayCol] - normlzdDefaultParamValsRow[0,arrayCol] )

        elif normlzdSens2ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol]:

            normlzd_pMid[0,arrayCol] = normlzdDefaultParamValsRow[0,arrayCol]
            normlzdLeftSensMatrix[arrayCol,:] = \
                ( normlzdDefaultMetricValsMatrix[arrayCol,:] - normlzdSens2MetricValsMatrix[arrayCol,:] ) \
                / ( normlzdDefaultParamValsRow[0,arrayCol] - normlzdSens2ParamValsRow[0,arrayCol] )
            normlzdRightSensMatrix[arrayCol,:] = \
                ( normlzdSens1MetricValsMatrix[arrayCol,:] - normlzdDefaultMetricValsMatrix[arrayCol,:] ) \
                / ( normlzdSens1ParamValsRow[0,arrayCol] - normlzdDefaultParamValsRow[0,arrayCol] )

        elif normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol]:

            normlzd_pMid[0,arrayCol] = normlzdSens1ParamValsRow[0,arrayCol]
            normlzdLeftSensMatrix[arrayCol,:] = \
                ( normlzdSens1MetricValsMatrix[arrayCol,:] - normlzdDefaultMetricValsMatrix[arrayCol,:] ) \
                / ( normlzdSens1ParamValsRow[0,arrayCol] - normlzdDefaultParamValsRow[0,arrayCol] )
            normlzdRightSensMatrix[arrayCol,:] = \
                ( normlzdSens2MetricValsMatrix[arrayCol,:] - normlzdSens1MetricValsMatrix[arrayCol,:] ) \
                / ( normlzdSens2ParamValsRow[0,arrayCol] - normlzdSens1ParamValsRow[0,arrayCol] )

        elif normlzdDefaultParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol]:

            normlzd_pMid[0,arrayCol] = normlzdSens2ParamValsRow[0,arrayCol]
            normlzdLeftSensMatrix[arrayCol,:] = \
                ( normlzdSens2MetricValsMatrix[arrayCol,:] - normlzdDefaultMetricValsMatrix[arrayCol,:] ) \
                / ( normlzdSens2ParamValsRow[0,arrayCol] - normlzdDefaultParamValsRow[0,arrayCol] )
            normlzdRightSensMatrix[arrayCol,:] = \
                ( normlzdSens1MetricValsMatrix[arrayCol,:] - normlzdSens2MetricValsMatrix[arrayCol,:] ) \
                / ( normlzdSens1ParamValsRow[0,arrayCol] - normlzdSens2ParamValsRow[0,arrayCol] )

        elif normlzdSens1ParamValsRow[0,arrayCol] < normlzdSens2ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol]:

            normlzd_pMid[0,arrayCol] = normlzdSens2ParamValsRow[0,arrayCol]
            normlzdLeftSensMatrix[arrayCol,:] = \
                ( normlzdSens2MetricValsMatrix[arrayCol,:] - normlzdSens1MetricValsMatrix[arrayCol,:] ) \
                / ( normlzdSens2ParamValsRow[0,arrayCol] - normlzdSens1ParamValsRow[0,arrayCol] )
            normlzdRightSensMatrix[arrayCol,:] = \
                ( normlzdDefaultMetricValsMatrix[arrayCol,:] - normlzdSens2MetricValsMatrix[arrayCol,:] ) \
                / ( normlzdDefaultParamValsRow[0,arrayCol] - normlzdSens2ParamValsRow[0,arrayCol] )

        elif normlzdSens2ParamValsRow[0,arrayCol] < normlzdSens1ParamValsRow[0,arrayCol] < normlzdDefaultParamValsRow[0,arrayCol]:

            normlzd_pMid[0,arrayCol] = normlzdSens1ParamValsRow[0,arrayCol]
            normlzdLeftSensMatrix[arrayCol,:] = \
                ( normlzdSens1MetricValsMatrix[arrayCol,:] - normlzdSens2MetricValsMatrix[arrayCol,:] ) \
                / ( normlzdSens1ParamValsRow[0,arrayCol] - normlzdSens2ParamValsRow[0,arrayCol] )
            normlzdRightSensMatrix[arrayCol,:] = \
                ( normlzdDefaultMetricValsMatrix[arrayCol,:] - normlzdSens1MetricValsMatrix[arrayCol,:] ) \
                / ( normlzdDefaultParamValsRow[0,arrayCol] - normlzdSens1ParamValsRow[0,arrayCol] )

        else:
            print("Error: Sensitivity parameter values are equal to each other or the default value in constructNormlzdSensCurvMatrices.")
            print( "normlzdSens1ParamValsRow=   ",
                  np.array2string(normlzdSens1ParamValsRow, formatter={'float_kind': lambda x: f"{x:12.8f}"}) )
            print( "normlzdSens2ParamValsRow=   ",
                  np.array2string(normlzdSens2ParamValsRow, formatter={'float_kind': lambda x: f"{x:12.8f}"}) )
            print( "normlzdDefaultParamValsRow= ",
                  np.array2string(normlzdDefaultParamValsRow, formatter={'float_kind': lambda x: f"{x:12.8f}"}) )
            quit()


    return ( normlzdCurvMatrixPoly, normlzdSensMatrixPoly, normlzdConstMatrixPoly,
             normlzdOrdDparamsMin, normlzdOrdDparamsMax,
             normlzd_pMid, normlzdLeftSensMatrix, normlzdRightSensMatrix)


def approxMatrixWithSvd( matrix , sValsRatio, sValsNumToKeep,
                         beVerbose):
    """
    Input: A matrix
    Output: A possibly lower-rank approximation of the matrix,
            with a max ratio of singular values specified by sValsRatio.
    """

    # vh = V^T = transpose of right-singular vector matrix, V.
    #  matrix = u @ np.diag(sVals) @ vh = (u * sVals) @ vh
    u, sVals, vh = np.linalg.svd( matrix, full_matrices=False )

    # Delete the small singular values in order to show just the most important patterns.
    # After this deletion, store inverse singular values in sValsTrunc
    sValsTrunc = np.copy(sVals)
    if (sValsNumToKeep == None):
        for idx, sVal in np.ndenumerate(sVals):
            # If a singular value is much smaller than largest singular value,
            #     then zero it out.
            if np.divide(sVals[0],np.maximum(sVal,np.finfo(float).eps)) > sValsRatio:
                sValsTrunc[idx] = 0.
    else:
        for idx, sVal in np.ndenumerate(sVals):
            if idx+1 > sValsNumToKeep:
                sValsTrunc[idx] = 0.

    if beVerbose:
        print("\nOriginal singular values =")
        print(sVals)
        print("\nsValsTrunc =")
        print(sValsTrunc)

    matrixApprox = u @ np.diag(sValsTrunc) @ vh
    #matrixApprox = (u * sVals) @ vh

    if beVerbose:
        print("\nstd/mean of (matrixApprox-matrix) = ")
        print(np.std(np.subtract(matrixApprox, matrix))/np.mean(matrix))

    return matrixApprox

#def constructNormlzd2ndOrderTensor(numParams, numMetrics,
#                                   normlzdCurvMatrix, normlzdSensMatrixPoly, normlzdConstMatrix, corrParam):
#    '''Constructs a numParams x numMetrics x numParams tensor of second derivatives, d2metrics/dparams1dparams2.
#    The derivatives are normalized by observed metric values and max param values.'''
#
#    normlzdCurvMatrix, normlzdSensMatrixPoly, normlzdConstMatrix, \
#    normlzdOrdDparamsMin, normlzdOrdDparamsMax = \
#        constructNormlzdSensCurvMatrices(metricsNames, paramsNames, transformedParamsNames, \
#                                   metricsWeights, obsMetricValsCol, normMetricValsCol, magParamValsRow, \
#                                   sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)
#
#    for i in np.arange(numMetrics):
#        # Construct matrix with M(j,k)=0 if curv opposite for i,j; otherwise, M(j,k)=1.
#        curvParamMatrix = np.outer( normlzdCurvMatrix[i,:], np.ones(numParams) )
#        corrParamMatrix = 0.5 * np.abs( np.sign( curvParamMatrix ) + np.sign( curvParamMatrix.T ) )
#        corrParamMatrix = np.fill_diagonal(corrParamMatrix, 1.0)
#        diagSqrtCurvMatrix = np.diag(np.sqrt(normlzdCurvMatrix[i,:]))
#
#    # d4 = np.einsum('j,jik,k->i', params, offdiagtensor, params)
#
#    return ( normlzd2ndOrderTensor  )

def calcNormlzdRadiusCurv(metricsNames, paramsNames, transformedParamsNames, paramsScales,
                          metricsWeights, obsMetricValsCol,
                          sensNcFilenames, sensNcFilenamesExt, defaultNcFilename):
    """
    Calculate radius of curvature of output from 2 sensitivity simulations plus the default
    simulation.
    """

    from set_up_inputs import setupDefaultParamVectors, \
                                           setupSensArrays
    from set_up_inputs import setUpDefaultMetricValsCol

    if ( len(paramsNames) != len(sensNcFilenames)   ):
        print("Number of parameters must equal number of netcdf files.")
        quit()

    # Number of tunable parameters
    numParams = len(paramsNames)

    # Number of metrics
    numMetrics = len(metricsNames)

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setUpDefaultMetricValsCol(metricsNames, defaultNcFilename)

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultParamValsRow, defaultParamValsOrigRow = \
        setupDefaultParamVectors(paramsNames, transformedParamsNames,
                                 numParams,
                                 defaultNcFilename)

    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sensMetricValsMatrix, sensParamValsRow, sensParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sensNcFilenames,
                        beVerbose=False)

    # Set up sensitivity-simulation matrices from the extended sensitivity simulation
    sensMetricValsMatrixExt, sensParamValsRowExt, sensParamValsOrigRowExt = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sensNcFilenamesExt,
                        beVerbose=False)

    normlzd_radius_of_curv = np.full_like(sensMetricValsMatrix, 0.0)

    # Calculate differences in parameter values between default, sensitivity,
    #    and extended sensitivity runs.
    delta_params_def_sens = sensParamValsRow - defaultParamValsRow
    delta_params_def_sensExt = sensParamValsRowExt - defaultParamValsRow
    delta_params_sens_sensExt = sensParamValsRowExt - sensParamValsRow

    # Calculate numMetrics x numParams matrix of metric values.
    delta_metrics_def_sens = sensMetricValsMatrix - defaultMetricValsMatrix
    delta_metrics_def_sensExt = sensMetricValsMatrixExt - defaultMetricValsMatrix
    delta_metrics_sens_sensExt = sensMetricValsMatrixExt - sensMetricValsMatrix
    for col in np.arange(numParams):
        for row in np.arange(numMetrics):
            # Distance between points in simulations = sqrt(dparam**2 + dmetric**2)
            length_def_sens = np.linalg.norm([delta_params_def_sens[0][col],
                                              delta_metrics_def_sens[row][col]])
            length_def_sensExt = np.linalg.norm([delta_params_def_sensExt[0][col],
                                                 delta_metrics_def_sensExt[row][col]])
            length_sens_sensExt = np.linalg.norm([delta_params_sens_sensExt[0][col],
                                                  delta_metrics_sens_sensExt[row][col]])
            semi_perim = 0.5 * ( length_def_sens + length_def_sensExt + length_sens_sensExt )
            # area of triangle formed by points.  Use Heron's formula.
            area = np.sqrt( semi_perim *
                           (semi_perim-length_def_sens) *
                           (semi_perim-length_def_sensExt) *
                           (semi_perim-length_sens_sensExt)
                          )
            if (area == 0.0):
                print( '\nIn calcNormlzdRadiusCurv, area == 0.0 for param ', paramsNames[col],
                        'and metric ', metricsNames[row] )

            # Greatest distance between parameter values in the 3 simulations:
            max_params_width = \
            np.max(np.abs([delta_params_def_sens[0][col],
                        delta_params_def_sensExt[0][col],
                        delta_params_sens_sensExt[0][col]]))
            if (max_params_width == 0.0):
                print( '\nIn calcNormlzdRadiusCurv, max_params_width == 0.0 for param ', paramsNames[col],
                        'and metric ', metricsNames[row] )

            # Calculate Menger curvature from triangle area and distance between points:
            normlzd_radius_of_curv[row][col] = 0.25 * length_def_sens*length_def_sensExt*length_sens_sensExt \
                                                / area / max_params_width

    #pdb.set_trace()
    fig, axs = plt.subplots(numMetrics, numParams, figsize=(24,36))
    for col in np.arange(numParams):
        for row in np.arange(numMetrics):

            paramVals = [defaultParamValsRow[0][col], sensParamValsRow[0][col], sensParamValsRowExt[0][col]]
            metricVals = [defaultMetricValsMatrix[row][col], sensMetricValsMatrix[row][col],
                  sensMetricValsMatrixExt[row][col]]

            axs[row, col].plot( paramVals, metricVals, marker=".", ls="" )
            axs[row, col].plot( paramVals, obsMetricValsCol[row][0] * np.ones((3,1)), color="r" )
            axs[row, col].set_xlabel(paramsNames[col])
            axs[row, col].set_ylabel(metricsNames[row])
            #fig.show()

    plt.show()
    plt.savefig('param_metric_scatter.png')
    #pdb.set_trace()

    return

def findOutliers(normlzdSensMatrix, normlzdWeightedSensMatrix,
                 defaultBiasesCol, normMetricValsCol, magParamValsRow, defaultParamValsOrigRow):
    """Find outliers in bias-senstivity scatterplot based on the RANSAC method."""




    #    ransac = linear_model.RANSACRegressor(max_trials=1000,random_state=0,
#                                          base_estimator=linear_model.LinearRegression(), residual_threshold=None)
#    ransac.fit(normlzdSensMatrix, -defaultBiasesCol / np.abs(normMetricValsCol) )
#    defaultBiasesApproxRansac = ransac.predict(normlzdSensMatrix) * np.abs(normMetricValsCol)
#    dnormlzdParamsSolnRansac = np.transpose( ransac.estimator_.coef_ )
#    ransac = linear_model.RANSACRegressor(max_trials=10000,random_state=0,
#             base_estimator=linear_model.ElasticNet(fit_intercept=False, random_state=0, tol=1e-3,
#                                                    l1_ratio=0.0, alpha=5),
#                                          residual_threshold=None)
#    ransac.fit(normlzdSensMatrix, -defaultBiasesCol / np.abs(normMetricValsCol) )
#    defaultBiasesApproxRansac = np.transpose( np.atleast_2d( ransac.predict(normlzdSensMatrix) ) ) \
#                                * np.abs(normMetricValsCol)
#    dnormlzdParamsSolnRansac = np.transpose( np.atleast_2d( ransac.estimator_.coef_ ) )
#    inlier_mask = ransac.inlier_mask_
#    outlier_mask = np.logical_not(inlier_mask)


    ransac = linear_model.HuberRegressor(fit_intercept=False)
    ransac.fit(normlzdSensMatrix, -defaultBiasesCol / np.abs(normMetricValsCol) )
    defaultBiasesApproxRansac = np.transpose( np.atleast_2d( ransac.predict(normlzdSensMatrix) ) ) \
                                * np.abs(normMetricValsCol)
    dnormlzdParamsSolnRansac = np.transpose( np.atleast_2d( ransac.coef_ ) )
#    defaultBiasesApproxRansac = ransac.predict(normlzdSensMatrix) * np.abs(normMetricValsCol)
#    dnormlzdParamsSolnRansac = np.transpose( ransac.coef_ )

    dparamsSolnRansac = dnormlzdParamsSolnRansac * np.transpose(magParamValsRow)
    paramsSolnRansac = np.transpose(defaultParamValsOrigRow) + dparamsSolnRansac

    outlier_mask = ransac.outliers_
    inlier_mask = np.logical_not(outlier_mask)

    print( "paramsSolnRansac = ", paramsSolnRansac )
    print( "dparamsSolnRansac = ", dparamsSolnRansac )

    # If the solution were perfect, this variable would equal
    #     the normalized, weighted right-hand side.
    normlzdWeightedDefaultBiasesApproxRansac = \
            normlzdWeightedSensMatrix @ dnormlzdParamsSolnRansac

    #pdb.set_trace()

    return (outlier_mask, defaultBiasesApproxRansac, normlzdWeightedDefaultBiasesApproxRansac,
            dnormlzdParamsSolnRansac, paramsSolnRansac)



def findParamsUsingElastic(normlzdSensMatrix, normlzdWeightedSensMatrix,
                 defaultBiasesCol, normMetricValsCol, metricsWeights,
                 magParamValsRow, defaultParamValsOrigRow,
                 normlzdCurvMatrix,
                 beVerbose):
    """Do linear regression with L1 (lasso or elastic net) regularization"""



    #regr = ElasticNet(fit_intercept=True, random_state=0, tol=1e-10, l1_ratio=0.5, alpha=0.01)
    #regr =linear_model.Lasso(fit_intercept=True, random_state=0, tol=1e-10, alpha=0.01) # don't fit intercept!;use line below
    regr = linear_model.Lasso(fit_intercept=False, random_state=0, tol=1e-10, alpha=0.01)
    #regr = linear_model.LassoCV(fit_intercept=True, random_state=0, eps=1e-5, tol=1e-10, cv=metricsWeights.size)
    #print( "alpha_ = ", regr.alpha_ )
    regr.fit(normlzdWeightedSensMatrix, -metricsWeights * defaultBiasesCol / np.abs(normMetricValsCol) )
    #regr.fit(normlzdSensMatrix, -defaultBiasesCol / np.abs(normMetricValsCol) )

    defaultBiasesApproxElastic = np.transpose( np.atleast_2d(regr.predict(normlzdWeightedSensMatrix)) ) \
                          * np.reciprocal(metricsWeights) * np.abs(normMetricValsCol)
    dnormlzdParamsSolnElastic = np.transpose( np.atleast_2d(regr.coef_) )
    dparamsSolnElastic = dnormlzdParamsSolnElastic * np.transpose(magParamValsRow)
    paramsSolnElastic = np.transpose(defaultParamValsOrigRow) + dparamsSolnElastic

    if beVerbose:
        print( "paramsSolnElastic = ", paramsSolnElastic )
        print( "dparamsSolnElastic = ", dparamsSolnElastic )

    #pdb.set_trace()

    # If the solution were perfect, this variable would equal
    #     the normalized, weighted right-hand side.
    defaultBiasesApproxElasticNonlin = \
            normlzdWeightedSensMatrix @ dnormlzdParamsSolnElastic \
                        * np.reciprocal(metricsWeights) * np.abs(normMetricValsCol) \
            + 0.5 * normlzdCurvMatrix @ (dnormlzdParamsSolnElastic**2) * np.abs(normMetricValsCol)

    #pdb.set_trace()

    return (defaultBiasesApproxElastic, defaultBiasesApproxElasticNonlin,
            dnormlzdParamsSolnElastic, paramsSolnElastic)



if __name__ == '__main__':
    main()
#        sensMatrixDashboard.run_server(debug=True)
