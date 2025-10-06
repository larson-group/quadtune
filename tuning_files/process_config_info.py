import numpy as np
import pandas as pd

from set_up_inputs import setUp_x_MetricsList, setUpDefaultMetricValsCol

def process_paramsnames_scales_and_filesuffixes(paramsNamesScalesAndFilenames,folder_name):
    """
    Split up paramsNamesScalesAndFilenames into: \n 
    - paramsNames
    - paramsScales
    - sensNcFilenames
    - sensNcFilenamesExt
    """


    dfparamsNamesScalesAndFilenames = \
        pd.DataFrame(paramsNamesScalesAndFilenames,
                    columns=['paramsNames', 'paramsScales',
                            'sensNcFilenamesSuffix', 'sensNcFilenamesSuffixExt'])
    paramsNames = dfparamsNamesScalesAndFilenames[['paramsNames']].to_numpy().astype(str)[:, 0]

    paramsScales = dfparamsNamesScalesAndFilenames[['paramsScales']].to_numpy().astype(float)[:, 0]

    sensNcFilenamesSuffix = dfparamsNamesScalesAndFilenames[['sensNcFilenamesSuffix']].to_numpy().astype(str)[:, 0]

    sensNcFilenames = np.char.add(folder_name, sensNcFilenamesSuffix)

    

    sensNcFilenamesSuffixExt = dfparamsNamesScalesAndFilenames[['sensNcFilenamesSuffixExt']].to_numpy().astype(str)[:, 0]

    sensNcFilenamesExt = np.char.add(folder_name, sensNcFilenamesSuffixExt)

    

    return (paramsNames, paramsScales, sensNcFilenames, sensNcFilenamesExt)



def process_prescribed_paramsnames(prescribedParamsNamesScalesAndValues,folder_name):
    """
    Split up prescribedParamsNamesScalesAndValues into: 
    - prescribedParamsNames
    - prescribedParamsScales
    - prescribedParamValsRow
    - prescribedSensNcFilenames
    - prescribedSensNcFilenamesExt
    - prescribedTransformedParamsNames

    """

    dfprescribedParamsNamesScalesAndValues = \
        pd.DataFrame(prescribedParamsNamesScalesAndValues,
                     columns=['prescribedParamsNames',
                              'prescribedParamsScales',
                              'prescribedParamVals',
                              'prescribedSensNcFilenamesSuffix', 'prescribedSensNcFilenamesSuffixExt'
                              ]
                     )
    prescribedParamsNames = dfprescribedParamsNamesScalesAndValues[['prescribedParamsNames']].to_numpy().astype(str)[:,0]


    # Extract scaling factors of parameter values from user-defined list paramsNamesScalesAndFilenames.
    # The scaling is not used for any calculations, but it allows us to avoid plotting very large or small values.
    prescribedParamsScales = dfprescribedParamsNamesScalesAndValues[['prescribedParamsScales']].to_numpy().astype(float)[:, 0]

    prescribedParamVals = dfprescribedParamsNamesScalesAndValues[['prescribedParamVals']].to_numpy().astype(float)[:, 0]
    prescribedParamValsRow = prescribedParamVals


    prescribedSensNcFilenamesSuffix = dfprescribedParamsNamesScalesAndValues[
                                    ['prescribedSensNcFilenamesSuffix']].to_numpy().astype(str)[:, 0]
    prescribedSensNcFilenames = np.char.add(folder_name, prescribedSensNcFilenamesSuffix)

    prescribedSensNcFilenamesSuffixExt = dfprescribedParamsNamesScalesAndValues[
                                       ['prescribedSensNcFilenamesSuffixExt']].to_numpy().astype(str)[:, 0]
    prescribedSensNcFilenamesExt = np.char.add(folder_name, prescribedSensNcFilenamesSuffixExt)

    prescribedTransformedParamsNames = np.array([''])


    return (prescribedParamsNames, prescribedParamsScales, prescribedParamValsRow, prescribedSensNcFilenames, prescribedSensNcFilenamesExt, prescribedTransformedParamsNames)





def process_metrics_names_weights_norms_custom(metricsNamesWeightsAndNormsCustom,metricsNames,metricsWeights,metricsNorms):
    """
    Split up metricsNamesWeightsAndNormsCustom into: 
    - metricsNames
    - metricsWeights
    - metricsNorms

    """


    # Split metricsNamesWeightsAndNormsCustom into metric names and the corresponding weights.
    dfMetricsNamesWeightsAndNormsCustom = \
        pd.DataFrame(metricsNamesWeightsAndNormsCustom,
                     columns=['metricsNamesCustom', 'metricsWeightsCustom', 'metricsNormsCustom'])
    metricsNamesCustom = dfMetricsNamesWeightsAndNormsCustom[['metricsNamesCustom']].to_numpy().astype(str)[:, 0]
    metricsWeightsCustom = dfMetricsNamesWeightsAndNormsCustom[['metricsWeightsCustom']].to_numpy().astype(float)
    metricsNormsCustom = dfMetricsNamesWeightsAndNormsCustom[['metricsNormsCustom']].to_numpy().astype(float)

    # Include custom regions in metricsNames:
    metricsNames = np.append(metricsNames, metricsNamesCustom)

    metricsWeights = np.vstack((metricsWeights, metricsWeightsCustom))


    #numMetricsCustom = len(metricsNames) - numMetricsNoCustom
    metricsNorms = np.vstack((metricsNorms, metricsNormsCustom))


    #metricsNamesNoprefix = np.char.replace(metricsNames, "SWCF_", "")
    metricsNamesNoprefix = metricsNames

    return (metricsNames, metricsWeights, metricsNorms, metricsNamesNoprefix)


def process_metrics_names_weights_norms(defaultNcFilename,varPrefixes):


    metricsNamesWeightsAndNorms, metricGlobalValsFromFile \
        = setUp_x_MetricsList(varPrefixes, defaultNcFilename)

    dfMetricsNamesWeightsAndNorms = \
        pd.DataFrame(metricsNamesWeightsAndNorms, columns=['metricsNames', 'metricsWeights', 'metricsNorms'])
    
    # Here we initialize metricsNames to all the standard (e.g., 20x20) regions for all varPrefixes
    #    Later down, we'll append any custom regions, e.g., DYCOMS
    metricsNames = dfMetricsNamesWeightsAndNorms[['metricsNames']].to_numpy().astype(str)[:, 0]
    metricsWeights = dfMetricsNamesWeightsAndNorms[['metricsWeights']].to_numpy().astype(float)



    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setUpDefaultMetricValsCol(metricsNames, defaultNcFilename)

    metricGlobalAvgs = np.diag(np.dot(metricsWeights.reshape(-1, len(varPrefixes), order='F').T,
                                      defaultMetricValsCol.reshape(-1, len(varPrefixes), order='F')))
    

    # Any special "custom" regions, e.g. DYCOMS, will be tacked onto the end of
    # the usual metrics vectors.  But we exclude those regions from numMetricsNoCustom.

    numMetricsNoCustom = len(metricsNames)

    return metricsNames, metricsWeights, metricGlobalAvgs, numMetricsNoCustom









