import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


def get_folder_name(numMetrics):
    # Determine the subfolder name based on numMetrics
    folder_name = ""
    if numMetrics == 162:
        folder_name = "20x20"
    elif numMetrics == 72:
        folder_name = "30x30"
    elif numMetrics == 32:
        folder_name = "45x45"
    else:
        folder_name = "unknownGrid"
    return folder_name


def create_folder(folder_name):
    # Create the 'Outputs' folder if it doesn't exist
    outputs_folder = "Outputs"
    os.makedirs(outputs_folder, exist_ok=True)

    # Create the subfolder
    subfolder_path = os.path.join(outputs_folder, folder_name)
    os.makedirs(subfolder_path, exist_ok=True)

    # Create inner folders for biases, losses, and params
    inner_folders = ["Residuals", "Parameters"]
    for folder in inner_folders:
        inner_folder_path = os.path.join(subfolder_path, folder)
        os.makedirs(inner_folder_path, exist_ok=True)


def bootstrapPlots(numMetricsToTune,
                   metricsNames,
                   residualsBootstrapMatrix,
                   residualsTunedCol,
                   residualsDefaultCol,
                   paramsNames,
                   paramsBoot,
                   paramsTuned,
                   defaultParamValsOrigRow,
                   paramBoundsBoot,
                   normResidualPairsMatrix,
                   tradeoffBinaryMatrix):
    """
    Generates and saves a series of bootstrap-based diagnostic plots for parameter tuning and residual analysis.

    Parameters:
        numMetricsToTune (int): Number of metrics selected for tuning.
        metricsNames (np.ndarray): Names of all metrics.
        residualsBootstrapMatrix (np.ndarray): Bootstrap residuals, with shape (n_samples, n_metrics).
        residualsTunedCol (np.ndarray): Residuals from tuned parameters, with shape (n_metrics, 1).
        residualsDefaultCol (np.ndarray): Residuals from default parameters, with shape (n_metrics, 1).
        paramsNames (np.ndarray): Names of the parameters.
        paramsBoot (np.ndarray): Bootstrap samples of parameters, with shape (n_samples, n_params).
        paramsTuned (np.ndarray): Tuned parameter values, with shape (n_params,).
        defaultParamValsOrigRow (np.ndarray): Original default parameter values, with shape (1, n_params).
        paramBoundsBoot (np.ndarray): Bootstrap confidence interval bounds,
                                      with shape (2, 3, n_params): [percentile/BCa] × [lower, tuned, upper] × parameters.
        normResidualPairsMatrix (np.ndarray): Pairwise min norm matrix of residuals, with shape (n_metrics, n_metrics).
        tradeoffBinaryMatrix (np.ndarray): Binary matrix indicating trade-offs, with shape (n_metrics, n_metrics).
    """
    print("Creating bootstrap plots . . .")
    # create folders to save in
    folderName = get_folder_name(numMetricsToTune)
    create_folder(folderName)

    plot_param_distributions(paramsBoot, paramsNames, defaultParamValsOrigRow, paramsTuned, folderName)
    plot_confidence_intervals(paramsNames, paramBoundsBoot, folderName)
    plot_residuals_distributions(residualsBootstrapMatrix, residualsDefaultCol, residualsTunedCol, metricsNames,
                                 folderName)
    plot_regional_variances(residualsBootstrapMatrix, folderName)
    plot_regional_msr(residualsBootstrapMatrix, folderName)
    plot_residual_map(normResidualPairsMatrix, metricsNames, False, "Residual_norm_heatmap", folderName)
    plot_residual_map(tradeoffBinaryMatrix, metricsNames, True, "Tradeoff_map", folderName)
    plot_tradeoff_scatter(tradeoffBinaryMatrix, residualsBootstrapMatrix, metricsNames, folderName)


def plot_param_distributions(paramsBoot, paramsNames, defaultParamValsOrigRow, paramsTuned, folderName):
    """
    Plots and saves the bootstrap distributions of each parameter.

    Parameters:
        paramsBoot (np.ndarray): Bootstrap samples of parameters, with shape (n_samples, n_params).
        paramsNames (np.ndarray): Names of the parameters.
        defaultParamValsOrigRow (np.ndarray): Original default parameter values, with shape (1, n_params).
        paramsTuned (np.ndarray): Tuned parameter values, with shape (n_params,).
        folderName (str): Subdirectory within "Outputs" where the plots will be saved.
    """
    for i, paramName in enumerate(paramsNames):
        paramSample = paramsBoot[:, i].flatten()

        fig, ax = plt.subplots(figsize=(6, 4))
        sns.histplot(paramSample, kde=True, stat="density", ax=ax)

        ax.axvline(defaultParamValsOrigRow[0, i], color='red', linestyle='--', linewidth=2, label='Default')
        ax.axvline(paramsTuned[i], color='green', linestyle='--', linewidth=2, label='Tuned')

        ax.set_xlabel(paramName)
        ax.set_title(f"{paramName}: mean = {np.mean(paramSample):.2g}, std = {np.std(paramSample):.2g}")
        ax.legend()
        ax.grid(True, linestyle=':', alpha=0.6)

        plt.tight_layout()
        plt.savefig(f"Outputs/{folderName}/Parameters/Distribution_{paramName}.png", dpi=300)
        plt.close("all")


def plot_confidence_intervals(paramsNames, paramBoundsBoot, folderName):
    """
    Plots and saves BCa and Percentile confidence intervals for tuned parameters.

    Parameters:
        paramsNames (np.ndarray): Names of the parameters.
        paramBoundsBoot (np.ndarray): Bootstrap confidence interval bounds,
                                      with shape (2, 3, n_params) where axis 0 indexes [Percentile, BCa],
                                      axis 1 indexes [lower, tuned, upper], and axis 2 is over parameters.
        folderName (str): Subdirectory within "Outputs" where the plot will be saved.
    """
    # Get intervals
    bcaLower, bcaUpper = paramBoundsBoot[1, 0], paramBoundsBoot[1, 2]
    percLower, percUpper = paramBoundsBoot[0, 0], paramBoundsBoot[0, 2]
    
    # Get tuned parameter values
    paramsTuned = paramBoundsBoot[0, 1]

    # Plot
    fig, axes = plt.subplots(len(paramsTuned), 1, figsize=(8, 2.2 * len(paramsTuned)), sharex=False)

    for i, ax in enumerate(axes):
        # Plot BCa interval
        ax.errorbar(paramsTuned[i], 0.2,
                    xerr=[[paramsTuned[i] - bcaLower[i]], [bcaUpper[i] - paramsTuned[i]]],
                    fmt='o', capsize=4, color='C0', label='BCa')

        # Plot Percentile interval
        ax.errorbar(paramsTuned[i], -0.2,
                    xerr=[[paramsTuned[i] - percLower[i]], [percUpper[i] - paramsTuned[i]]],
                    fmt='o', capsize=4, color='C1', label='Percentile')

        # Decorations
        ax.axvline(paramsTuned[i], color='green', linestyle='--', linewidth=0.8)
        ax.set_yticks([-0.2, 0.2])
        ax.set_yticklabels(['Percentile', 'BCa'])
        ax.set_title(paramsNames[i])
        ax.set_xlabel("Parameter Value")
        ax.grid(True, axis='x', linestyle=':', alpha=0.5)

    plt.tight_layout()
    plt.savefig(f"Outputs/{folderName}/Parameters/ErrorBars.png", dpi=300)
    plt.close("all")


def plot_residuals_distributions(residualsBootstrapMatrix, residualsDefaultCol, residualsTunedCol, metricsNames,
                                 folderName):
    """
    Plots and saves the distribution of bootstrap residuals for each metric.

    Parameters:
        residualsBootstrapMatrix (np.ndarray): Bootstrap residuals with shape (n_samples, n_metrics).
        residualsDefaultCol (np.ndarray): Residuals from default parameters, with shape (n_metrics, 1).
        residualsTunedCol (np.ndarray): Residuals from tuned parameters, with shape (n_metrics, 1).
        metricsNames (np.ndarray): Names of each metric.
        folderName (str): Subdirectory within "Outputs" where the plots will be saved.
    """
    for i in range(len(metricsNames)):
        residualsSample = residualsBootstrapMatrix[:, i]
        plt.figure(figsize=(12, 4))
        sns.displot(data=residualsSample, kde=True, stat="density", height=4, aspect=2.5)
        plt.axvline(residualsDefaultCol[i], color='red', linestyle='--', label="default residuals")
        plt.axvline(residualsTunedCol[i], color='green', linestyle='--', label="tuned residuals")
        plt.axvline(0, color='black', linestyle='-', linewidth=0.8, alpha=0.7)
        plt.xlabel(metricsNames[i])
        plt.title(f"Mean: {np.mean(residualsSample):.2g}, Variance: {np.var(residualsSample):.2g}")
        plt.tight_layout()
        plt.savefig(f"Outputs/{folderName}/Residuals/Distribution_{metricsNames[i]}.png", dpi=300)
        plt.close("all")


def plot_regional_variances(residualsBootstrapMatrix, folderName):
    """
    Computes and saves a plot of the variance of residuals per region.

    Parameters:
        residualsBootstrapMatrix (np.ndarray): Bootstrap residuals with shape (n_samples, n_metrics).
        folderName (str): Subdirectory within "Outputs" where the plot will be saved.
    """
    variances = np.var(residualsBootstrapMatrix, axis=0).reshape(9, 18)

    plt.figure(figsize=(12, 6))
    plt.imshow(variances, aspect='equal')
    plt.colorbar(label='Variance')

    # Set x-ticks to start at 1 instead of 0
    plt.xticks(ticks=np.arange(18), labels=np.arange(1, 19))
    plt.yticks(ticks=np.arange(9), labels=np.arange(1, 10))

    plt.tight_layout()
    plt.savefig(f"Outputs/{folderName}/Residuals/Variances.png", dpi=300)
    plt.close("all")


def plot_regional_msr(residualsBootstrapMatrix, folderName):
    """
    Computes and saves a plot of the Mean Squared Residual (MSR) per region.

    Parameters:
        residualsBootstrapMatrix (np.ndarray): Bootstrap residuals with shape (n_samples, n_metrics).
        folderName (str): Subdirectory within "Outputs" where the plot will be saved.
    """
    msr = np.mean(np.power(residualsBootstrapMatrix, 2), axis=0).reshape(9, 18)

    plt.figure(figsize=(12, 6))
    plt.imshow(msr, aspect='equal')
    plt.colorbar(label='MSR')

    # Set x-ticks to start at 1 instead of 0
    plt.xticks(ticks=np.arange(18), labels=np.arange(1, 19))
    plt.yticks(ticks=np.arange(9), labels=np.arange(1, 10))

    plt.tight_layout()
    plt.savefig(f"Outputs/{folderName}/Residuals/MSR.png", dpi=300)
    plt.close("all")


def plot_residual_map(matrix, metricsNames, binary, plotName, folderName):
    """
    Creates and saves a heatmap or binary map of a residual matrix.

    Parameters:
        matrix (np.ndarray): A square matrix (n_metrics x n_metrics) (e.g., of correlations, trade-offs, or binary flags) to visualize.
        metricsNames (np.ndarray): List of metric names used to label axes.
        binary (bool): If True, use a binary grayscale colormap ("Greys"); otherwise, use default colormap.
        plotName (str): Filename (without extension) for saving the plot.
        folderName (str): Subdirectory within "Outputs" where the plot will be saved.
    """
    n = matrix.shape[0]
    # Plot the matrix
    tickPositions = np.arange(len(metricsNames))
    tickLabels = [metricsNames[i] if i % 9 == 0 else "" for i in tickPositions]

    fig, ax = plt.subplots(figsize=(14, 14))
    if binary:
        ax.imshow(matrix, cmap='Greys', interpolation='none')
    else:
        ax.imshow(matrix)

    # Set ticks and partial labels
    ax.set_xticks(tickPositions)
    ax.set_yticks(tickPositions)
    ax.set_xticklabels(tickLabels, rotation=90, fontsize=12)
    ax.set_yticklabels(tickLabels, fontsize=12)

    # Draw light grid lines to help trace
    ax.set_xticks(np.arange(-.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n, 1), minor=True)
    ax.grid(which='minor', color='lightgray', linewidth=0.3)

    # Optional: thicker ticks every 9th label
    for i in range(0, n, 9):
        ax.axvline(i - 0.5, color='gray', linewidth=0.5)
        ax.axhline(i - 0.5, color='gray', linewidth=0.5)

    plt.tight_layout()
    plt.savefig(f"Outputs/{folderName}/Residuals/{plotName}.png", dpi=300)
    plt.close("all")


def plot_tradeoff_scatter(tradeoffBinaryMatrix, residualsBootstrapMatrix, metricsNames, folderName):
    """
    Create and save scatter plots of metric residuals where trade-offs are detected.

    Parameters:
        tradeoffBinaryMatrix (np.ndarray): Binary matrix (n_metrics x n_metrics) indicating trade-offs.
        residualsBootstrapMatrix (np.ndarray): Bootstrap residuals (n_samples x n_metrics).
        metricsNames (np.ndarray): Names of the metrics.
        folderName (str): Subfolder name for saving outputs inside "Outputs".
    """
    n_metrics = len(metricsNames)

    for i in range(n_metrics):
        for j in range(i):
            if tradeoffBinaryMatrix[i, j]:
            #if (i+1==81 and j+1==67) or (i+1==77 and j+1==67) or (i+1==67 and j+1==62) or \
            #   (i+1==105 and j+1==81) or (i+1==105 and j+1==77) or (i+1==105 and j+1==62):
                x = residualsBootstrapMatrix[:, i]
                y = residualsBootstrapMatrix[:, j]

                plt.figure(figsize=(6, 6))
                plt.scatter(x, y, alpha=0.5, s=10)
                plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
                plt.axvline(0, color='black', linestyle='--', linewidth=0.8)
                plt.xlabel(metricsNames[i])
                plt.ylabel(metricsNames[j])
                plt.title(f"Joint residual distribution: {metricsNames[i]} vs {metricsNames[j]}")
                plt.grid(True, linestyle=':', alpha=0.6)
                plt.tight_layout()
                plt.savefig(f"Outputs/{folderName}/Residuals/Tradeoff_{metricsNames[i]}_vs_{metricsNames[j]}.png", dpi=300)
                plt.close("all")
