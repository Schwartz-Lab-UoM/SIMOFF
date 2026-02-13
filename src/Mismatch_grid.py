import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import seaborn as sns

def mismatch(function_result, model, num_columns=100): #function_result: the result of optimize_objective_from_qualitative_criteria_agreement function
    """
    Plots heatmaps of agreement (binary 0/1) for the first and last `num_columns` iterations.

    Parameters:
        function_result : tuple
        Output of SIMOFF
        Assumes: function_result must be set to be equal to the successful SIMOFF run, e.g. function result = result['reaction_id1',
        'reaction_id2'], where an 100% accurate solution was found using reaction_id1 and reaction_id2 as input
        model: cobra model
        num_columns (int): Number of columns (iterations) to plot per subplot (default: 100).
    """
    
    agreement_df = function_result[3]
    id_to_name = {rxn.id: rxn.name for rxn in model.reactions}
    agreement_df.index = agreement_df.index.map(lambda rid: id_to_name.get(rid, rid))  # fallback to ID if name missing
    num_cols = agreement_df.shape[1]
    cols_to_plot = min(num_columns, num_cols)

    first_chunk = agreement_df.iloc[:, :cols_to_plot]
    last_chunk = agreement_df.iloc[:, -cols_to_plot:]

    # Define custom colormap
    cmap = mcolors.ListedColormap(["#d63384", "#20c997"])  # red for 0, green for 1
    bounds = [-0.5, 0.5, 1.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Create subplots
    fig_height = max(6, len(agreement_df) * 0.4)
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(cols_to_plot / 4, fig_height))

    # Plot first chunk
    sns.heatmap(
        first_chunk,
        ax=axes[0],
        cmap=cmap,
        norm=norm,
        cbar=False,
        linewidths=0.2,
        linecolor='white',
        square=False,
        xticklabels=max(1, cols_to_plot // 10),
        yticklabels=True
    )
    axes[0].set_title(f"First {cols_to_plot} Iterations", fontsize=10)
    axes[0].set_xlabel("Iteration")
    axes[0].set_ylabel("Reaction ID")

    # Plot last chunk
    sns.heatmap(
        last_chunk,
        ax=axes[1],
        cmap=cmap,
        norm=norm,
        cbar=False,
        linewidths=0.2,
        linecolor='white',
        square=False,
        xticklabels=max(1, cols_to_plot // 10),
        yticklabels=True
    )
    axes[1].set_title(f"Last {cols_to_plot} Iterations", fontsize=10)
    axes[1].set_xlabel("Iteration")
    axes[1].set_ylabel("Reaction ID")

    axes[1].set_xticklabels(
        axes[1].get_xticklabels(),
        rotation=45,
        ha='right',
        fontsize=10
    )

    plt.tight_layout()
    plt.show()
